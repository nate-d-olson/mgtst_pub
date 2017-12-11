require(tidyverse)
require(metagenomeSeq)
mrexp_obj <- readRDS("~/Projects/mgtst_pipelines/mothur/mothur_mrexp.rds")

pdat <- pData(mrexp_obj)  %>% 
      mutate(t_fctr = factor(t_fctr, level = c(0:5,10,15,20)))
t1_pdat <- pdat %>% dplyr::select(biosample_id, id, t_fctr) %>% 
      dplyr::rename(T1 = t_fctr, T1_id = id) 
t2_pdat <- pdat %>% dplyr::select(biosample_id, id, t_fctr) %>% 
      dplyr::rename(T2 = t_fctr, T2_id = id) 
titration_comp_df <- left_join(t1_pdat, t2_pdat) %>% 
      filter(as.numeric(T1) < as.numeric(T2)) %>% 
      group_by(biosample_id, T1, T2) %>% 
      summarise(sam_names = c(T1_id,T2_id) %>% unique() %>% list(),
                n_sams = c(T1_id,T2_id) %>% unique() %>% length())

make_titration_comp_subset_df <- function(mrexp_obj, titration_comp_df){
      
      ## Subsetting mrexp for titrations
      subset_mrexp <- function(sam_names){
            mrexp_obj %>% {.[,which(colnames(.) %in% sam_names)]}
      }
      
      ## Dataframe with list of subsetted mrexp
      titration_comp_df %>% 
            mutate(mrexp_sub = map(sam_names, subset_mrexp))
}

titration_comp_mothur_df <- mrexp_obj %>%
      make_titration_comp_subset_df(titration_comp_df)

mrexp_to_edgeR <- function(mrexp_obj, group, method = "RLE", ...){
      require(edgeR)
      ## Extracting count data - no scaling or transformation
      x <- mrexp_obj %>% 
            metagenomeSeq::MRcounts(norm = FALSE, log = FALSE, sl = 1) %>% 
            as.matrix()
      x <- x + 1 # add 1 to prevent log(0) issues
      
      ## Check `group` argument
      if (identical(all.equal(length(group), 1), TRUE) & ncol(mrexp_obj) > 1) {
            ## Assumes grouop is a categorical sample variable name
            group <- pData(mrexp_obj) %>% .[,group]
            group <- as.numeric(as.character(group))
            T1 <- min(group)
            T2 <- max(group)
            group <- factor(group, levels = c(as.character(T1), as.character(T2)))
      }
      
      ## Use taxonomy information at gene annotations 
      ## - Where OTUname is incorporated into the results
      taxonomy <- fData(mrexp_obj)
      if (!is.null(taxonomy)) {
            taxonomy <- taxonomy %>% as.matrix() %>% data.frame()
      }
      
      
      ## Convert into a DGEList
      y <- DGEList(counts = x, group = group, genes = taxonomy,
                   remove.zeros = TRUE)
      
      
      ## Calc normalization factors
      if (method == "CSS") {
            z <- edgeR::calcNormFactors(y, method = "none")
            mg_nf <- metagenomeSeq::calcNormFactors(mrexp_obj, p = 0.75)
    
            z$samples$norm.factors <- mg_nf$normFactors[which(rownames(mg_nf) == rownames(z$samples))]
            
            ## EdgeR multiplies the library size * norm factor - diving by
            ## library size accounts for this
            z$samples$norm.factors <- z$samples$norm.factors/z$samples$lib.size
      } else if (method == "TSS") {
            z <- edgeR::calcNormFactors(y, method = "none")
            ## EdgeR multiplies the library size * norm factor - using the
            ## squared library size accounts for this
            z$samples$norm.factors <- z$samples$lib.size^2
      } else {
            z <- edgeR::calcNormFactors(y, method = method)
      }
      
      
      ## Check for division by zero inside `calcNormFactors`
      if ( !all(is.finite(z$samples$norm.factors))) {
            stop("Something wrong with edgeR::calcNormFactors on this data,
                 non-finite $norm.factors, consider changing `method` argument.")
      }
      
      ## Estimate dispersions
      z %>% estimateCommonDisp() %>% 
          estimateTagwiseDisp()
}


calc_logfc_norm_comp <- function(comp_df, method) {
  logFC_edgeR_mothur_df <- comp_df %>% 
        mutate(fit = map(mrexp_sub, mrexp_to_edgeR, group = "t_fctr", method = method),
               fit = map(fit, exactTest)) 
  
  logFC_edgeR_mothur_coefs_df <- logFC_edgeR_mothur_df %>% 
        mutate(fit_coefs = map(fit, topTags, n = Inf, adjust.method = "BH")) 
  
  logFC_edgeR_mothur_coefs_df %>% 
        select(biosample_id, T1, T2, fit_coefs) %>% 
        mutate(fit_coefs = map(fit_coefs, ~.@.Data[[1]])) %>% 
        unnest()
  
}

# raw_logFC <- calc_logfc_norm_comp(titration_comp_mothur_df, method = "none")
# rle_logFC <- calc_logfc_norm_comp(titration_comp_mothur_df, method = "RLE")
# tmm_logFC <- calc_logfc_norm_comp(titration_comp_mothur_df, method = "TMM")
# uq_logFC <- calc_logfc_norm_comp(titration_comp_mothur_df, method = "upperquartile")
# css_logFC <- calc_logfc_norm_comp(titration_comp_mothur_df, method = "CSS")
# tss_logFC <- calc_logfc_norm_comp(titration_comp_mothur_df, method = "TSS")
# 
# norm_logFC_df <- list(RAW = raw_logFC, RLE = rle_logFC, TMM = tmm_logFC, 
#      UQ = uq_logFC, CSS = css_logFC, TSS = tss_logFC) %>% 
#       map_df(bind_rows,.id = "norm") 
# 
# saveRDS(norm_logFC_df, "~/Desktop/norm_logFC_df.RDS")
