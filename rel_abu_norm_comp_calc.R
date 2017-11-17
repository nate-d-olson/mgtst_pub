library(metagenomeSeq)
library(tidyverse)

mrexp_obj <- readRDS("~/Projects/mgtst_pipelines/mothur/mothur_mrexp.rds")

mrexp_to_edgeR <- function(mrexp_obj, method){
    require(edgeR)
    ## Extracting count data - no scaling or transformation
    x <- mrexp_obj %>% 
        metagenomeSeq::MRcounts(norm = FALSE, log = FALSE, sl = 1) %>% 
        as.matrix()
    x <- x + 1 # add 1 to prevent log(0) issues
    
    ## Use taxonomy information at gene annotations 
    ## - Where OTUname is incorporated into the results
    taxonomy <- fData(mrexp_obj)
    if (!is.null(taxonomy)) {
        taxonomy <- taxonomy %>% as.matrix() %>% data.frame()
    }
    
    
    ## Convert into a DGEList
    y <- DGEList(counts = x, genes = taxonomy,
                 remove.zeros = TRUE)
    
    
    ## Calc normalization factors
    if (method == "CSS") {
        z <- edgeR::calcNormFactors(y, method = "none")
        mg_nf <- metagenomeSeq::calcNormFactors(mrexp_obj, p = 0.75)
        z$samples$norm.factors <- mg_nf$normFactors[which(rownames(mg_nf) == rownames(z$samples))]
    } else if (method == "TSS") {
        z <- edgeR::calcNormFactors(y, method = "none")
        z$samples$norm.factors <- z$samples$lib.size
    } else {
        z <- edgeR::calcNormFactors(y, method = method)
    }
    
    ## Check for division by zero inside `calcNormFactors`
    if ( !all(is.finite(z$samples$norm.factors))) {
        stop("Something wrong with edgeR::calcNormFactors on this data,
             non-finite $norm.factors, consider changing `method` argument.")
    }
    z 
}

get_norm_count_df <- function(method, mrexp){
    ## extract normalizaed counts
    count_mat <- mrexp_to_edgeR(mrexp, method = method) %>% 
        cpm(normalized.lib.sizes = TRUE)

    ## Scalling CSS and TSS so they are on the same scale as other normalization method. 
    if (method == "CSS") {
        count_mat = count_mat * 1000
    } else if (method == "TSS") {
        count_mat = count_mat * 100000
    } 

    ## Tidy counts
    count_df <- count_mat %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "feature_id") %>% 
        gather("id", "count", -feature_id)

    ## Adding sample data
    sample_dat <- pData(mrexp) %>% 
        select(biosample_id, id, t_fctr) 
    
    ## Count data and sample data 
    count_df <- left_join(count_df, sample_dat)

    ## Calculating replicate mean and variance 
    count_df %>% 
        group_by(biosample_id, feature_id, t_fctr) %>% 
        summarise(mean_count = mean(count),
                  var_count = var(count))
}

norm_count_df <- list(RAW = "none", 
                      RLE = "RLE", TMM = "TMM",
                      UQ = "upperquartile", 
                      CSS = "CSS", TSS = "TSS") %>% 
    map(get_norm_count_df, mrexp = mrexp_obj) %>% 
    map_df(bind_rows,.id = "norm_method") %>% 
    mutate(pipe = "mothur")


saveRDS(norm_count_df, "~/Desktop/norm_count_df.RDS")
