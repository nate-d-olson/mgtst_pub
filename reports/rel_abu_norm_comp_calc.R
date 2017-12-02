library(metagenomeSeq)
library(tidyverse)

mrexp_obj <- readRDS("~/Projects/mgtst_pipelines/mothur/mothur_mrexp.rds")

get_norm_count_df <- function(method, mrexp){
    require(matrixStats)
    ## Extract count matrix from MRexperiment
    count_mat <- mrexp@assayData$counts 
    
    ## extract normalizaed counts
    if (method == "RAW") {
        ## Raw counts no normalization applied
        norm_factors = 1
    } else if (method == "UQ") {
        ## Upper quartile normalization
        count_mat[count_mat == 0] <- NA
        norm_factors <- colQuantiles(count_mat, p = 0.75 ,na.rm = TRUE)
    } else if (method == "CSS") {
        ## Cumulative sum scaling Paulson et al. 2013
        norm_factors <- metagenomeSeq::calcNormFactors(count_mat, p = 0.75)
        norm_factors <- norm_factors$normFactors
    } else if ( method == "TSS") {
        ## Total sum scaling 
        norm_factors <- colSums(count_mat)
    } else if (method %in% c("RLE","TMM")) {
        ## EdgeR RLE and TMM normalization methods
        norm_factors <- edgeR::calcNormFactors(count_mat, method = method) 
        norm_factors <- norm_factors * colSums(count_mat)
    } else {
        warning("Normalization method not defined")
    }
    
    ## Normalizing counts
    norm_mat <- sweep(count_mat, 2, norm_factors,'/')

    ## Tidy counts
    count_df <- norm_mat %>% 
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

norm_count_df <- list(RAW = "RAW", RLE = "RLE", TMM = "TMM",
                      UQ = "UQ", CSS = "CSS", TSS = "TSS") %>% 
    map(get_norm_count_df, mrexp = mrexp_obj) %>% 
    map_df(bind_rows, .id = "norm_method") %>%
    mutate(pipe = "mothur")


saveRDS(norm_count_df, "~/Desktop/norm_count_df.RDS")
