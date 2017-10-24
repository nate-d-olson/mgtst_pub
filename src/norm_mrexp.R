## Normalize MRexperiment Matricies and return data frame

# raw counts - no normalization or transformation
calc_raw_counts <- function(mrexp){
      mrexp@assayData$counts %>% as_tibble() %>% 
            rownames_to_column(var = "otuID") %>% 
            gather("samID","count",-otuID) %>% 
            left_join(sam_dat) 
}

# CSS normalization, default upper quartile
calc_css_counts <- function(mrexp, norm = TRUE,log = TRUE,sl = 1000, p = 0.75){
      mrexp %>% cumNorm(p = p) %>% 
            MRcounts(norm, log, sl) %>% as_tibble() %>% 
            rownames_to_column(var = "otuID") %>% 
            gather("samID","count",-otuID) %>% 
            left_join(sam_dat) 
}

# TSS - proportion scaling
# TSS from http://mixomics.org/mixmc/normalisation/ 
calc_tss_counts <- function(mrexp){
      mrexp@assayData$counts %>% {apply(., 2, function(x){ x/sum(x) })} %>% 
            as_tibble() %>% rownames_to_column(var = "otuID") %>%
            gather("samID","count",-otuID) %>% 
            left_join(sam_dat) 
}

# TSS with log2 transformation
calc_tsslog_counts <- function(mrexp){
      mrexp@assayData$counts %>% 
      {apply(., 2, function(x){ x/sum(x) })} %>% {log2(. + 1)} %>%
            as_tibble() %>% rownames_to_column(var = "otuID") %>%
            gather("samID","count",-otuID) %>% 
            left_join(sam_dat) 
}

## DESeq method median of ratios
calc_dsq_counts <- function(mrexp){
      mrexp@assayData$counts %>% {./estimateSizeFactorsForMatrix(.)} %>% 
            as_tibble() %>% rownames_to_column(var = "otuID") %>%
            gather("samID","count",-otuID) %>% 
            left_join(sam_dat) 
}

## DESeq with log2 transformation
calc_dsqlog_counts <- function(mrexp){
      mrexp@assayData$counts %>% 
      {./estimateSizeFactorsForMatrix(.)} %>% {log2(. + 1)} %>%
            as_tibble() %>% rownames_to_column(var = "otuID") %>%
            gather("samID","count",-otuID) %>% 
            left_join(sam_dat) 
}
