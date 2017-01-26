## MA plot functions

get_ma_df <- function(i, mrexp, present, depth){
      ## subset pre and post samples
      mrexp_pre_post <- mrexp[, which(pData(mrexp)$sampleID == i & 
                                            pData(mrexp)$dilution %in% c(-1,0))]
      
      mrexp_pre_post_filt <- filterData(mrexp_pre_post, present = present, depth = depth)
      
      ## Normalized count table
      count_tbl <- cumNormMat(mrexp_pre_post_filt, p = 0.75)
      
      ## Pre and Post sample ids
      pre_sams <- mrexp_pre_post_filt %>% pData() %>% {.[.$dilution == 0, ]} %>% rownames()
      post_sams <- mrexp_pre_post_filt %>% pData() %>% {.[.$dilution == -1, ]} %>% rownames()
      
      ## Mean normalized counts by OTU
      rowMeans_pre <- count_tbl %>% {.[,colnames(.) %in% pre_sams]} %>%  rowMeans()
      rowMeans_post <- count_tbl %>% {.[,colnames(.) %in% post_sams]} %>%  rowMeans()
      
      
      ## dataframe with MA values
      ## A - mean counts across all samples
      ## logFC (M) - log10 fold change pre/post
      data_frame(sampleID = i, 
                 otu = rownames(count_tbl), 
                 A = rowMeans(count_tbl),
                 pre_means = rowMeans_pre,
                 post_means = rowMeans_post) %>% 
            mutate(logFC = log2(pre_means/post_means))
      
}

get_ma_df_by_sample <- function(mrexp, sample_ids, present=4, depth=1){
      sample_ids %>% map_df(get_ma_df, mrexp,present, depth)
}

plot_by_sample_ma <- function(ma_df){
      ggplot(ma_df) + geom_point(aes(x = A, y = logFC, group = otu), alpha = 0.5) + 
            geom_hline(aes(yintercept = -1), linetype = 2) + 
            geom_hline(aes(yintercept = 1), linetype = 2) + 
            theme_bw() +
            scale_x_log10() +
            facet_wrap(~sampleID)
}