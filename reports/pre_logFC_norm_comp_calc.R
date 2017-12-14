library(tidyverse)
library(stringr)
## Only looking at pre specific and pre dominant features

## Normalized counts summarized across PCR replicates
norm_count_df <- readRDS("~/Desktop/norm_count_df.RDS")



## Presence absence info
pa_summary_df <- readRDS("~/Desktop/to_file/mgtst_RDS/pa_summary_anno_df.RDS") %>% 
    select(pipe, biosample_id, feature_id, T00, T20, pa_mixed)

## Pre - Post logFC - for features observed in all four pre PCR replicates
prepost_df <- norm_count_df %>% left_join(pa_summary_df) %>% 
    filter(T20 == 4, t_fctr %in% c(0,20)) %>% 
    select(-var_count, -sd_count, -T00, -T20, -pa_mixed) 

prepost_logFC <- prepost_df %>% 
    mutate(t_fctr = paste0("T",str_pad(t_fctr, 2, "left",0))) %>% 
    spread(t_fctr, mean_count) %>% 
    mutate(logFC = log2(T20/T00))

saveRDS(prepost_logFC,"~/Desktop/norm_logFC_prepost.RDS")

## Find features with logFC > 5 or undefined....
pre_logFC <- prepost_logFC %>% filter(logFC > 5)

### Similar numbers of features across normalization method except for UQ
## pre_logFC %>% group_by(norm_method, biosample_id) %>% summarise(count = n()) %>% spread(norm_method, count)


## Calculate logFC for titration comparisons ###################################

##### Titration comparisons 
meta_dat <- norm_count_df %>% ungroup() %>% 
    select(biosample_id, t_fctr) %>% 
    distinct() %>% 
    mutate(t_fctr = factor(t_fctr, level = c(0:5,10,15,20)))

t1_dat <- rename(meta_dat, T1 = t_fctr) 
t2_dat <- rename(meta_dat, T2 = t_fctr) 

## Generate titration comparison sets and reduce redundancy
titration_comp_df <- left_join(t1_dat, t2_dat) %>% 
    filter(as.numeric(T1) < as.numeric(T2))

##### Adding count data for pre specific and pre-dominant features
pre_features <- pre_logFC %>% select(norm_method, biosample_id, feature_id)
pre_norm_count_df <- left_join(pre_features, norm_count_df)

## Selecting columns needed to calculate logFC titration comparison 
pre_count <- pre_norm_count_df %>% 
    select(norm_method, biosample_id, feature_id, t_fctr, mean_count)

titration_comp_count_df <- titration_comp_df %>% 
    ## Titaration 1 counts 
    rename(t_fctr = T1) %>% left_join(pre_count) %>% 
    rename(T1_count = mean_count, T1 = t_fctr) %>% 
    ## Titration 2 counts 
    rename(t_fctr = T2) %>% left_join(pre_count) %>% 
    rename(T2_count = mean_count, T2 = t_fctr)
    

### Calculating logFC
titration_comp_logFC <- titration_comp_count_df %>% 
    mutate(logFC = log2(T2_count / T1_count))

## Saving logFC 
saveRDS(titration_comp_logFC, "~/Desktop/norm_logFC_df.RDS")
