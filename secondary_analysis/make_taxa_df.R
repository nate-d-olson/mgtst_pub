## Creating taxa_df for set and relative abundance plots
library(phyloseq)
library(tidyverse)

pipeline_dir <- "~/Projects/mgtst_pipelines"
ps_list <- list(
      dada2 = file.path(pipeline_dir, "dada2/dada_ps.rds"),
      mothur =  file.path(pipeline_dir, "mothur/mothur_ps.rds"),
      qiime =  file.path(pipeline_dir, "qiime/qiime_ps.rds")
) %>% 
      map(readRDS)

get_taxa_df <- function(ps){
      taxa_df <- tax_table(ps) %>% 
            as.data.frame() %>%
            tibble::rownames_to_column(var = "feature_id") %>% 
            mutate_all(funs(str_remove(.,".__")))
      
      taxa_df <- tibble(feature_id = taxa_names(ps), 
                        f_counts = taxa_sums(ps)) %>% 
            left_join(taxa_df) 
      
      taxa_df %>% 
            group_by(Rank1, Rank2, Rank3, Rank4, Rank5, Rank6) %>% 
            summarise(total_count = sum(f_counts))
}

taxa_df <- ps_list %>% map_dfr(get_taxa_df, .id = "pipe")

taxa_df <- taxa_df %>% 
      group_by(pipe) %>% 
      mutate(rel_abu = total_count / sum(total_count))

saveRDS(taxa_df, "data/taxa_df.RDS")