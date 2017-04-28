## Tidy MRexperiment
get_taxa <- function(mrexp_list){
      mrexp_list %>% map_df(fData, .id = "pipe") %>% 
            dplyr::rename(feature_id = OTUname)
}