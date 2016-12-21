## Extract ERCC qPCR data

## Standard Curve 2016-09-19 ---------------------------------------------------
get_ercc <- function(xcl_file){
      cols <- c("ercc_plate","ercc_std",
                paste0("std_rep",rep(1:3, each = 2),"_",rep(c("Ct","quant"),3)),
                "sampleID",
                paste0("sam_rep",rep(1:3, each = 2),"_",rep(c("Ct","quant"),3)))
      
      ## ERCC plasmid id
      ercc_assays <- rep(c(84,12,34,157,57,108,130,2,92,35), each = 8)
      
      read_excel(path = xcl_file, sheet = "ERCC_Quant_20161208",
                 skip = 3, col_names = FALSE, na = "Undetermined") %>%
            set_colnames(cols) %>% filter(!is.na(ercc_plate)) %>% 
            add_column(ercc = ercc_assays) %>% 
            gather("key","value",-ercc_plate,-ercc,-ercc_std, -sampleID) %>% 
            mutate(sampleID = if_else(grepl("std",key), ercc_std, sampleID)) %>% 
            select(-ercc_std) %>% separate(key, c("sample_type","rep","key"), sep = "_") %>% 
            spread(key,value) %>% 
            mutate(Ct = as.numeric(Ct), quant = as.numeric(quant))
}


qpcrERCC <- get_ercc("data/raw/MixStudy_Nate_20161209.xls") 
      
ProjectTemplate::cache("qpcrERCC")