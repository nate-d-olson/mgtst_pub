## Extract bacterial qPCR standard curves

## Standard Curve 2016-09-19 ---------------------------------------------------
get_std_0919 <- function(xcl_file){
      cols <- c("well","sample_name","conc","plate1","plate2","plate3")
      
      read_excel(path = xcl_file, sheet = "QDNA_20160919", 
                 skip = 3, col_names = FALSE, na = "Undetermined") %>% 
            select(-X12, -X13, -X5,-X7,-X9) %>% set_colnames(cols) %>% 
            filter(sample_name %in% paste0("Std",1:7)) %>% 
            gather("plate","Ct",-well, -sample_name, -conc) %>% 
            mutate(conc = as.numeric(conc), Ct = as.numeric(Ct), 
                   std = "zymo", date = "2016-09-19")
}

## Standard Curve 2016-09-19 ---------------------------------------------------
get_std_1209 <- function(xcl_file){
      cols <- c("well","sample_name","conc","plate1","plate2","plate3")
      full_cols <- c(paste("shan",cols,sep = "_"), paste("zymo",cols,sep = "_"))
      
      bac_std <- read_excel(path = xcl_file,
                            sheet = "ReDo_QDNA_20161209",skip = 3, 
                            na = "Undetermined", col_names = FALSE) %>% 
            select(-X5,-X7,-X9,-X15,-X17,-X19) %>% 
            filter(X2 %in% paste0("Std",1:7)) %>% set_colnames(full_cols)
      
      shan_std <- bac_std %>% select(starts_with("shan")) %>% 
            set_colnames(cols) %>% mutate(std = "shan")
      
      bac_std <- bac_std %>% select(starts_with("zymo")) %>% 
            set_colnames(cols) %>% mutate(std = "zymo") %>% 
            bind_rows(shan_std)
      
      bac_std %>% gather("plate","Ct",-well, -sample_name, -conc, -std) %>% 
            mutate(conc = as.numeric(conc), Ct = as.numeric(Ct), date = "2016-12-09")
}

## Generating full dataset and caching
qpcrBacStd <- "data/raw/MixStudy_Nate_20161209.xls" %>% 
      {full_join(get_std_0919(.), get_std_1209(.))}

ProjectTemplate::cache("qpcrBacStd")