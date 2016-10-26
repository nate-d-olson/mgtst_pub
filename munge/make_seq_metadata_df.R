## Raw seq stats and metadata
library(Rqc)
project_dir <- "~/Projects/16S_etec_mix_study/"

fastq_dir <- file.path(project_dir, "data/fastq/jhu_run2/")
dat_files <- list.files(path = fastq_dir,pattern = "001.fastq.gz$",
                        recursive = TRUE, full.names = TRUE)

seq_ds_id <- dat_files %>% str_split("/") %>% flatten_chr() %>% 
      grep(pattern = "fastq.gz", .,value = TRUE) %>%
      str_replace(".fastq.gz","") %>% paste0("Fq_",.) %>% 
      str_replace("-",".")

names(dat_files) <- seq_ds_id


read_groups <- rep(NA,length(dat_files))
read_groups[grepl("/1-.*_R1", dat_files)] <- "plate1_R1"
read_groups[grepl("/1-.*_R2", dat_files)] <- "plate1_R2"
read_groups[grepl("/2-.*_R1", dat_files)] <- "plate2_R1"
read_groups[grepl("/2-.*_R2", dat_files)] <- "plate2_R2"


## Run 
read_groups <- rep(NA,length(dat_files))
read_groups[grepl("/1-.*_R1", dat_files)] <- "plate1_R1"
read_groups[grepl("/1-.*_R2", dat_files)] <- "plate1_R2"
read_groups[grepl("/2-.*_R1", dat_files)] <- "plate2_R1"
read_groups[grepl("/2-.*_R2", dat_files)] <- "plate2_R2"
grp_df <- data_frame(read_group = read_groups, 
                     seq_ds_id, 
                     filename = basename(dat_files)) %>%
      separate(read_group,c("plate","Read")) %>% 
      mutate(ill_id = str_replace(filename, "_.*",""))

## read count data 
qa_list <- readRDS("../data/rqcQA_list.rds")
qa_file_info <- perFileInformation(qa_list) %>% 
      select(-format,-path)

## study metadata
data(sample_sheet)
meta_df <- sample_sheet %>% 
      mutate(pos_ns = str_replace(pos, "_",""),
             ill_id = paste(pcr_16S_plate, pos_ns, sep = "-")) %>% 
      filter(seq_lab == "JHU", barcode_lab == "JHU") %>% 
      mutate(pcr_16S_plate = as.character(pcr_16S_plate)) %>% 
      left_join(grp_df) %>% left_join(qa_file_info)

saveRDS(meta_df,"../data/seq_metadata_df.RDS")