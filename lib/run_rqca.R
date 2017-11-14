library(stringr)
library(purrr)

## Fastq file list -------------------------------------------------------------

project_dir <- "~/Projects/16S_etec_mix_study"
#fastq_dir <- file.path(project_dir, "data/raw_seq/nist_run1/170209_16s_MiSeq_Run_Folder")
fastq_dir <- file.path(project_dir, "data/raw_seq/nist_run2")
# fastq_dir <- file.path(project_dir, "data/raw_seq/nist_run1")
dat_files <- list.files(path = fastq_dir,pattern = "NIST.*001.fastq.gz$",
                        recursive = TRUE, full.names = TRUE)

## Seq dataset IDs
seq_ds_id <- dat_files %>% str_split("/") %>% flatten_chr() %>%
    grep(pattern = "fastq.gz", .,value = TRUE) %>%
    str_replace(".fastq.gz","") %>% paste0("Fq_",.) %>%
    str_replace("-",".")
names(dat_files) <- seq_ds_id

## Read groups to pass as argument to RqcA
read_groups <- rep(NA,length(dat_files))
read_groups[grepl("/NIST-1-.*_R1", dat_files)] <- "plate1_F"
read_groups[grepl("/NIST-1-.*_R2", dat_files)] <- "plate1_R"
read_groups[grepl("/NIST-2-.*_R1", dat_files)] <- "plate2_F"
read_groups[grepl("/NIST-2-.*_R2", dat_files)] <- "plate2_R"


## Running RqcA ----------------------------------------------------------------
library(Rqc)
qa_list <- list()
step_size <- 1 #for all data

n_files <- length(dat_files)
for(i in 0:((n_files/step_size) - 1)){
    print(i)
    qa_list <- c(qa_list,
                 rqcQA(dat_files[(step_size * i+1):(step_size*(i + 1))],
                       group = read_groups[(step_size *i+1):(step_size *(i+1))],
                       workers = step_size))
}

names(qa_list) <- names(dat_files)

# setwd("~/Projects/16S_etec_mix_study/analysis/nist_run1_qa/")
saveRDS(qa_list,"rqcQA_list.rds")



