## Biological Sample Info
get_biosampleInfo <- function(){
      ## sample selection data frame
      biosample_ids <- c("E01JH0004","E01JH0011","E01JH0016","E01JH0017","E01JH0038")
      lab_id = c(3,8,12,13,28,138,115,999,177,105)
      treatment = c(rep("Pre",5),rep("Post",5))
      timepoint = c(rep(-1,5), c(4,2,2,5,2))
      biosample_df <- data_frame(biosample_id= rep(biosample_ids,2), lab_id, treatment, timepoint)
      
      sampling_info <- read_excel("data/raw/1st 16S_PCR protocol.xls", skip = 1)
      colnames(sampling_info) <- c("biosample_id","Initials",
                                   colnames(sampling_info)[-c(1,2)])
      sampling_info <- sampling_info[sampling_info$biosample_id %in% biosample_ids,]
      
      
      subject_df <- sampling_info[,1:2]
      
      
      subject_info <- data_frame()
      j <- 3 # column index
      j_colnames <- c("stool_weight","processing_date","lab_id")
      for(i in c(-1:7,9,28,84)){
          ## day specific df
          df <- sampling_info[,j:(j+2)]
          colnames(df) <- j_colnames
          df <- df %>% mutate(stool_weight = as.character(stool_weight))
          j <- j + 3
          subject_info <- subject_df %>% mutate(timepoint = i) %>%
              bind_cols(df) %>% bind_rows(subject_info)
      }
      
      biosampleInfo <- subject_info %>%
          mutate(lab_id = as.numeric(lab_id),
                 stool_weight = as.numeric(stool_weight)) %>%
          select(-Initials) %>%
          left_join(biosample_df, .)
      
      ### Sample DNA concentration
      DNA_conc <- read_excel("data/raw/1st 16S_PCR protocol.xls",sheet = 2)
      DNA_conc <- DNA_conc[,1:9]
      colnames(DNA_conc) <- c("treatment","Tube","SampleID","lab_id","conc_ngul",
                              "vol_ul","total_ug","sample_ul","tris_ul")
      sample_dna_con <- DNA_conc %>%
          separate(SampleID, c("BioRep","biosample_id","DayCol","timepoint")) %>%
          select(-DayCol, -BioRep, -Tube) %>%
          mutate(timepoint = ifelse(timepoint == "1","-1", timepoint) %>% as.numeric(),
                 treatment = ifelse(timepoint == -1,"Pre","Post"))
      
      biosampleInfo <- left_join(biosampleInfo, sample_dna_con)
      
      
      write_csv(subject_info,"data/raw/biosample_info.csv")
      biosampleInfo
}

biosampleInfo <- get_biosampleInfo()

ProjectTemplate::cache("biosampleInfo")

rm(get_biosampleInfo)