## Two sample titatration 16S study experimental design sample sheet

get_sampleSheet <- function(){
      ## biological sample ids
      biosample_id <- c("E01JH0004","E01JH0011","E01JH0016","E01JH0017","E01JH0038")
      labID = c(3,8,12,13,28,138,115,999,177,105)
      treatment = c(rep("pre",5),rep("post",5))
      timepoint = c(rep(20,5), c(4,2,2,5,2))
      biosample_df <- data_frame(biosample_id= rep(biosample_id,2), labID, treatment, timepoint)
      
      
      
      ## Experimental Design
      make_plan_a_96_df <- function(biosample_id, titration){
            ## experimental design for 96 sample layout
            bio_replicates <- rep(biosample_id,each = 2*length(titration))
            titration_96 <- rep(titration, times = 2*length(biosample_id))
            
            ntc <- data_frame(biosample_id = rep("NTC",6),
                              sample_type = "control", titration = NA)
            
            df_96 <- data_frame(biosample_id = bio_replicates, 
                                sample_type = "titration", 
                                titration = titration_96) %>% 
                  bind_rows(ntc) %>% 
                  mutate(sample_type = ifelse(titration %in% c(0,20), "unmixed",sample_type)) 
      }
      
      make_plan_a_df <- function(biosample_id, titration){
          ## PCRs
          plate_1 <- make_plan_a_96_df(biosample_id, titration) %>%
              mutate(pcr_16S_plate = 1, pcr_16S_id = 1:n())
          plate_2 <- make_plan_a_96_df(biosample_id, titration) %>%
              mutate(pcr_16S_plate = 2, pcr_16S_id = (n() + 1):(2*n()))
          pcr_plates <- bind_rows(plate_1, plate_2)
      
          ## Barcode
          barcode_jhu <- pcr_plates %>% mutate(barcode_lab = "JHU",
                                               barcode_id = 1:n())
          barcode_nist <- pcr_plates %>% mutate(barcode_lab = "NIST",
                                                barcode_id = (n() + 1):(2*n()))
          seq_plates <- bind_rows(barcode_jhu, barcode_nist)
      
          seq_jhu <- seq_plates %>% mutate(seq_lab = "JHU")
          seq_nist <- seq_plates %>% mutate(seq_lab = "NIST")
          return(list(pcr_sample_sheet = pcr_plates,
                      seq_sample_sheet = bind_rows(seq_jhu, seq_nist)))
      }
      
      titration <- c(20,0:4,5,10,15)
      plan_a_sample_sheets <- make_plan_a_df(biosample_id, titration)
      plan_a_pcr <- plan_a_sample_sheets$pcr_sample_sheet
      plan_a <- plan_a_sample_sheets$seq_sample_sheet
      
      ## Format PCR plate layout -----------------------------------------------------
      pcr_plate_layout<- plan_a_pcr %>%
          filter(pcr_16S_plate == 1, titration != 20 | is.na(titration)) %>%
          mutate(half = c(rep(c(rep(0,8),rep(6,8)),5),rep(c(0,6), each = 3)),
                 col = half + as.numeric(factor(biosample_id)),
                 row = c(rep(c("A","B","C","D","E","F","G","H"),10),
                         rep(c("A","D","H"), 2)))
      
      pcr_plate_layout <- plan_a_pcr %>% filter(pcr_16S_plate == 1, titration == 20) %>%
          mutate(half = rep(c(0,6), 5), col = rep(c(6,12), 5),
                 row = rep(c("F","B","C","E","G"),each = 2)) %>%
          bind_rows(pcr_plate_layout) %>% select(-pcr_16S_plate)
      
      pcr_plates <- pcr_plate_layout %>% mutate(pcr_16S_id = pcr_16S_id + 96) %>%
          bind_rows(pcr_plate_layout, .)
      
      
      ## Annotating with barcode information -----------------------------------------
      illumina_index <- read_csv("data/raw/illumina_index.csv", comment = "#") %>%
          dplyr::rename(kit_version = `Kit Version`, index_name = `i7 index name`)
      assign_plate <- c(A= 1, B= 1, C= 2, D = 2)
      assign_lab <- c(A= "JHU", B= "NIST", C= "NIST", D = "JHU")
      
      ## split and assign plate position
      forward_index <- illumina_index %>% filter(Index == "i5") %>%
          group_by(kit_version) %>% mutate(row = LETTERS[1:8]) %>%
          mutate(pcr_16S_plate = assign_plate[kit_version],
                 barcode_lab = assign_lab[kit_version]) %>%
          dplyr::rename(For_Index = Index, For_Index_ID = index_name,
                 For_sample_sheet = sample_sheet, For_barcode_seq = barcode_seq) %>%
          select(kit_version, For_Index_ID, pcr_16S_plate, barcode_lab, row)
      
      
      reverse_index <- illumina_index %>% filter(Index == "i7")%>%
          group_by(kit_version) %>% mutate(col = 1:12) %>%
          mutate(pcr_16S_plate = assign_plate[kit_version],
                 barcode_lab = assign_lab[kit_version]) %>%
          dplyr::rename(Rev_Index = Index, Rev_Index_ID = index_name,
                 Rev_sample_sheet = sample_sheet, Rev_barcode_seq = barcode_seq) %>%
          select(kit_version, Rev_Index_ID, pcr_16S_plate, barcode_lab, col)
      
      annotated_index <- full_join(forward_index, reverse_index)
      
      full_sample_sheet <- left_join(plan_a, pcr_plates) %>% left_join(annotated_index)
      
      ## Sample ID format ------------------------------------------------------------
      sam_id <- paste0("B", 0:5)
      names(sam_id) <- c("NTC","E01JH0004","E01JH0011","E01JH0016",
                         "E01JH0017","E01JH0038")
      
      mix_id <- paste0("M", 0:9)
      mix_id_df <- data_frame(mix_id, titration = c(NA, titration))
      pcr_id_df <- data_frame(pcr_id = paste0("P", 1:4),
                              pcr_16S_plate = rep(1:2, each = 2), half = rep(c(0,6), 2))
      lib_id <- c("JHU" = "L1","NIST" = "L2")
      seq_id <- c("JHU" = "S1","NIST" = "S2")
      
      sample_sheet <- full_sample_sheet %>%
          left_join(pcr_id_df) %>% left_join(mix_id_df) %>%
          mutate(sam_id = sam_id[biosample_id],
                 lib_id = lib_id[barcode_lab], seq_id = seq_id[seq_lab],
                 ID = paste(sam_id, mix_id, pcr_id, lib_id, seq_id, sep ="_")) %>%
          unite(pos, row, col,sep = "") %>%
          select(ID, biosample_id, titration, pcr_16S_plate, pos, barcode_lab,kit_version,
                 For_Index_ID, Rev_Index_ID, seq_lab)
      
      
      ## adding corrected titrations
      sample_sheet <- sample_sheet %>% 
            mutate(t_fctr = factor(titration),
                  t_fctr = fct_recode(t_fctr, `0` = "20", `20` = "0"))
      write_csv(sample_sheet, "data/raw/sample_sheet.csv")
      sample_sheet
}

sampleSheet <- get_sampleSheet()

ProjectTemplate::cache("sampleSheet")

rm(get_sampleSheet)

tstPcrSampleMetaData <- sampleSheet %>% filter(barcode_lab == "JHU",seq_lab == "JHU") %>% 
            select(-ID, -barcode_lab, -kit_version, -For_Index_ID, -Rev_Index_ID, -seq_lab)

ProjectTemplate::cache("tstPcrSampleMetaData")