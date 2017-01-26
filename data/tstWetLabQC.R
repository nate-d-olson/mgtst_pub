## DNA quant picogreen - first step 16S PCR
get_tstWetLabQC <- function(){
     pcr_16S_pico <- read_excel("data/raw/160224_pico_firstPCR.xlsx",
           col_names = c("pcr_16S_plate","pos","excite_emiss",
                         "conc_ngml", "conc_ngul","qubit"), skip = 9) %>% 
           filter(!is.na(pos), pos != "Well")

      write_csv(pcr_16S_pico,"data/raw/tstWetLabQC.csv")
      pcr_16S_pico
} 

tstWetLabQC <- get_tstWetLabQC()


ProjectTemplate::cache("tstWetLabQC") 

rm(get_tstWetLabQC)
