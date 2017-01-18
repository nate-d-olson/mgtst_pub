## 1st stage PCR concentration after clean-up

## Two sample titration prep
get_tstWetLabQC <- function(){
      tstWetLabQC <- read_excel("data/raw/1st 16S_PCR protocol.xls",
                             sheet = 4, skip = 41)
      colnames(tstWetLabQC) <- c("Titration","ERCC_ul","Prep_DNA",
                              "Prep_vol","Final_vol")
      
      write_csv(tstWetLabQC,"data/raw/tstWetLabQC.csv")
}
tstWetLabQC <- get_tstWetLabQC()

ProjectTemplate::cache("tstWetLabQC")

rm(get_tstWetLabQC)