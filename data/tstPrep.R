## 1st stage PCR concentration after clean-up

## Two sample titration prep
get_tstPrep <- function(){
      tstPrep <- read_excel("data/raw/1st 16S_PCR protocol.xls",
                             sheet = 4, skip = 41)
      colnames(tstPrep) <- c("Titration","ERCC_ul","Prep_DNA",
                              "Prep_vol","Final_vol")
      
      write_csv(tstPrep,"data/raw/tstWetLabQC.csv")
}
tstPrep <- get_tstPrep()

ProjectTemplate::cache("tstPrep")

rm(get_tstPrep)