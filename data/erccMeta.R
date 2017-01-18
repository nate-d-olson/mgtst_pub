# ERCC Spike-in ----------------------------------------------------------------

## loading source data files
get_erccMeta <- function(){
      ercc_assay <- read_tsv("data/raw/thermo_fisher_ercc_qPCR_assays.txt")
      ercc_spike <- read_tsv("data/raw/ERCC-Spike-ins.tsv",comment = "#") %>%
          rename(`Gene Symbol` = Control)
      
      
      ## Assays selected to minimize differences in amplicon length between the pre
      ## and post qPCR assays.
      selected_assays <- c("Ac03459877_a1", "Ac03459922_a1", "Ac03459958_a1",
                           "Ac03459987_a1", "Ac03460000_a1","Ac03460028_a1",
                           "Ac03460000_a1","Ac03460028_a1","Ac03459872_a1",
                           "Ac03460039_a1","Ac03459892_a1","Ac03459925_a1")
      
      ercc_spike_ins <- ercc_assay %>%
          select(`Gene Symbol`, `Assay ID`,`Amplicon Length`) %>%
          filter(`Assay ID` %in% selected_assays) %>% full_join(ercc_spike, .) %>%
          rename(ercc_id = `Gene Symbol`,
                 biosample_id = Sample_spike,
                 treatment = Treatment,
                 assay_id = `Assay ID`,
                 GC = `GC (excluding polyA tail)`,
                 amplicon_length = `Amplicon Length`)
      
      write_csv(ercc_spike_ins,"data/raw/ercc_spike_ins.csv")
      ercc_spike_ins
}
erccMeta <- get_erccMeta()

ProjectTemplate::cache("erccMeta")

rm(get_erccMeta)
