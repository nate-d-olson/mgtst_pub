## Generating rds for phyloseq and MRexp objects
library(metagenomeSeq)
library(phyloseq)
library(mgtst)

proj_dir <- "~/Projects/16S_etec_mix_study"
pipe_dir <- file.path(proj_dir, "analysis","pipelines")
qiime_dir <- file.path(pipe_dir, "qiime")




## Munging mothur data files
mrexp_filenames <- list(qiime_denovo_chimerafilt = "../data/mrexp_qiime_denovo_chimera_filt.RDS",
                        qiime_denovo_nochimerafilt = "../data/mrexp_qiime_denovo_nochimera.RDS",
                        qiime_openref_chimerafilt = "../data/mrexp_qiime_refclus_chimera_filt.RDS",
                        qiime_openref_nochimerafilt = "../data/mrexp_qiime_refclus_nochimera.RDS")

mrexp_obj <- mrexp_filenames %>% map(readRDS)

#fvarLabels(mrexp_obj$qiime_openref_chimerafilt) <- paste0("taxonomy",1:7)

## Rename Taxonomy Levels
for(i in 1:4){
      fvarLabels(mrexp_obj[[i]]) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")    
}

mrexp_obj[[1]]@phenoData@data
