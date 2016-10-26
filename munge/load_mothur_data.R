## Load mothur
library(phyloseq)
library(mgtst)
library(dplyr)

mothur_dir <- file.path("~/Projects/16S_etec_mix_study/analysis/pipelines/mothur/data/process")
mothur_file_root <- "mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list"
mothur_list_filename <- paste0(mothur_file_root,".list")
mothur_share_filename <- paste0(mothur_file_root,".shared")
mothur_tax_filename <- paste0(mothur_file_root,".0.03.cons.taxonomy")

ps_mothur <- phyloseq::import_mothur(mothur_list_file = file.path(mothur_dir,mothur_list_filename),
                        mothur_shared_file = file.path(mothur_dir,mothur_share_filename),
                        mothur_constaxonomy_file = file.path(mothur_dir,mothur_tax_filename))

## add metadata
mothur_sam_names <- ps_mothur %>% sample_names()

meta_df <- sample_sheet %>%
      mutate(pos_ns = str_replace(pos, "_",""),
             sam_names = paste(pcr_16S_plate, pos_ns, sep = "-")) %>%
      filter(seq_lab == "JHU", 
             barcode_lab == "JHU", 
             sam_names %in% mothur_sam_names) %>% 
      as.data.frame()

rownames(meta_df) <- meta_df$sam_names
sample_data(ps_mothur) <- meta_df[match(mothur_sam_names,meta_df$sam_names),]

## Save RDS
saveRDS(ps_mothur, "../data/phyloseq_mothur.RDS")

## Convert to MRexperiment and save RDS
mrexp_mothur <- phyloseq_to_metagenomeSeq(ps_mothur)
saveRDS(mrexp_mothur, "../data/mrexp_mothur.RDS")

########## NOTES ########### 

# Mothur taxa table has 8 ranks, based on cons.tax.summary file it looks like
# the 8th rank is a subspecies, this is may vary by org...

