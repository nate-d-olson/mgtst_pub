## Generating rds for phyloseq and MRexp objects
library(metagenomeSeq)
library(phyloseq)
library(mgtst)

proj_dir <- "~/Projects/16S_etec_mix_study"
pipe_dir <- paste0(proj_dir, "/analysis/pipelines")
qiime_dir <- paste0(pipe_dir, "/qiime/")
dada2_dir <- paste0(pipe_dir, "/dada2/")
mothur_dir <- paste0(pipe_dir, "/mothur/")

save_ps_and_mrexp <- function(biom_file, out_dir, pipeline){
    ### Phyloseq
    ps <- load_biom_jhu(paste0(out_dir,biom_file), pipeline)
    saveRDS(ps, paste0(out_dir, "phyloseq_obj.rds"))

    ### MetagenomeSeq
    mrexp <- phyloseq_to_metagenomeSeq(ps)
    saveRDS(mrexp, paste0(out_dir,"mrexp_obj.rds"))
}

## QIIME -----------------------------------------------------------------------
###     Standard Pipeline
biom_file <- "otu_table_mc2_w_tax_no_pynast_failures.biom"
out_dir <- paste0(qiime_dir,"/otus_uc_fast/")
save_ps_and_mrexp(biom_file, out_dir, pipeline = "qiime")

###     Standard Pipeline with chimera filter
## Error with reading biom
#biom_file <- "otu_table_mc2_w_tax_no_pynast_failures.biom"
#out_dir <- paste0(qiime_dir,"/otus_uc_fast_no_chimera/")
#save_ps_and_mrexp(biom_file, out_dir, pipeline = "qiime")

###      De Novo Clustering
biom_file <- "otu_table.biom"
out_dir <- paste0(qiime_dir,"/otus_uc_fast_denovo/")
save_ps_and_mrexp(biom_file, out_dir, pipeline = "qiime")

###      De Novo Clustering with chimera filter
# biom_file <- "otu_table.biom"
# out_dir <- paste0(qiime_dir,"/otus_uc_fast_denovo_no_chimera/")
# save_ps_and_mrexp(biom_file, out_dir, pipeline = "qiime")


## DADA2 -----------------------------------------------------------------------
ps <- readRDS(paste0(dada2_dir, "/processed_data/ps_obj-2016-04-25.rds"))
### issue with taxa Seq3692
ps <- prune_taxa(!is.na(taxa_sums(ps)),ps)
ps <- prune_samples(sample_sums(ps)> 0, ps)
mrexp <- phyloseq_to_metagenomeSeq(ps)
saveRDS(mrexp, paste0(dada2_dir, "/processed_data/mrexp_obj-2016-04-25.rds"))


## Mothur ----------------------------------------------------------------------
out_dir <- paste0(mothur_dir, "/data/process/")
biom_file <- "mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.biom"
save_ps_and_mrexp(biom_file, out_dir, pipeline = "mothur")
