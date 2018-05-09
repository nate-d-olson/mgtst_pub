## Loading MRexperiment objects ------------------------------------------------
## argument - path to mgtst_pipelines
## returns a list pipelines objects
get_mrexp <- function(pipeline_dir){
      mrexp_files <- list(
            dada2 = file.path(pipeline_dir, "dada2/dada_mrexp.rds"),
            mothur =  file.path(pipeline_dir, "mothur/mothur_mrexp.rds"),
            qiime =  file.path(pipeline_dir, "qiime/qiime_mrexp.rds"),
            unclustered =  file.path(pipeline_dir, "mothur/unclustered_mrexp.rds")
      )
      
      mrexp_files %>% map(readRDS)
}


## Loading Phyloseq objects ------------------------------------------------
## argument - path to mgtst_pipelines
## returns a list pipelines objects
get_phyloseq <- function(pipeline_dir){
      ps_files <- list(
            dada2 = file.path(pipeline_dir, "dada2/dada_ps.rds"),
            mothur =  file.path(pipeline_dir, "mothur/mothur_ps.rds"),
            qiime =  file.path(pipeline_dir, "qiime/qiime_ps.rds"),
            unclustered =  file.path(pipeline_dir, "mothur/unclustered_ps.rds")
      )
      
      ps_files %>% map(readRDS)
}

