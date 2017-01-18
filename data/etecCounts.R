## ETEC counts from Pop et al 2016

get_etecCounts <- function(){
      load("data/raw/etecExp.rda")
      
      ## experiment count data
      pd <- pData(etecExp)
      
      ## log2 qPCR data --------------------------------------------------------------
      qpcr <- log2(as.numeric(as.character(pd$QPCRLT)) + 1)
      o <- order(pd$Day,pd$Subj)
      days <- unique(pd$Day)
      sampIds <- unique(pd$Subj)
      
      tmp <- matrix(NA, nr=length(sampIds), nc=length(days))
      rownames(tmp) <- sampIds
      colnames(tmp) <- days
      for (d in seq(along=days)) {
          ii <- pd$Day == days[d]
          m <- match(pd$Subj[ii], sampIds)
          tmp[m,d] <- qpcr[ii]
      }
      qpcr_df <- tmp %>% as.data.frame() %>%
          mutate(ID = rownames(.), method = "qPCR")
      
      #maxday <- apply(qpcr,1,which.max)
      
      ## log2 qPCR data --------------------------------------------------------------
      ecoliFeatures <- which(fData(etecExp)$species == "Escherichia coli")
      tmp <- log2(MRcounts(etecExp[ecoliFeatures,], norm=TRUE) + 1)
      
      # # changes day to -1 to match ngs day
      ind <- sapply(qpcr_df$ID, function(x) which(pd$Subj == x & pd$Day == -1))
      drop <- sapply(ind, function(x) length(x)==0)
      ind <- unlist(ind[!drop])
      
      ind <- sapply(names(ind),
                    function(x) which(pd$Subj == x))# & pd$Day == colnames(qpcr)[maxday[x]]))
      
      cnts <- tmp
      aggTmp <- colSums(cnts,na.rm=TRUE)
      
      tmp <- matrix(NA, nr=length(ind), nc=length(days))
      rownames(tmp) <- names(ind)
      colnames(tmp) <- days
      for (d in seq(along=days)) {
          ii <- which(pd$Day == days[d])
          m <- match(pd$Subj[ii], names(ind))
          tmp[m[!is.na(m)],d] <- aggTmp[ii[!is.na(m)]]
      }
      ngs_df <- tmp %>% as.data.frame() %>%
          mutate(ID = rownames(.), method = "NGS")
      
      
      
      sample_etec_counts <- bind_rows(qpcr_df, ngs_df) %>%
          gather(key= "time",value = "count",-ID, -method) %>%
          mutate(time = as.numeric(time)) %>%
          spread(method, count)
      
      write_csv(sample_etec_counts, "data/raw/sample_etec_counts.csv")
      
      sample_etec_counts
}

etecCounts <- get_etecCounts()

ProjectTemplate::cache("etecCounts") 

rm(get_etecCounts)
