## dataframe with values for rarefaction curve

## Modification of rarecurve function in vegan - only return dataframe and not plot
`rarecurve_dat` <-
      function(x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
               label = TRUE, col, lty, ...)
      {
            ## matrix is faster than data.frame
            x <- as.matrix(x)
            ## check input data: must be counts
            if (!identical(all.equal(x, round(x)), TRUE))
                  stop("function accepts only integers (counts)")
            ## sort out col and lty
            if (missing(col))
                  col <- par("col")
            if (missing(lty))
                  lty <- par("lty")
            tot <- rowSums(x)
            S <- specnumber(x)
            ## remove empty rows or we fail
            if (any(S <= 0)) {
                  message("empty rows removed")
                  x <- x[S > 0,, drop =FALSE]
                  tot <- tot[S > 0]
                  S <- S[S > 0]
            }
            nr <- nrow(x)
            ## rep col and lty to appropriate length
            col <- rep(col, length.out = nr)
            lty <- rep(lty, length.out = nr)
            ## Rarefy
            out <- lapply(seq_len(nr), function(i) {
                  n <- seq(1, tot[i], by = step)
                  if (n[length(n)] != tot[i])
                        n <- c(n, tot[i])
                  drop(rarefy(x[i,], n))
            })
            Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
            Smax <- sapply(out, max)
            ## set up plot
            # plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
            #      type = "n", ...)
            # ## rarefied richnesses for given 'sample'
            # if (!missing(sample)) {
            #       abline(v = sample)
            #       rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), y = z,
            #                                              xout = sample, rule = 1)$y)
            #       abline(h = rare, lwd=0.5)
            # }
            # ## rarefaction curves
            # for (ln in seq_along(out)) {
            #       N <- attr(out[[ln]], "Subsample")
            #       lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
            # }
            # ## label curves at their endpoitns
            # if (label) {
            #       ordilabel(cbind(tot, S), labels=rownames(x), ...)
            # }
            invisible(out)
      }

get_curve_df <- function(mrexp){
      count_tbl <- mrexp[,pData(mrexp)$dilution %in% c(-1,0)]@assayData$counts %>% t() 
      samID <- rownames(count_tbl)
      
      curve_dat <- rarecurve_dat(count_tbl, step = 500)
      
      curve_df <- curve_dat %>% set_names(samID) %>% map(data.frame) %>%
            map_df(rownames_to_column, var = "sample", .id = "samID") %>%
            mutate(sample = str_replace(sample, "N","") %>% as.numeric())
      
      colnames(curve_df)[3] <- "OTUs"
      
      pData(mrexp) %>% rownames_to_column(var = "samID") %>%
            mutate(pcr_16S_plate = as.numeric(pcr_16S_plate)) %>% 
            select(sampleID, samID, dilution, pcr_16S_plate) %>%
            right_join(curve_df)
}