## Extra tidy functions

## for defining column names in a pipe
set_colnames <- function(df, col_names){
      colnames(df) <- col_names
      df
}

## Not in opperator
"%!in%" <- function(x,y){
      !(x %in% y)
}

