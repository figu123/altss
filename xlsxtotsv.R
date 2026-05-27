library(readxl)

setwd("~/Desktop/quarto/test/results/altss/")

files <- list.files(".", pattern="\\.xlsx$", full.names=TRUE)

lapply(files, function(f){
  df <- read_excel(f)
  write.table(df, sub("\\.xlsx$", ".tsv", f),
              sep="\t", row.names=FALSE, quote=FALSE)
})
