library(readxl)

files <- list.files("./altss", pattern="\\.xlsx$", full.names=TRUE)

lapply(files, function(f){
  df <- read_excel(f)
  write.table(df, sub("\\.xlsx$", ".tsv", f),
              sep="\t", row.names=FALSE, quote=FALSE)
})
