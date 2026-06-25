library(readxl)

files <- list.files("./altss", pattern="\\.xlsx$", full.names=TRUE)

for (i in 1:length(files)) {
  fn <- files[i]
  exceltable <- read_xlsx(fn)
  write_tsv(exceltable, sub("\\.xlsx", "\\.tsv", fn))
  
}

