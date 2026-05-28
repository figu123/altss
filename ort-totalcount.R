

library(readxl)
library(writexl)
library(tidyverse)

for (i in 1:nrow(dados)) {
  df <- read_excel(paste0("./altss/", dados[i,]$Filenames, "-alTSS.xlsx"))
  df$total_ortho <- 0
  
  for (k in 1:nrow(dados)) {
    # if (dados[k,]$Filenames == dados[i,]$Filenames) { next; }
    tmpfoo <- paste0("ortho_", dados[k,]$Filenames)
    if (tmpfoo %in% colnames(df)) {
      df$total_ortho = df$total_ortho + !is.na({df %>% select(tmpfoo)}[[tmpfoo]]);
    }
  }
  
  write_xlsx(df, paste0("./altss/", dados[i,]$Filenames, "-alTSS.xlsx"))
}
