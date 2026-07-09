
library(tidyverse)

# just a simple script to list the amount of alTSS/pTSS on different p-values.

dados <- read_tsv("./data/dados.tsv")

# MAKE SURE THIS IS LOWER THAN OR EQUAL TO THE THRESHOLD USED IN ALTSS.R
pvalue_threshold <- 1e-6 

summary <- data.frame(
  Filenames = as.character(),
  Species = as.character(),
  Pvalue = as.numeric(),
  TSS = as.numeric(),
  pTSS = as.numeric(),
  gTSS = as.numeric(),
  alTSS = as.numeric(),
  Genes = as.numeric()
)
for (i in 1:nrow(dados)) {
  fn <- dados[i,]$Filenames
  
  print(paste0("For ", fn))
  
  merged <- read_tsv(paste0("./results/merged/", fn, ".tsv")) %>% 
    filter(Pvalue < pvalue_threshold)
  
  ptss <- read_tsv(paste0("./results/altss/", fn, "-pTSS.tsv")) %>% 
    filter(Pvalue < pvalue_threshold)
  
  summary <- summary %>% add_row(
    Filenames = fn,
    Species = dados[i,]$Species,
    Pvalue = pvalue_threshold,
    TSS = nrow(merged),
    pTSS = nrow({ptss}),
    gTSS = nrow({ptss %>% filter(SubClass == "gTSS")}),
    alTSS = nrow({ptss %>% filter(SubClass == "alTSS")}),
    Genes = nrow({ptss %>% group_by(Gene) %>% slice_head(n = 1)})
  )
}
write_tsv(summary, "./results/summary2.tsv")
