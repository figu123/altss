
# run multipetssorthogroups.R first.


library(tidyverse)

seqs <- read_tsv("all.tsv") %>% group_by(Orthogroup) %>% filter(n_distinct(filename) > 2) %>% ungroup()

dir.create("./iqtreeinput")

orthos <- sort(unique(seqs$Orthogroup))
orthos <- orthos[!is.na(orthos)]

print("Should I collapse the duplicates, or not?")
# to collapse, just instead of group_by(Gene), do group_by(filename)

print("Also, can each .fa have only 2 species, or do I need at least 3?")
# if so, then filter for group_by(Orthogroup)

for (i in 1:length(orthos)) {
  mydf <- seqs %>% filter(Orthogroup == orthos[i]) %>%
    group_by(filename) %>% arrange(d) %>% slice_head(n = 1)
  
  myfile <- file(paste0("./iqtreeinput/", orthos[i], ".fa"), open = "a")
  
  for (j in 1:nrow(mydf)) {
    writeLines(paste0("> ", mydf[j,]$filename), myfile)
    writeLines(mydf[j,]$Sequence, myfile)
  }
  
  close(myfile)
}






