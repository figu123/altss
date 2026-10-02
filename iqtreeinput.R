
# run orthofinder.sh first, and parse its results for the Orthogroup table to exist.

library(tidyverse)

seqs <- read_tsv("./results/candidates-big.tsv")
seqs <- seqs %>% filter(actual_n_orthogroups >= 3)

dir.create("./results/iqtreeinput")

orthos <- sort(unique(seqs$Orthogroup))
orthos <- orthos[!is.na(orthos)]

print("Should I use all, candidates-big or candidates?")
# all == all merged TSSs, no filters.
# candidates filters for genes that have multiple TSSs (alternative promoter candidates, 2+ TSSs for 1 gene)
#   which have orthologues that ALSO HAVE multiple TSSs.
# candidates-big filters for genes that have multiple TSSs, but is more inclusive.
#   it includes all orthologues found for those genes, regardless if these orthologues have 2+ TSSs or not.
# So for example, for a cold regulating protein present in haloarchea, i can see:
#   H. mediterranei: HFX_RS02250 - 1 pTSS (only visible in candidates-big.tsv), 
#   H. volcanii: HVO_RS07105 - many pTSSs (visible on both)
#   H. salinarum: VNG_RS07115 - 2 pTSSs (visible on both)

print("Should I collapse duplicate TSSs for the same gene, or not?")
# to collapse, just instead of group_by(Gene), do group_by(filename)

print("Also, can each .fa have only 2 species, or do I need at least 3?")
# if so, then filter(actual_n_orthogroups >= 3)
# After trying to run IQ-Tree, this is what I got: ERROR: Alignment must have at least 3 sequences

for (i in 1:length(orthos)) {
  # generates .fa file for each orthogroup.
  mydf <- seqs %>% filter(Orthogroup == orthos[i]) %>%
    group_by(filename) %>% arrange(d) %>% slice_head(n = 1)
  
  myfile <- file(paste0("./results/iqtreeinput/", orthos[i], ".fa"), open = "a")
  
  for (j in 1:nrow(mydf)) {
    writeLines(paste0("> ", mydf[j,]$filename), myfile)
    writeLines(mydf[j,]$Sequence, myfile)
  }
  
  close(myfile)
}


