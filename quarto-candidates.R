

# script created to filter "candidates".
# here i'm fetching for all TSSs that:
#  1. have orthogroups
#  2. has more than >= 2 TSSs
#  3. at least one of their orthogroup counterpart also has >= 2 TSS
# my final list has 282 candidates = 282 genes that have >= 2 pTSS in more than 1 species

library(tidyverse)

dados <- read_tsv("./data/dados.tsv")

df <- tibble(filename = as.character())
for (i in 1:nrow(dados)) {
  fn <- dados[i,]$Filenames
  
  ptss <- read_tsv(paste0("./results/altss/", fn, "-pTSS.tsv"))
  ptss$filename <- fn
  
  df <- bind_rows(df, ptss)
}

# df contains all pTSS observations thus far.
# you can apply further filters to it. example:
#df <- df %>% filter(Pvalue < 1e-6) %>% 
#  filter(gTSSDistance == 0 | abs(gTSSDistance) > 6) %>%
#  group_by(Gene) %>% mutate(TSSCount = n()) %>% ungroup()

# candidates filters for pTSSs that: 
#    (1) has alternative TSSs (TSS Count > 1 for that gene);
#    (2) has Orthogroup (has orthology with the gene of another species);
#    (3) within that orthogroup, you can find at least 1 other species that has alternative TSSs (that passed through the first filter)

filtered_candidates <- df %>% 
  filter(TSSCount > 1) %>% 
  filter(!is.na(Orthogroup)) %>% 
  group_by(Orthogroup) %>% filter(n_distinct(filename) >= 2) %>% ungroup()

uniqs <- filtered_candidates %>% group_by(Orthogroup) %>% slice_head(n = 1)
# write_tsv(uniqs, "./uniqs.tsv")
paste0("Number of ORTHOGROUPS with 2 or more organisms with pTSS > 1: ", nrow(uniqs))
paste0("Or about ", nrow({filtered_candidates %>% group_by(Gene) %>% slice_head(n = 1)}), " genes.")

# unshared_alternatives are pTSS that:
#    (1) has alternative TSSs for that gene;
#    (2) has Orthogroup (has orthology with the gene of another species)
#    (3) within that Orthogroup, they are the only species that have alternative TSSs.
unshared_alternatives <- df %>% 
  filter(TSSCount > 1) %>% 
  filter(!is.na(Orthogroup)) %>% 
  group_by(Orthogroup) %>% filter(n_distinct(filename) == 1) %>% ungroup()

uniqs2 <- unshared_alternatives %>% group_by(Orthogroup) %>% slice_head(n = 1)
paste0("Number of ORTHOGROUPS with only 1 organism with pTSS > 1: ", nrow(uniqs2))
paste0("Or about ", nrow({unshared_alternatives %>% group_by(Gene) %>% slice_head(n = 1)}), " genes.")

# finally, there's the "non-orthogroupped alTSSs"
# i.e. it's about genes that have alternatives, but they don't have orthologues.
ungrouped_alternatives <- df %>%
  filter(TSSCount > 1) %>% 
  filter(is.na(Orthogroup))

uniqs3 <- ungrouped_alternatives %>% group_by(Gene) %>% slice_head(n = 1)
paste0("Number of GENES with no orthology and pTSS > 1: ", nrow(uniqs3))

#negative_candidates <- df %>% filter(!is.na(Orthogroup)) %>% group_by(Orthogroup) %>%
#        filter(TSSCount == 1) %>% filter(n_distinct(filename) >= 2) %>% ungroup()

write_tsv(filtered_candidates, "./results/candidates.tsv")
# write_tsv(df, "./allmerged.tsv")

