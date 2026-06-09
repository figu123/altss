

# script created to filter "candidates".
# here i'm fetching for all TSSs that:
#  1. have orthogroups
#  2. has more than >= 2 TSSs
#  3. at least one of their orthogroup counterpart also has >= 2 TSS
# my final list has 282 candidates = 282 genes that have >= 2 pTSS in more than 1 species

dados <- read_tsv("./dados.tsv")

df <- tibble(filename = as.character())
for (i in 1:nrow(dados)) {
  fn <- dados[i,]$Filenames
  
  ptss <- read_xlsx(paste0("./altss/", fn, "-pTSS.xlsx"))
  ptss$filename <- fn
  
  df <- bind_rows(df, ptss)
}

dff <- df %>% filter(!is.na(Orthogroup)) %>% group_by(Orthogroup) %>%
       filter(TSSCount > 1) %>%
       filter(n_distinct(Chromossome) >= 2) %>% ungroup()

uniqs <- dff %>% group_by(Orthogroup) %>% slice_head(n = 1)

paste0("Number of Orthogrous with 2 or more organisms with pTSS > 1: ", nrow(uniqs))

View(dff)
View(uniqs)

write_xlsx(uniqs, "./uniqs.xlsx")
write_xlsx(dff, "./all.xlsx")



