

# script created to filter "candidates".
# here i'm fetching for all TSSs that:
#  1. have orthogroups
#  2. has more than >= 2 TSSs
#  3. at least one of their orthogroup counterpart also has >= 2 TSS

# should make two outputs

# candidates-big i only expanded to include all genes of that shared orthogroup.
# So for example, for a cold regulating protein present in haloarchea, i can see:
#   H. mediterranei: HFX_RS02250 - 1 pTSS (only visible in candidates-big.tsv), 
#   H. volcanii: HVO_RS07105 - many pTSSs (visible on both)
#   H. salinarum: VNG_RS07115 - 2 pTSSs (visible on both)


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

# and lastly, a consideration about orthogroups.
# total_ortho has the count of orthologues for the protein considering the WHOLE GENOME of species.
#  (same goes for other mentionings of orthologues).
# if I were to describe the ACTUAL COUNT of orthologues for my dataset, the number will be different.
# for example, consider a hypothetical sugar metabolism gene that has 5 orthologues between species.
# maybe in the mapped TSSs, only 2 orthologues were mapped, because the dRNA-seq only captured the expression of these 2 orthologues (only they had reads)
#   [[idk why im explaining this in this detail, you get it]]

df <- df %>% group_by(Orthogroup) %>% 
  mutate(actual_n_orthogroups = n_distinct(filename)) %>% 
  ungroup() %>%
  relocate(actual_n_orthogroups, .after = total_ortho)


filtered_candidates <- df %>% 
  filter(TSSCount > 1) %>% 
  filter(!is.na(Orthogroup)) %>% 
  group_by(Orthogroup) %>% filter(n_distinct(filename) >= 2) %>% ungroup()

uniqs <- filtered_candidates %>% group_by(Orthogroup) %>% slice_head(n = 1)
# write_tsv(uniqs, "./uniqs.tsv")
paste0("Number of ORTHOGROUPS with 2 or more organisms with pTSS > 1: ", nrow(uniqs))
paste0("Or about ", nrow({filtered_candidates %>% group_by(Gene) %>% slice_head(n = 1)}), " genes.")
# == nrow({filtered_candidates %>% filter(SubClass == "gTSS")})

orthogroups <- unique(filtered_candidates$Orthogroup)
expanded_filtered_candidates <- df %>% filter(Orthogroup %in% orthogroups)


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
# == nrow({unshared_alternatives %>% filter(SubClass == "gTSS")})



# 916 orthogroups inferred containing alTSS
print(paste0("There were ",
             nrow({df %>% group_by(Orthogroup) %>% filter(!is.na(Orthogroup) & SubClass == "alTSS") %>% slice_head(n = 1)}),
             " orthogroups inferred containing alTSSs.",
             " This should be the sum of the two previously defined datatables: ", 
             nrow(uniqs2), " + ", nrow(uniqs)))


# finally, there's the "non-orthogroupped alTSSs"
# i.e. it's about genes that have alternatives, but they don't have orthologues.
ungrouped_alternatives <- df %>%
  filter(TSSCount > 1) %>% 
  filter(is.na(Orthogroup))

uniqs3 <- ungrouped_alternatives %>% group_by(Gene) %>% slice_head(n = 1)
paste0("Number of GENES with no orthology and pTSS > 1: ", nrow(uniqs3))

#negative_candidates <- df %>% filter(!is.na(Orthogroup)) %>% group_by(Orthogroup) %>%
#        filter(TSSCount == 1) %>% filter(n_distinct(filename) >= 2) %>% ungroup()
# write_tsv(df, "./allmerged.tsv")



# 	Primary transcript start sites were mapped for 7692 genes across all eight datasets, of which 5150 genes (67%) had coding proteins that could be infered to orthogroups coded by OrthoFinder. Out of the initial 7692 genes, 1613 have multiple pTSSs (21%) being potentially regulated by alternative promoters. Lastly, 1100 genes that presented multiple pTSSs and were mapped into orthogroups (14%). The distribution of multiple TSSs per annotated gene aligns with previously described data on secondary transcription start sites (Zhou M et al., 2025).

print(paste0(nrow({df %>% filter(Class == "P") }),
             " primary transcript sites (pTSSs) were mapped for ",
             nrow({df %>% group_by(Gene) %>% slice_head(n = 1)}),
             " genes across all ", 
             nrow(dados), 
             " datasets. Of which, ",
             nrow({df %>% group_by(Gene) %>% slice_head(n = 1) %>% filter(!is.na(Protein))}),
             " have annotated proteins. Of these, ",
             nrow({df %>% group_by(Gene) %>% slice_head(n = 1) %>% filter(!is.na(Orthogroup))}),
             " proteins were inferred to orthogroups by Orthofinder. ",
             "Out of the initial ", nrow({df %>% group_by(Gene) %>% slice_head(n = 1)}), " genes, ",
             nrow({df %>% group_by(Gene) %>% filter(SubClass == "alTSS") %>% slice_head(n = 1)}),
             " have multiple TSSs being potentially regulated by alternative promoters.",
             " Lastly, ",
             nrow({df %>% group_by(Gene) %>% filter(SubClass == "alTSS") %>% slice_head(n = 1) %>% filter(!is.na(Orthogroup))}),
             " genes that presented multiple TSSs were mapped into orthogroups.",
             " Out of these ",
             nrow({df %>% group_by(Gene) %>% filter(SubClass == "alTSS") %>% slice_head(n = 1) %>% filter(!is.na(Orthogroup))}),
             " multiple pTSS orthogroup-mapped genes, ",
             nrow({filtered_candidates %>% group_by(Gene) %>% slice_head(n = 1)}),
             " are part of orthogroups where at least genes of two species present multiple pTSSs",
             " (that is, the presence of alternative TSSs is conserved in genes across these organisms == ARE LISTED IN CANDIDATES.TSV).",
             " The remaining ",
             nrow({unshared_alternatives %>% group_by(Gene) %>% slice_head(n = 1)}),
             " genes despite having multiple TSSs, are the only gene with multiple TSSs in their orthogroup.",
             " (that is, only one species presents a gene with multiple TSSs)"
             ))

print(paste0("The list of ", 
             nrow({filtered_candidates %>% group_by(Gene) %>% slice_head(n = 1)}), 
             " candidate genes with alternative TSSs as a conserved mechanism across species is available in ./results/candidates.tsv. ",
             " They were mapped to ",
             nrow(uniqs),
             " orthogroups (so that's the count of 'unique gene functions' you need to analyse)."
      ))

# Exporting candidates
write_tsv(expanded_filtered_candidates, "./results/candidates-big.tsv")
write_tsv(filtered_candidates, "./results/candidates.tsv")



# Exporting "list of all pTSSs"
write_tsv(df %>% 
            rename("Filenames" = "filename") %>% 
            left_join(dados), 
          "./results/all-pTSSs.tsv")

# Exporting "list of all orthogrouped genes with multiple pTSSs"
write_tsv(df %>% 
            group_by(Gene) %>% 
            filter(SubClass == "alTSS") %>% 
            slice_head(n = 1) %>% 
            filter(!is.na(Orthogroup)) %>% 
            rename("Filenames" = "filename") %>% 
            left_join(dados),
          "./results/all-orthogrouped-genes-with-multiple-pTSS.tsv")



# further verification that tables "filtered_candidates" and "unshared_alternatives"
#    include all multiple TSS orthogroup-mapped genes:

# abcdef <- c({filtered_candidates %>% group_by(Gene) %>% slice_head(n = 1)}$Gene,
#            {unshared_alternatives %>% group_by(Gene) %>% slice_head(n = 1)}$Gene)
# length(abcdef) == nrow({df %>% group_by(Gene) %>% filter(SubClass == "alTSS") %>% slice_head(n = 1) %>% filter(!is.na(Orthogroup))})

# okay it seems to be alright


