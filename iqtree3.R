library(tidyverse)


# Generating iqtreeinput ------
dados <- read_tsv("./data/dados.tsv")
dados <- dados |> mutate(label = paste0(toupper(substr(Filenames, 1, 1)), ". ", substr(Filenames, 2, nchar(Filenames))))

df <- tibble(filename = as.character(), label = as.character())
for (i in 1:nrow(dados)) {
  fn <- dados[i,]$Filenames
  ptss <- read_tsv(paste0("./results/altss/", fn, "-pTSS.tsv"))
  ptss$filename <- fn
  ptss$label <- dados[i,]$label
  df <- bind_rows(df, ptss)
}

# Orthogroups of all genes mapped by dRNA-seq.
orthogroups_all <- df %>% 
  group_by(Orthogroup, filename) %>% 
  slice_head(n = 1) %>%
  select(filename, Orthogroup) %>%
  filter(!is.na(Orthogroup))

length(unique(orthogroups_all$Orthogroup))

# Orthogroups of multiple-TSSs-genes (So only genes with multiple TSSs are considered)
orthogroups_with_altss <- df %>% 
  group_by(Orthogroup, filename) %>% 
  filter(SubClass == "alTSS") %>%
  slice_head(n = 1) %>%
  select(filename, Orthogroup) %>%
  filter(!is.na(Orthogroup))

length(unique(orthogroups_with_altss$Orthogroup))
# Here, there will be Orthogroups that only contain 1 multiple-TSS-gene. Check it out:
# # View(orthogroups_with_altss %>% ungroup() %>% group_by(Orthogroup) %>% summarize(count = n()))

# Orthogroups of multiple-TSSs-genes (So only genes with multiple TSSs are considered)
# where each orthogroup has at least 2 members == CANDIDATES.TSV
orthogroups_with_shared_altss <- df %>%
  group_by(Orthogroup, filename) %>% 
  filter(SubClass == "alTSS") %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(Orthogroup) %>%
  # summarize(count = n())
  filter(n() > 1) %>%
  select(filename, Orthogroup) %>%
  filter(!is.na(Orthogroup))

length(unique(orthogroups_with_shared_altss$Orthogroup))
# Here, at minimum Orthogroups will contain 2 multiple-TSS-gene. Check it out:
# # View(orthogroups_with_shared_altss %>% ungroup() %>% group_by(Orthogroup) %>% summarize(count = n()))


# Same amount of orthogroups as before, but now extend for genes that don't have multiple TSSs.
# So you only increased the number of genes, now it includes genes that only have 1 TSS.
# == CANDIDATES_BIG.TSV
orthogroups_with_shared_altss_big <- df %>%
  group_by(Orthogroup, filename) %>% 
  filter(Orthogroup %in% unique(orthogroups_with_shared_altss$Orthogroup)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(filename, Orthogroup) %>%
  filter(!is.na(Orthogroup))

length(unique(orthogroups_with_shared_altss_big$Orthogroup))

# IMPORTANT: 
# Select one from the 4 previous orthogroup variables to generate a file.
selected <- orthogroups_all
selected_fileout <- "input_orthogroups_all.fasta"

selected <- orthogroups_with_altss
selected_fileout <- "input_orthogroups_with_altss.fasta"

selected <- orthogroups_with_shared_altss
selected_fileout <- "input_candidates.fasta"

selected <- orthogroups_with_shared_altss_big
selected_fileout <- "input_candidates_big.fasta"



# Each column corresponds to an orthogroup.
# Each line corresponds to a species.
# Each cell can have TRUE or FALSE, if the gene for that corresponding line-species exists in that column-orthogroup.
# Further filters can exist, depending if your starting collection of genes comes from "orthogroups_with_shared_altss" or not.
# (Which filters only for genes that contain alTSS, that is, multiple TSSs)
species_og_df <- selected %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = Orthogroup, values_from = value)
species_og_df <- species_og_df %>% as.data.frame()
rownames(species_og_df) <- species_og_df$filename
species_og_df$filename <- NULL

# The table below maps every single orthogroup to either "PRESENT" or "ABSENT" for all species.
# It's the same as before, but its functionally equivalent to a list, 
# since the table is just a very long mapping of 1 key to 1 value:
# What species (filename) maps to what Orthogroup (Orthogroup) and is it present or not? (value)
sp_og_presence <- selected %>%
  mutate(value = TRUE) %>% 
  ungroup() %>%
  complete(filename, Orthogroup, fill = list(value = FALSE)) %>%
  arrange(Orthogroup)

# From the table above, you can create the .fasta that is used by iqtree to generate the tree.
# You're just converting the information to a string, where 1 is presence, 0 is absence.
# Then running IQTree binary to generate the tree.
seqs <- sp_og_presence %>%
  group_by(filename) %>%
  summarise(sequence = paste0(as.integer(value), collapse = ""))

file <- file(paste0("./results/iqtreeinput4/", selected_fileout), "w")

for (i in 1:nrow(seqs)) {
  writeLines(paste0(">", seqs$filename[i]), file)
  writeLines(seqs$sequence[i], file)
}

close(file)





# OLD
# i can try this in 3 ways
# 1. check all orthologues for the data i have (don't filter for pTSS)
# 2. check orthologues for all orthogroups with at least 1 case of multiple pTSS
# 3. check orthologues for all orthogroups with at least 2 caes of multiple pTSS

# this is really convoluted, not sure if i needed to create a new table from scratch
# but at least it'll be safe

# for starters, i need a table with all the genes that my sample has
# (or well, all orthogroups)
# correspondence_orthogroups1 <- df %>% 
#   group_by(Orthogroup) %>% # filter(TSSCount > 1) %>%
#   slice_head(n = 1) %>%
#   pull(Orthogroup)
# 
# valid_orthogroups <- list()
# 
# for (espec in dados$Filenames) {
#   valid_espec <- unique(df %>% filter(filename == espec) %>% pull(Orthogroup))
#   valid_orthogroups[[espec]] <- valid_espec[!is.na(valid_espec)]
# }
# 
# 
# orthotable1 <- tibble(
#   Orthogroup = as.character(),
#   filename = as.character(),
#   exists = as.logical()
# )
# for (ortho in correspondence_orthogroups1) {
#   for (espec in dados$Filenames) {
#     orthotable1 <- orthotable1 %>% add_row(
#       Orthogroup = ortho,
#       filename = espec,
#       exists = ortho %in% valid_orthogroups[[espec]]
#     )
#   }
# }
# 
# orthotable1_final <- orthotable1 %>% 
#   pivot_wider(names_from = "filename", values_from = "exists")
# 





# For IQTree Model selection, this is quite a hard task.

# There are two models for binary data, JC2 and GTR2.
# https://iqtree.github.io/doc/Substitution-Models#binary-and-morphological-models

# JC2 has some assumptions about the data, mainly, that transitions between states have equal probabilities.
# It seems like this model was used when studying SNPs, hence this assumption.
# https://www.puc-rio.br/ensinopesq/ccpg/pibic/relatorio_resumo2008/resumos/ctc/mat/mat_a.pdf

# But my current data is very "weird", and by that, I mean:

# For each orthogroup, I am assigning:
# 1 = if the species has a gene in this Orthogroup which was mapped by dRNA-Seq + the gene has multiple TSSs
# 0 = otherwise.

# Therefore, each "0" in my data can be:
# I. Not detected as an expressed gene in the dRNA-Seq dataset;
# II. A gene that couldn't be mapped to an Orthogroup (eg. if its absent for that species);
# III. Gene that doesn't contain alternative TSSs

# Note: We have filtered for orthogroups that only have alternative TSSs, 
#       so this limits the total number of groups we're looking at.
#       All other genes of the species (that don't have alternatives)
#       have been discarded.

# Therefore, I think it is out of the question to add assumptions on 
# my tree, since I'm not expecting transformations to be symmetrical.

# # e.g. Transformations (absence of a gene mapped to an orthogroup) could 
# # merely be a consequence from the sequencing depth of the lab...

# That's why we're using GTR2 model.

#
# #!/bin/bash
# source ~/miniconda3/bin/activate
# conda activate iqtree
# iqtree -s input.fasta -st BIN -m GTR2

# from here, im trying to install ggtree on my R version but its not working
# so now im trying to use renv

# # for things in CRAN, looks like this works.
# renv::install("BiocManager")
# renv::install("tidyverse")
# renv::install("circlize")
# renv::install("ggthemes")
# renv::install("gridExtra")
# renv::install("seqinr")

# # otherwise, i specify the lib

# lib <- .libPaths()[1]
# # make sure that lib points to your renv
# BiocManager::install("ggtree", lib = lib)
# BiocManager::install("ComplexHeatmap", lib = lib)
# BiocManager::install("Rsubread", lib = lib)
# BiocManager::install("rtracklayer", lib = lib)
# BiocManager::install("UniProt.ws", lib = lib)

# Load libs and functions
library(ggtree) # ggtree 3.14.0 requires ggplot2 3.5.2
library(treeio)
library(dplyr)
library(ggplot2)
library(reshape)
library(tidyr)
library(viridisLite)
library(RColorBrewer)
library(ggthemes)



# From here on out, make sure the tree you're reading (eg. input_candidates_big)
# and the species-orthogroup-dataframe generated in the previous bit of code (species_og_df)
# correspond to each other. this allows you to both visualize the tree made by ggtree, 
# and the table with all the data used to generate the tree.

# TODO: Pass these through a function, instead of manually assigning it

# read tree
tre = read.tree("./results/iqtreeinput4/input_candidates_big.fasta.treefile")
# data = read.table("../VF241020_parsed_table.csv", sep=";", header=T) %>% as_tibble()

# Plot tree, left side
p = ggtree(tre) + geom_tiplab(size=3, align=TRUE, linesize=.5)

# Parse dataframe
# lbls = data$label
# df = data %>% select(-label) %>% as.data.frame()
# rownames(df) = lbls

# Plot table
gheatmap(p, species_og_df, offset=1.5, width=3.6, color="black",
         colnames=TRUE, colnames_angle = 90, font.size = 2, colnames_offset_y = 0,
         legend_title="genotype") +
  scale_x_ggtree() +
  scale_fill_manual(
    values = c(
      rev(ggthemes::tableau_color_pal(palette = "Tableau 20")(11)),
      "gray"
    ), 
    na.value = "white"
  )





# Writing down the results.

#
# 1. orthogroups_with_shared_altss_big, or candidates_big.
# It's okay-ish!
# 
# The ones that clustered together:
#   hvolcanii and hsalinarum
#   sacidocaldarius and mjannaschii
# The ones further away:
#   tonnurineus and hmediterranei
#
# It makes sense for these two to be further away, I mean, there's barely
# any data for them (they're the ones with the lowest amount of pTSSs really)
#
# It's just a bit beyond me that pfuriosus and tkodakarensis don't cluster here,
# since they're both thermococcus, and both had good sequencing depth.
# 
# hvolcanii and hsalinarum being close together is ok.
#
#

# 
# 2. orthogroups_with_altss
# nope. this one is too polluted. don't use this one.
#


# 
# 3. orthogroups_all
# it's exactly what you'd expect from the taxonomy!
# but at this point you're not filtering for alternative TSSs, so whats the point.
# it's basically the same as candidates_big, but... the thermococcus are clustered.
#

#
# 4. candidates
# not good.
# i think some species like hsalinarum, tkodakarensis and tonnurineus don't have 
# enough genes to be able to cluster with other species here.
#

