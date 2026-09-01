
library(tidyverse)

faas <- list.files(path = file.path("./data/genomes"), pattern = "\\.faa$")
orthos <- list.files(path = file.path("./results/orthofinder/output/"), pattern = "\\.tsv$", recursive = TRUE)
orthos <- orthos[grepl("/Orthologues/", orthos)]

dados <- read_tsv("./data/dados.tsv", show_col_types = FALSE)
resultSummary <- read_tsv("./data/dados.tsv", show_col_types = FALSE)

resultSummaryGenes <- read_tsv("./data/dados.tsv", show_col_types = FALSE)
resultSummaryAltss <- read_tsv("./data/dados.tsv", show_col_types = FALSE)
resultSummaryAltss2 <- read_tsv("./data/dados.tsv", show_col_types = FALSE)

for (i in 1:nrow(dados)) {
  filename <- dados$Filenames[i]
  assembly <- dados$Assembly[i]
  print(paste("For file", filename))
  
  if (!file.exists(paste0("./results/altss/", filename, "-pTSS.tsv"))) {
    print(paste0("Could not find pTSS table for ", filename))
    next;
  }
  if (!file.exists(paste0("./data/genomes/", filename, "-features.tsv"))) {
    print(paste0("Could not find features table (genome) for ", filename))
    next;
  }
  
  features <- read_tsv(paste0("./data/genomes/", filename, "-features.tsv"), show_col_types = FALSE)
  names(features)[1] <- "feature"
  features <- features[features$feature == "CDS",]
  
  ptss <- read_tsv(paste0("./results/altss/", filename, "-pTSS.tsv"), show_col_types = FALSE)
  
  if ("Orthogroup" %in% colnames(ptss)) {
    print(paste0("This script was already run for ", filename))
    # running it multiple times would duplicate the number of columns each time. 
    # failsafe:
    ptss <- ptss %>% select(-total_ortho, -Orthogroup, -starts_with("ortho_"))
  }
  
  # ptssuq <- ptss |> group_by(Gene) |> slice_head() |> ungroup()
  # a. Defining "allowed" df:
  # allowed <- unique(ptss$Gene)
  # protList <- unique(ptss$Protein[!is.na(ptss$Protein)])
  
  print("Listing Results_*/Orthologues/Orthologues_{filename}/*.tsv")
  speciesOrthologuesFiles <- orthos[grepl(paste0(filename, "(.*)\\/"), orthos)]
  
  ptss$total_ortho <- 0
  ptss$Orthogroup <- NA
  for (j in 1:length(speciesOrthologuesFiles)) {
    otherSpecies <- sub("(.*)__v__", "", speciesOrthologuesFiles[j])
    speciesOrthologuesTable <- read_tsv(paste0("./results/orthofinder/output/",
                                               speciesOrthologuesFiles[grepl(otherSpecies, speciesOrthologuesFiles)]),
                                        show_col_types = FALSE)
    otherSpecies <- sub(".tsv", "", otherSpecies)
    colnames(speciesOrthologuesTable) <- c("Orthogroup", "Protein", paste0("ortho_", otherSpecies))
    
    # 2nd column is the Protein from CurrentSpecies == Filename.
    # This column may have multiple proteins separated by commas, so 
    # below I make the table longer -- each row gets one protein.
    speciesOrthologuesTable <- speciesOrthologuesTable %>%
      separate_rows(Protein, sep = ",") %>%
      mutate(Protein = trimws(Protein))
    
    # combine tables feature from dplyr 
    # https://nyu-cdsc.github.io/learningr/assets/data-transformation.pdf
    # + coalesce(data.y, data.x) inside mutate(). 
    #   == take data.y (new table value) if it exists otherwise keep data.x (original value)
    # + select (to remove old data.)
    ptss <- ptss %>%
      left_join(speciesOrthologuesTable, by = "Protein") %>%
      mutate(Orthogroup = coalesce(Orthogroup.y, Orthogroup.x)) %>%
      select(-Orthogroup.y, -Orthogroup.x)
    
    ptss$total_ortho <- ptss$total_ortho + ptss$Protein %in% speciesOrthologuesTable$Protein
  }
  
  ptss <- ptss %>% relocate(Orthogroup, .after = total_ortho)
  
  summarySpecies <- ptss %>% 
    summarise(across(everything(), ~ sum(!is.na(.)))) %>%
    mutate(Filenames = filename) %>%
    select(Filenames, Protein, Orthogroup, starts_with("ortho_"))
  # na_summary <- ptss %>% summarise(across(everything(), ~ sum(!is.na(.))/n()))
  
  # filtered for genes ===  1 gene = 1 tss  === n(gTSS)
  summarySpeciesGenes <- ptss %>% 
    filter(SubClass == "gTSS") %>%
    summarise(across(everything(), ~ sum(!is.na(.)))) %>%
    mutate(Filenames = filename) %>%
    select(Filenames, Protein, Orthogroup, starts_with("ortho_"))
  
  # filtered only for genes that have multiple TSS
  summarySpeciesAltss <- ptss %>% 
    filter(TSSCount > 1) %>%
    summarise(across(everything(), ~ sum(!is.na(.)))) %>%
    mutate(Filenames = filename) %>%
    select(Filenames, Protein, Orthogroup, starts_with("ortho_"))
  
  # filtered only for genes that have multiple TSS, 1 gene == 1 tss
  summarySpeciesAltss2 <- ptss %>% 
    filter(TSSCount > 1 & SubClass == "gTSS") %>% 
    summarise(across(everything(), ~ sum(!is.na(.)))) %>%
    mutate(Filenames = filename) %>%
    select(Filenames, Protein, Orthogroup, starts_with("ortho_"))
  
  # the mutate below is quite verbose, but it merges columns & 
  resultSummary <- resultSummary %>%
    left_join(summarySpecies, by = "Filenames") %>%
    mutate(across(ends_with(".x"),
                  ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
                  .names = "{sub('.x$','',.col)}")) %>%
    select(-ends_with(".x"), -ends_with(".y"))
  
  resultSummaryGenes <- resultSummaryGenes %>%
    left_join(summarySpeciesGenes, by = "Filenames") %>%
    mutate(across(ends_with(".x"),
                  ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
                  .names = "{sub('.x$','',.col)}")) %>%
    select(-ends_with(".x"), -ends_with(".y"))
  
  resultSummaryAltss <- resultSummaryAltss %>%
    left_join(summarySpeciesAltss, by = "Filenames") %>%
    mutate(across(ends_with(".x"),
                  ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
                  .names = "{sub('.x$','',.col)}")) %>%
    select(-ends_with(".x"), -ends_with(".y"))
  
  resultSummaryAltss2 <- resultSummaryAltss2 %>%
    left_join(summarySpeciesAltss2, by = "Filenames") %>%
    mutate(across(ends_with(".x"),
                  ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
                  .names = "{sub('.x$','',.col)}")) %>%
    select(-ends_with(".x"), -ends_with(".y"))
  
  write_tsv(ptss, paste0("./results/altss/", filename, "-pTSS.tsv"))
}
resultSummary <- resultSummary %>% 
  select(Filenames, Protein, Orthogroup, starts_with("ortho_"))

# resultSummary2 <- resultSummary %>%
#   mutate(Orthogroup = Orthogroup/Protein) %>%
#   mutate(across(starts_with("ortho_"), ~ .x/Protein))

resultSummary3 <- resultSummaryGenes %>% 
  select(Filenames, Protein, Orthogroup, starts_with("ortho_"))

resultSummary4 <- resultSummaryAltss %>% 
  select(Filenames, Protein, Orthogroup, starts_with("ortho_"))

resultSummary5 <- resultSummaryAltss2 %>% 
  select(Filenames, Protein, Orthogroup, starts_with("ortho_"))

write_tsv(resultSummary, "./results/orthofinder-summary-TSSs.tsv")
write_tsv(resultSummary3, "./results/orthofinder-summary-genes.tsv")
write_tsv(resultSummary4, "./results/orthofinder-summary-alTSSs.tsv")
write_tsv(resultSummary5, "./results/orthofinder-summary-genes-alTSS.tsv")


