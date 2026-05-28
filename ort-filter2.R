library(readxl)
library(writexl)
library(rtracklayer)
library(tidyverse)

# Used in extracting upstream regions (promoter regions)
library(seqinr)

# Used for downloading genomes.
library(rentrez)
library(R.utils)

# Used for checking read count in each gene.
library(Rsubread)

# Used for 5' UTR graphs
library(ggplot2)

# Used for parsing .faa extension file (format: fasta, context: seqs, parser: read_seq_len)
# see documentation at https://cran.rstudio.com/web/packages/gggenomes/gggenomes.pdf
# library(gggenomes)
# nevermind, it didn't work, ill be using seqinr instead
library(seqinr)

start_wd <- "~/files/filesaguia-orthofinder/orthofinder/analise"
setwd(start_wd)

faa_folder <- "faa"
orthologues_folder <- "Orthologues"
genome_folder <- "genomes"

dados <- read_tsv("dados.tsv")
faas <- list.files(path = file.path(paste0("./", faa_folder)), pattern = "\\.faa$")
orthos <- list.files(path = file.path(paste0("./", orthologues_folder)), pattern = "\\.tsv$", recursive = TRUE)

# we have 4 datasets:
# [1] TSS list acquired from TSSAR, containing gene names
# [2] .gbff file containing annotations about genes and proteins
# [3] .faa file containing protein sequences
# [4] Orthologues between every pair, as a .tsv

# first, i'll be mapping genes of [1] that have proteins listed in [2]
# then, i'll be mapping proteins of [1+2] that have described sequences in [3]

filesuffix <- "pTSS" # altss / ptss

for (i in 1:nrow(dados)) {
  filename <- dados$Filenames[i]
  assembly <- dados$Assembly[i]
  print(paste("For file", filename))
  
  allowed <- data.frame(
    protein = as.character(),
    gene = as.character()
  )
  
  gbff <- read_gbk(file.path(start_wd, genome_folder, paste0(filename, ".gbff")))
  gbff <- gbff |> filter(type == "CDS")
  altss <- read_excel(file.path(start_wd, "altss", paste0(filename, "-", filesuffix, ".xlsx")))
  altssuq <- altss |> group_by(Gene) |> slice_head() |> ungroup()
  # uniq genes
  
  for (j in 1:nrow(altssuq)) {
    np <- {gbff |> filter(feat_id == paste0("cds-", altssuq$Gene[j])) |> select(protein_id)}$protein_id[1]
    # if (is.na(np)) {
    #   print(paste0("Warning: Gene ", altss$Gene[j], " doesn't have annotated protein"))
    # }
    allowed <- allowed |>
      add_row(
        protein = np,
        gene = altssuq$Gene[j]
      )
  }
  
  altss$Protein <- NA
  for (j in 1:nrow(altss)) {
    prot <- allowed |> filter(gene == altss[j,]$Gene)
    if (nrow(prot) > 0) {
      altss[j,]$Protein <- prot$protein
    }
  }
  
  # lists through all given TSSs, listing those that have proteins in .gbff files.
  # it is expected for some TSSs to not have proteins.
  allowed <- tibble({allowed |> filter(!is.na(allowed$protein))})
  print(paste0("Mapped proteins for ", 
               nrow(allowed), "/", 
               nrow(altssuq), " genes (", 
               round(100*nrow(allowed)/nrow(altssuq), 1), "%)" ))
  print(paste0("Mapped proteins for ",
               sum(!is.na(altss$Protein), na.rm = TRUE), "/",
               nrow(altss), " TSSs"))
  print(paste0("So in total, ",
               length(unique(altss$Protein[!is.na(altss$Protein)])),
               " unique proteins were found."))
  
  protList <- unique(altss$Protein[!is.na(altss$Protein)])
  
  orthoFiles <- orthos[grepl(paste0(filename, "(.*)\\/"), orthos)]
  
  altss$Orthogroup <- NA
  for (j in 1:length(orthoFiles)) {
    compSpecies <- sub("_(.*)", "", sub("(.*)__", "", orthoFiles[j]))
    compProts <- read_tsv(file.path(start_wd, orthologues_folder, orthoFiles[j]))
    colnames(compProts) <- c("Orthogroup", filename, compSpecies)
    
    altss[[paste0("ortho_", compSpecies)]] <- NA
    for (k in 1:nrow(altss)) {
      # searches for a matching protein from altss[k,] inside the .tsv from Orthologues
      if (is.na(altss[k,]$Protein)) {
        next
      }
      myProtein <- compProts |> filter(grepl(altss[k,]$Protein, compProts[[2]]))
      
      
      if (nrow(myProtein) > 0) {
        # print(paste0("Found protein for", altss[k,]$Protein))
        altss[k,]$Orthogroup <- myProtein[[1]]
        altss[k,][[paste0("ortho_", compSpecies)]] <- myProtein[[3]]
      } # else {
      #  print(paste0("No orthologues for", altss[k,]$Protein))
      #}
      
      
    }
  }
  
  print("So, how many alTSSs have orthologues?")
  print(sum(!is.na(altss$Orthogroup)))
  print(paste0("Out of all TSSs, ", sum(!is.na(altss$Orthogroup)),
               "/", nrow(altss), " have orthologues. "))
  
  print(paste0("Oh well, that's ", 
               length(unique(altss$Orthogroup[!is.na(altss$Orthogroup)])), 
               " Orthogroups out of ", length(unique(altss$Gene)), " genes."))
  # i dont trust this to work, because "WP_010870676.1" isnt working. idk why.
  
  write_xlsx(altss, path = file.path(start_wd, "altss", paste0(filename, "-", filesuffix, ".xlsx")))

  
  
  # now we check for which ones exist and dont excist and all that
  # now we read the .faa file.
  # faafile <- file.path(start_wd, paste0("protein_", assembly, ".faa"))
  # if (!file.exists(faafile)) {
  #   print(paste0("File .faa does not exist: ", faafile))
  #   next
  # }
  # 
  # foobar <- seqinr::read.fasta(file = faafile, seqtype = "AA")
  # 
  # tmp <- data.frame(
  #   name = as.character(),
  #   annot = as.character(),
  #   class = as.character(),
  #   seq = as.character(),
  #   mapped = as.logical()
  # )
  # for (j in 1:length(foobar)) {
  #   tmp <- tmp |> add_row(
  #     name = attr(foobar[[j]], "name"),
  #     annot = attr(foobar[[j]], "Annot"),
  #     class = attr(foobar[[j]], "class"),
  #     seq = paste(foobar[[j]], collapse = ""),
  #     mapped = attr(foobar[[j]], "name") %in% allowed$protein
  #   )
  # }
  # 
  # tmp <- tibble(tmp)
  # 
  # ptablemapped <- tmp |> filter(mapped)
  # 
  # print(paste0("Read .faa file, finding ", nrow(tmp), " proteins."))
  # 
  # print(paste0("Out of all alTSSs proteins, ", length(tmp$mapped[tmp$mapped]), "/", nrow(allowed), " were mapped to existing proteins in .faa file."))
  # 
  # fileout <- file.path(start_wd, paste0("filtered_", assembly, ".faa"))
  # 
  # outputlist <- list()
  # for (j in 1:nrow(ptablemapped)) {
  #   outputlist[[ sub(">", "", ptablemapped$annot[j]) ]] <- ptablemapped$seq[j]
  # }
  # 
  # 
  # write.fasta(outputlist, names(outputlist), fileout, as.string = TRUE, nbchar = 80)
  # 
}
