# last edit: 2026 - march - 05

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
library(gridExtra)
library(ggplot2)

# Used for parsing gbff file
library(gggenomes)

setwd("~/files/analises/full")
start_wd <- getwd()
print(paste0("Running at ", getwd()))

dados <- read_tsv("./dados.tsv")
genome_folder <- "./genomes"


dados <- dados |> mutate(label = paste0(toupper(substr(Filenames, 1, 1)), ". ", substr(Filenames, 2, nchar(Filenames))))



# OLD
# Making a new column of the distance to the associated gene for each pTSS.
# I suppose this region corresponds to 5'UTR region
# And transcripts that have small lengths (0-5nc) are leaderless transcripts.
# The size of this region depends on the organism you're working with.
# More information at doi: 10.1186/s12864-016-2920-y


# print("Misc/ graphs")
# ggdf <- data.frame(
#   nextFeatureDistance = as.integer(sub(".*?(\\d+)nt upstream.*", "\\1", filtered$Location)),
#   intensityTSS = filtered$LibP + filtered$LibM, 
#   intensityAll = filtered$LibP_Total + filtered$LibM_Total
# ) %>% 
#   mutate(utr = ifelse(nextFeatureDistance > 5, "5UTR", "Leaderless"))
# 
# ggdf2 <- data.frame(
#   gene_count = c(
#     length(unique({filtered %>% filter(TSSCount == 1)}$Gene)),
#     length(unique({filtered %>% filter(TSSCount == 2)}$Gene)),
#     length(unique({filtered %>% filter(TSSCount == 3)}$Gene)),
#     length(unique({filtered %>% filter(TSSCount > 3)}$Gene))
#   ),
#   tss_count = c("1 TSS", "2 TSSs", "3 TSSs", "4+ TSSs")
# )
# fooplot <- ggplot(ggdf, aes(nextFeatureDistance)) + geom_bar() + xlim(-1, 21)
# ggsave(file.path(misc_folder, paste0(filename, "-leaderless.png")), plot = fooplot)
# fooplot <- ggplot(ggdf[ggdf$nextFeatureDistance > 5, ], aes(nextFeatureDistance)) + geom_histogram(binwidth=10)
# ggsave(file.path(misc_folder, paste0(filename, "-5utr.png")), plot = fooplot)
# fooplot <- ggplot(ggdf2, aes(tss_count, gene_count)) + geom_col()
# ggsave(file.path(misc_folder, paste0(filename, "-altsscount.png")), plot = fooplot)





##### Boxplotting the Coverage of Genes (RSubread) #####

values_unpaired <- list()
values_paired <- list()

for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  
  runfor <- function(pairtype="unpaired") {
    fp <- file.path(start_wd, "readcount", paste0(filename, "-", pairtype, ".csv"))
    if (!file.exists(fp)) {
      print(paste0(fp, " does not exist, skipping."))
      return(c(0))
    }
    mycsv <- t(read_csv(fp))
    colnames(mycsv) <- mycsv[1,]
    mycsv <- mycsv[-1,]
    mycsv <- tibble(data.frame(mycsv)) |> 
      mutate(across(everything(), as.numeric)) |>
      summarise(across(everything(), sum))
    mycsv <- data.frame(t(mycsv))
    colnames(mycsv) <- "Reads"
    mycsv$Gene <- sub("gene.", "", rownames(mycsv))
    mycsv <- tibble(mycsv) |> relocate(Gene, Reads)
    
    return(mycsv$Reads)
  }
  values_unpaired[[filename]] <- runfor("unpaired")
  values_paired[[filename]] <- runfor("paired")
}

df <- stack(values_unpaired)
colnames(df) <- c("GeneReadcount", "species")

ggplot(df, aes(x = species, y = GeneReadcount)) +
  geom_boxplot() + ggtitle("Each point is a gene. (unpaired)") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") + 
  ylim(0, 15000) # + geom_jitter(width = 0.3, size = 0.2)

df <- stack(values_paired)
colnames(df) <- c("GeneReadcount", "species")

ggplot(df, aes(x = species, y = GeneReadcount)) +
  geom_boxplot() + ggtitle("Each point is a gene. (paired)") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") + 
  ylim(0, 15000) # + geom_jitter(width = 0.3, size = 0.2)



###### TSSs Class x Gene Count x Sequence Depth ((DEPRECATED)) ######
# continuation: comparing gene amount for each species and stuff like that
dff <- data.frame(
  species = as.character(),
  totalReads = as.numeric(),
  totalGenes = as.numeric(),
  pTSS = as.numeric(),  # all TSSs that are classified with 250bp upstream of a gene. Equals to the sum gTSS + alTSS
  gTSS = as.numeric(),  # pTSS with highest intensity/signal for a gene (i.e. highest read count)
  alTSS = as.numeric(), # pTSS that don't have highest intensity/signal for a gene
  gsTSS = as.numeric(), # gTSS that don't have an alTSS. Single gTSSs.
  gcTSS = as.numeric()  # gTSS that have an alTSS. Coupled gTSSs.
)
for (i in 1:length(values_unpaired)) {
  filename <- names(values_unpaired)[i]
  total <- sum(c(values_unpaired[[i]], values_paired[[i]]))
  gff <- read_gff3(file.path(start_wd, genome_folder, paste0(filename, ".gff")))
  geneCount <- nrow({gff |> filter(type == "gene")})
  data <- read_xlsx(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx") ))
  nrow({data |> filter(TSSCount > 1 & SubClass == "gTSS")})
  nrow({data |> filter(TSSCount > 1 & SubClass == "alTSS")})
  nrow({data |> filter(SubClass == "gTSS")})
  
  dff <- dff |> add_row(
    species = filename,
    totalReads = total,
    totalGenes = geneCount,
    pTSS = nrow(data),
    gTSS = nrow({data |> filter(SubClass == "gTSS")}),
    alTSS = nrow({data |> filter(TSSCount > 1 & SubClass == "alTSS")}),
    gsTSS = nrow({data |> filter(TSSCount == 1 & SubClass == "gTSS")}),
    gcTSS = nrow({data |> filter(TSSCount > 1 & SubClass == "gTSS")})
  )
}

dff <- dff |> 
  mutate(label = paste0(toupper(substr(species, 1, 1)), ". ", substr(species, 2, nchar(species)))) |>
  mutate(size = round(totalReads/totalGenes))

tssToGenePlots <- list()

tssToGenePlots[[1]] <- ggplot(dff, aes(x = totalGenes, y = alTSS, label = label, color = label)) +
  geom_point(aes(size = size )) +
  #geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8), 
  #          vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1), 
  #          check_overlap = TRUE, fontface = "italic") + theme_bw() +
  # scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(1, 7)) +
  theme(legend.position = "none")  # Remove all legends


tssToGenePlots[[2]] <- ggplot(dff, aes(x = totalGenes, y = pTSS, label = label, color = label)) +
  geom_point(aes(size = size )) +
  #geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8), 
  #          vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1), 
  #          check_overlap = TRUE, fontface = "italic") + theme_bw() +
  # scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(1, 7))  +
  theme(legend.position = "none")

tssToGenePlots[[3]] <- ggplot(dff, aes(x = totalGenes, y = gTSS, label = label, color = label)) +
  geom_point(aes(size = size )) +
  #geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8), 
  #          vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1), 
  #          check_overlap = TRUE, fontface = "italic") + theme_bw() +
  # scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(1, 7))  +
  theme(legend.position = "none")

tssToGenePlots[[4]] <- ggplot(dff, aes(x = totalGenes, y = pTSS, label = label, color = species)) +
  geom_point(aes(size = size ), alpha = 0, show.legend = TRUE) +
  #geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8), 
  #          vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1), 
  #          check_overlap = TRUE, fontface = "italic") + theme_bw() +
  # scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(1, 7)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)),
         size = guide_legend(override.aes = list(alpha = 1))) +
  theme_void() +
  theme(legend.position = "right",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0))


tssToGeneGraph <- grid.arrange(
  tssToGenePlots[[1]],
  tssToGenePlots[[2]],
  tssToGenePlots[[3]],
  tssToGenePlots[[4]],
  widths = c(1, 1, 1, 1),
  nrow = 1,
  # nrow = 4,
  top = "Relationship between TSS count, gene count and sequencing depth"
)

ggsave(file.path(start_wd, "./misc/graphstsstogene.png"), plot = tssToGeneGraph, units = "px", width=1080, height=720, scale=3)




# v2

tssToGenePlots[[1]] <- ggplot(dff, aes(x = totalGenes, y = alTSS, label = label, color = label)) +
  geom_point(aes(size = size )) +
  geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8),
            vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1),
            check_overlap = TRUE, fontface = "italic") + theme_bw() +
  scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(4, 12)) +
  theme(legend.position = "none") + xlab("Número de Genes")


tssToGenePlots[[2]] <- ggplot(dff, aes(x = totalGenes, y = pTSS, label = label, color = label)) +
  geom_point(aes(size = size )) +
  geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8),
            vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1),
            check_overlap = TRUE, fontface = "italic") + theme_bw() +
  scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(4, 12))  +
  theme(legend.position = "none") + xlab("Número de Genes")



tssToGenePlots[[3]] <- ggplot(dff, aes(x = totalGenes, y = gTSS, label = label, color = species)) +
  geom_point(aes(size = size )) +
  geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8),
            vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1),
            check_overlap = TRUE, fontface = "italic") + theme_bw() +
  scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(4, 12)) + # Adjust point size range
  theme(legend.position = "none") + xlab("Número de Genes")

tssToGenePlots[[4]] <- get_legend(
  ggplot(dff, aes(x = totalGenes, y = gTSS, label = label, color = species)) +
    geom_point(aes(size = size )) +
    geom_text(hjust = c(.05, .05, .05, .8, .05, .05, .05, .8),
              vjust = c(-1, -1, -1, 1.55, -1, -1, -1, -1),
              check_overlap = TRUE, fontface = "italic") + theme_bw() +
    scale_size(breaks = c(1000, 5000, 10000, 20000, 50000), limits = c(1000, 60000), range = c(4, 12)) + # Adjust point size range
    guides(color = "none") + labs(color = "Espécie", size = "Reads / Gene")
)

tssToGeneGraph <- grid.arrange(
  tssToGenePlots[[3]],
  tssToGenePlots[[1]],
  tssToGenePlots[[2]],
  tssToGenePlots[[4]], 
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4)),
  # nrow = 4,
  top = "Relação entre número de genes e diferentes classes de TSS"
)


##### Better version of the one above #####

frtc_folder <- "~/files/frtc/"
tssarout_folder <- "tssarresult" # generated from running tssar.sh
readcount_folder <- "readcount" # generated from running RSubRead (altss.R)

mydf <- tibble(
  Output = as.character(),
  Depth = as.numeric(), # sum of both LibP and LibM read counts 
  filename = as.character(),
  Class = as.character(),
  Count = as.character(),
  CountsPerGene = as.numeric(),
  CountsPerGene80 = as.numeric(),
  LibP = as.character(),
  LibM = as.character(),
  DepthLibP = as.numeric(),
  DepthLibM = as.numeric(),
  cpgGenecount = as.numeric(),
  cpgPairedLibP = as.numeric(),
  cpgUnpairedLibP = as.numeric(),
  cpgPairedLibM = as.numeric(),
  cpgUnpairedLibM = as.numeric(),
  cpgGenecount80 = as.numeric(),
  cpgPairedLibP80 = as.numeric(),
  cpgUnpairedLibP80 = as.numeric(),
  cpgPairedLibM80 = as.numeric(),
  cpgUnpairedLibM80 = as.numeric()
)
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  print(paste0("Checking Sequencing Depth for ", filename))
  
  # Here, I use the same code of generating the "merged" table.
  # Because I need to extract data from TSSAR, BEFORE merging things. 
  files <- list.files(
    path = file.path(frtc_folder, filename, tssarout_folder),
    pattern = "\\.table")
  files <- files[!grepl("smear.table", files)]    # removes smear.table
  pathes <- file.path(frtc_folder, filename, tssarout_folder, files)
  
  df <- data.frame(File = files, Path = pathes) |>
    mutate(Id = sub("_([^_]+)$", "", sub("_([^_]+)$", "", File)))
  print("Found the following files | chromossomes:")
  head(df)
  
  # first, check for empty files
  invalid_df <- df |> filter(file.info(Path)$size == 0)
  if (nrow(invalid_df) > 0) {
    print("WARN There are files in tssarresult that are empty.")
    print(invalid_df$File)
    df <- df |> filter(file.info(Path)$size > 0)
  }
  
  merged <- df |>
    group_by(File) |>
    do(bind_rows(lapply(.$Path, read.table, header = FALSE, sep="\t"))) |>
    mutate(OriginalFile = sub("^[^_]*_", "", 
                              sub("^[^_]*_", "", 
                                  sub("_TSS.table", "", File)
                              ))) |>
    mutate(File = { chr <- sub("_([^_]+)$", "", sub("_([^_]+)$", "", File)) }) |>
    ungroup()
  colnames(merged) <- c("Chromossome", "Position", "Strand", "V4", "V5", "d", 
                        "Pvalue", "V8", "Class", "Location", "OriginalFile")

  # From here on out, you're no longer look at "merged" data generation.
  # Because we are not removing duplicate TSSs, instead we will make a graph
  
  # Samples
  samples <- read.csv(file.path("./samples/", paste0(filename, ".csv")))
  # samples <- samples %>% pivot_longer(cols = 1:2, values_to = "srr", names_to = "Library")
  samples$DepthLibP <- NA
  samples$DepthLibM <- NA
  
  # fetches sequencing depth from entrez. Output is in kilobase.
  getLibraryDepth <- function(libraryId) {
    csv <- entrez_fetch(db="sra", id = libraryId, rettype="runinfo", retmode="text")
    depthdf <- read.csv(text = csv)
    return(round(depthdf$bases/1000)) 
  }
  for (j in 1:nrow(samples)) {
    samples[j,]$DepthLibP <- getLibraryDepth(samples[j,]$LibP)
    samples[j,]$DepthLibM <- getLibraryDepth(samples[j,]$LibM)
  }
  samples <- samples %>% mutate(Depth = DepthLibM + DepthLibP)
  #samps <- samples %>% 
  #  pivot_wider(names_from = Library, values_from = srr) %>% 
  #  group_by(Output) %>% summarize(
  #    Depth = sum(Depth),
  #    LibP = first(na.omit(LibP)),
  #    LibM = first(na.omit(LibM))
  #  ) %>% ungroup()
  
  samples$TSS <- 0
  samples$pTSS <- 0
  for (j in 1:nrow(samples)) {
    samples[j,]$TSS <- nrow({merged |> filter(OriginalFile == samples[j,]$Output)})
    samples[j,]$pTSS <- nrow({merged |> filter(OriginalFile == samples[j,]$Output, Class == "P")})
  }

  # Lastly, I need to fetch RSubRead results
  samples$cpgPairedLibP <- 0
  samples$cpgUnpairedLibP <- 0
  samples$cpgPairedLibM <- 0
  samples$cpgUnpairedLibM <- 0
  samples$cpgGenecount <- 0
  samples$cpgPairedLibP80 <- 0
  samples$cpgUnpairedLibP80 <- 0
  samples$cpgPairedLibM80 <- 0
  samples$cpgUnpairedLibM80 <- 0
  samples$cpgGenecount80 <- 0
  samples$CountsPerGene <- 0
  samples$CountsPerGene80 <- 0
  
  unpairedfile <- file.path(start_wd, readfile_folder, paste0(filename, "-unpaired.csv"))
  pairedfile <- file.path(start_wd, readfile_folder, paste0(filename, "-paired.csv"))
  unpairedtable <- NA
  pairedtable <- NA
  if (file.exists(unpairedfile)) {
    unpairedtable <- read_csv(unpairedfile)
  }
  if (file.exists(pairedfile)) {
    pairedtable <- read_csv(pairedfile)
  }
  
  for (j in 1:nrow(samples)) {
    if (is.tibble(unpairedtable)) {
      LibPcol <- unpairedtable %>% select(starts_with(samples[j,]$LibP))
      LibMcol <- unpairedtable %>% select(starts_with(samples[j,]$LibM))
      genesCounted <- nrow(unpairedtable)
      
      # removes top10% and bottom10%, only getting the middle
      LibPEighty <- LibPcol %>% arrange(LibPcol[[1]]) %>%
        slice(ceiling(n() * 0.1):(floor(n() * 0.9)))
      LibMEighty <- LibMcol %>% arrange(LibMcol[[1]]) %>%
        slice(ceiling(n() * 0.1):(floor(n() * 0.9)))
      genesEighty <- nrow(LibPEighty)
      
      # sums everything up
      samples[j,]$cpgUnpairedLibP <- {LibPcol %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgUnpairedLibM <- {LibMcol %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgGenecount <- genesCounted
      samples[j,]$cpgUnpairedLibP80 <- {LibPEighty %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgUnpairedLibM80 <- {LibMEighty %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgGenecount80 <- genesEighty
    }
    
    # i will now copypaste the code for paired
    if (is.tibble(pairedtable)) {
      LibPcol <- pairedtable %>% select(starts_with(samples[j,]$LibP))
      LibMcol <- pairedtable %>% select(starts_with(samples[j,]$LibM))
      genesCounted <- nrow(pairedtable)
      
      # removes top10% and bottom10%, only getting the middle
      LibPEighty <- LibPcol %>% arrange(LibPcol[[1]]) %>%
        slice(ceiling(n() * 0.1):(floor(n() * 0.9)))
      LibMEighty <- LibMcol %>% arrange(LibMcol[[1]]) %>%
        slice(ceiling(n() * 0.1):(floor(n() * 0.9)))
      genesEighty <- nrow(LibPEighty)
      
      # sums everything up
      samples[j,]$cpgPairedLibP <- {LibPcol %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgPairedLibM <- {LibMcol %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgGenecount <- genesCounted
      samples[j,]$cpgPairedLibP80 <- {LibPEighty %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgPairedLibM80 <- {LibMEighty %>% summarise(across(everything(), ~ sum(., na.rm = TRUE)))}[[1]]
      samples[j,]$cpgGenecount80 <- genesEighty
    }
    
    samples[j,]$CountsPerGene <- (samples[j,]$cpgPairedLibP + samples[j,]$cpgUnpairedLibP + samples[j,]$cpgPairedLibM + samples[j,]$cpgUnpairedLibM)/samples[j,]$cpgGenecount
    samples[j,]$CountsPerGene80 <- (samples[j,]$cpgPairedLibP80 + samples[j,]$cpgUnpairedLibP80 + samples[j,]$cpgPairedLibM80 + samples[j,]$cpgUnpairedLibM80)/samples[j,]$cpgGenecount80
  }

  samples$filename <- filename
  samples <- samples %>% pivot_longer(cols = c("TSS", "pTSS"), names_to = "Class", values_to = "Count")
  
  mydf <- rbind(mydf, samples)
}

ggplot(mydf, aes(Depth, Count)) + geom_point(aes(colour = filename)) + scale_x_log10()

##### Checking how many reads that alTSS have less relative to gTSS #####
print("1. Checking how many reads that alTSS have less relative to gTSS")
# For example, if this number is 0.7, then every 
# altss has on average 70% of the readcount of their respective gtss. 
# if the number is very low (eg. 5%), then its an indicator that altss have a 
# very low expression relative to their gtss.)
plots <- list()
altssToGtss <- data.frame(
  species = as.factor(as.character()),
  altssReadsRelativeToGtss = as.numeric(),
  sd = as.numeric() # standard deviation
)

altssReadsRelativeToGtss <- list()

for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  
  data <- read_xlsx(file.path(start_wd, "altss", paste0(filename, "-alTSS.xlsx") ))
  
  # data |> select(all_of(c("SubClass", "Both_Total"))) |> group_by(SubClass) |> summarise(alignedReads = sum(Both_Total), n = n())
  
  foobar <- data |> 
    # mutate for relative alTSS intensity to its own gTSS (i.e. an alTSS with 5200 reads has 52% of the intensity of its gTSS with 10000 reads)
    group_by(Gene) |> arrange(desc(Both_Total)) |> mutate(relativ = Both_Total/max(Both_Total)) |> ungroup()
  
  altssReadsRelativeToGtss[[filename]] <- {foobar |> filter(SubClass != "gTSS")}$relativ
  
  foobar <- foobar |>
    # get avg read count
    group_by(SubClass) |> summarise(alignedReads = sum(Both_Total), avgReadRelative = mean(relativ), sd = sd(relativ), avgReadsPerGeneAbsolute = mean(Both_Total), n = n()) |> ungroup() |>
    mutate(label = paste0(SubClass, "\n(n = ", n, ")"))
  
  # foobar <- data |> group_by(SubClass) |> arrange(desc(Both_Total)) |> slice_head(n = 50) |> summarise(alignedReads = sum(Both_Total), avgReadsPerGene = mean(Both_Total), n = n(), sd = sd(Both_Total))

  # highest / lowest 10% more or less
  # data |> group_by(SubClass) |> arrange(desc(Both_Total)) |> slice_head(n = 30) |> summarise(alignedReads = sum(Both_Total), n = n(), sd = sd(Both_Total))
  # data |> group_by(SubClass) |> arrange(desc(Both_Total)) |> slice_tail(n = 100) |> summarise(alignedReads = sum(Both_Total), n = n(), sd = sd(Both_Total))
  
  altssToGtss <- altssToGtss |> add_row(
    species = filename,
    altssReadsRelativeToGtss = {foobar |> filter(SubClass == "alTSS")}$avgReadRelative,
    sd = {foobar |> filter(SubClass == "alTSS")}$sd
  )
  plots[[filename]] <- ggplot(foobar, aes(label, avgReadRelative)) +
    geom_col() + labs(title=filename) + xlab("TSS Class") + ylab("Average Reads Per Gene (Relative)")

}



df <- stack(altssReadsRelativeToGtss)
colnames(df) <- c("value", "species")

readPlot <- ggplot(df, aes(x = species, y = value)) +
  geom_boxplot() + ylim(0, 1) + ylab("alTSS reads / respective gTSS reads") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") + 
  geom_jitter(width = 0.3, size = 0.6)

ggsave(file.path(start_wd, "./misc/graphsaverageread.png"), plot = readPlot, units = "px", width=940, height=940, scale=3)



# 
# ggplot(altssToGtss, aes(x = paste0(species, "\nsd: ", round(sd, 2)), y = altssReadsRelativeToGtss)) +
#   geom_col()
# 
# readsgraphs <- grid.arrange(
#   plots[["mjannaschii"]], 
#   plots[["hsalinarum"]], 
#   plots[["pfuriosus"]], 
#   plots[["hvolcanii"]], 
#   plots[["tonnurineus"]],
#   plots[["tkodakarensis"]],
#   plots[["sacidocaldarius"]],
#   plots[["hmediterranei"]],
#   nrow = 4,
#   top = "Reads per Gene x TSS Class"
# )
# 
# ggsave(file.path(start_wd, "./misc/graphsaveragereads.png"), plot = readsgraphs, units = "px", width=940, height=940, scale=3)
# 


library(ggpubr)

##### Checking how many alTSS / gTSS exist for each p-value #####
print("2. Checking amount of alTSS for each p-value. Red is gTSS.")
pvalues = c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-15, 1e-18)
results <- list()
plotspval <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  results[[filename]] <- tibble(
    pval = numeric(),
    gtss = numeric(),
    altss = numeric()
  )
  
  data <- read_xlsx(file.path(start_wd, "altss", paste0(filename, "-alTSS.xlsx") ))
  
  for (j in 1:length(pvalues)) {
    currentp <- pvalues[j]
    g <- data |> filter(Pvalue < currentp) |>
      group_by(Gene) |> 
      mutate(TSSCount = n()) |>
      arrange(desc(Both_Total)) |>
      mutate(SubClass = ifelse(row_number() == 1, "gTSS", "alTSS")) |>
      mutate(gTSSDistance = first(Position) - Position) |>
      filter(TSSCount > 1) |>
      # now you just need to count the amount of gtss and altss
      group_by(SubClass) |> summarize(count = n())
    
    results[[filename]] <- results[[filename]] |> add_row(
      pval = currentp, 
      altss = {g |> filter(SubClass == "alTSS")}$count,
      gtss = {g |> filter(SubClass == "gTSS")}$count)
    
  }
  
  
  plotspval[[filename]] <- ggplot(results[[filename]], aes(x = -1*log(pval, base=10))) +
    geom_line(aes(x = -1*log(pval, base=10), y = gtss), color="red", size=1) + 
    geom_line(aes(x = -1*log(pval, base=10), y = altss), color="blue", size=1) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) + # increases number of ticks, see https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ylab("TSS Count") + xlab("p-value (log 1e-x)") + ggtitle( bquote(italic(.(dados[i,]$label))) ) + theme_bw()
  

}

plotspvalg <- lapply(dados$Filenames, function(name) plotspval[[name]])
pvalgraphs <- grid.arrange(
  grobs = plotspvalg,
  nrow = 4,
  top = "P-value x TSS Count."
)

leg <- get_legend(
  ggplot(data.frame(x=1:2, y=1:2, group = c("gTSS", "alTSS")), 
         aes(x, y, color = group)) +
    geom_point(size = 5) +
    scale_color_manual(values = c(gTSS = "red", alTSS = "blue"), 
                       guide = guide_legend(title = NULL)) +
    theme(legend.position = "bottom") + theme_bw()
)


p <- ggarrange(plotlist = plotspval, nrow = 4, ncol = 2,  
               common.legend = TRUE, legend = "bottom", legend.grob = leg)
annotate_figure(p, top = "P-value x TSS Count")

ggsave(file.path(start_wd, "./misc/graphspvalgraphs.png"), plot = pvalgraphs, units = "px", width=940, height=940, scale=3)

##### 5' UTR Length & TSS Count per Species #####
plotsfiveprime = list()
plotsleaderless = list()
plotstsscount = list()

for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  
  data <- read_xlsx(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx") ))
  
  ggdf <- data.frame(
    nextFeatureDistance = as.integer(sub(".*?(\\d+)nt upstream.*", "\\1", data$Location)),
    intensityTSS = data$LibP + data$LibM, 
    intensityAll = data$LibP_Total + data$LibM_Total
  ) %>% 
    mutate(utr = ifelse(nextFeatureDistance > 5, "5UTR", "Leaderless"))
  
  ggdf2 <- data.frame(
    gene_count = c(
      length(unique({data %>% filter(TSSCount == 1)}$Gene)),
      length(unique({data %>% filter(TSSCount == 2)}$Gene)),
      length(unique({data %>% filter(TSSCount == 3)}$Gene)),
      length(unique({data %>% filter(TSSCount > 3)}$Gene))
    ),
    tss_count = c("1 TSS", "2 TSSs", "3 TSSs", "4+ TSSs")
  )
  
  # ggdf$range <- cut(ggdf$nextFeatureDistance, breaks = seq(0, max(ggdf$nextFeatureDistance) + 10, by = 10))
  plotsfiveprime[[filename]] <- ggplot(ggdf[ggdf$nextFeatureDistance > 19, ], aes(nextFeatureDistance)) + 
    geom_histogram(binwidth=10) + ggtitle(paste0(filename, " (pTSS = ", nrow(data), ")")) # + scale_fill_gradientn(colors = c("white", "blue"))
  plotsleaderless[[filename]] <- ggplot(ggdf, aes(nextFeatureDistance)) + geom_bar() + xlim(-1, 21) +
    ggtitle(paste0(filename, " (pTSS = ", nrow(data), ")"))
  plotstsscount[[filename]] <- ggplot(ggdf2, aes(tss_count, gene_count)) + geom_col() + ggtitle(paste0(filename, " (pTSS = ", nrow(data), ")"))

  
  # ggplot(ggdf[ggdf$nextFeatureDistance > 5, ], aes(nextFeatureDistance, fill=intensityTSS)) + geom_histogram(binwidth=10) + scale_fill_gradientn(colors = c("white", "blue"), limits=c(0, 100000))
  # ggplot(ggdf[ggdf$nextFeatureDistance > 5, ], aes(nextFeatureDistance, fill=factor(intensityTSS))) + geom_histogram(binwidth=10) # + scale_fill_gradientn(colors = c("white", "blue"), limits=c(0, 100000))
  
}

# Arrange the plots
plotsfiveprimeg <- lapply(dados$Filenames, function(name) plotsfiveprime[[name]])
graphsfiveprime <- grid.arrange(
  grobs = plotsfiveprimeg,
  nrow = 4,
  top = "5' UTR length x Count of TSSs (With 5' UTR above 20nt)"
)


plotsleaderlessg <- lapply(dados$Filenames, function(name) plotsleaderless[[name]])
graphsleaderless <- grid.arrange(
  grobs = plotsleaderless,
  nrow = 4,
  top = "5' UTR length x Count of TSSs (With 5' UTR below 20nt)"
)

plotstsscountg <- lapply(dados$Filenames, function(name) plotstsscount[[name]])
graphstss <- grid.arrange(
  grobs = plotstsscountg,
  nrow = 4,
  top = "pTSSs with multiple TSSs (alTSSs) per species"
)

ggsave(file.path(start_wd, "./misc/graphstss.png"), plot = graphstss, units = "px", width=940, height=940, scale=3)
ggsave(file.path(start_wd, "./misc/graphsleaderless.png"), plot = graphsleaderless, units = "px", width=940, height=940, scale=3)
ggsave(file.path(start_wd, "./misc/graphsfiveprime.png"), plot = graphsfiveprime, units = "px", width=940, height=940, scale=3)


##### iTSS distribution (position relative to the start of the gene) #####
# The question for me was "what if pTSSs are being misclassified as iTSS because they are just a single nc in front of the start of a gene?"
# It's just a way of confirming that most pTSSs you caught have been classified correctly.
# Although the information of distribution of iTSSs alongside the genes could be interesting for other researchers.
itssPlots <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  gbff <- read_gbk(file.path(start_wd, genome_folder, paste0(filename, ".gbff")))
  
  tsstable <- read_excel(file.path(start_wd, "merged", paste0(filename, ".xlsx")))
  
  
  itss <- tsstable |> filter(grepl("within", Location)) |> 
    mutate(intraGene = sapply(
      strsplit(Location, "within gene\\(s\\) "),
      function(x) strsplit(x[2], ";")[[1]][1]
    ))
  
  itss <- itss |> mutate(intraGeneLocationStart = 0)
  itss <- itss |> mutate(intraGeneLocationEnd = 0)
  for (j in 1:nrow(itss)) {
    tmp <- gbff |> filter(feat_id == paste0("cds-", itss[j,]$intraGene))
    
    if (nrow(tmp) == 0) {
      print(paste0("Could not find gene information for gene ", itss[j,]$intraGene))
    } else {
      itss[j,]$intraGeneLocationStart = tmp$start[1]
      itss[j,]$intraGeneLocationEnd = tmp$end[1]
    }
  }
  
  itss <- itss |> filter(!(intraGeneLocationStart == 0 & intraGeneLocationEnd == 0))
  
  itss <- itss |> mutate(intraGeneDistance = 0)
  for (j in 1:nrow(itss)) {
    if (itss[j,]$Strand == "+") {
      itss[j,]$intraGeneDistance <- abs(itss[j,]$Position - itss[j,]$intraGeneLocationStart)
    } else {
      itss[j,]$intraGeneDistance <- abs(itss[j,]$intraGeneLocationEnd - itss[j,]$Position)
    }
  }
  
  itss <- itss |> mutate(intraGeneSize = abs(intraGeneLocationEnd - intraGeneLocationStart))
  
  itss <- itss |> mutate(intraGeneLocationRelative = abs(intraGeneDistance/intraGeneSize))
  
  itssPlots[[filename]] <- ggplot(itss, aes(x = intraGeneLocationRelative)) + 
    geom_histogram(binwidth = 0.01) + xlim(0, 1) + ggtitle(filename)
}


# Arrange the plots
itssPlotsg <- lapply(dados$Filenames, function(name) itssPlots[[name]])
itssGraph <- grid.arrange(
  grobs = itssPlotsg,
  nrow = 4,
  top = "iTSS position relative to start of the gene"
)

ggsave(file.path(start_wd, "./misc/graphsitsspos.png"), plot = itssGraph, units = "px", width=940, height=940, scale=3)


##### Distance each alTSS is from its respective gTSS #####
# I don't know what to say for expected results of this one.
plotsb <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  tsstable <- read_excel(file.path(start_wd, "altss", paste0(filename, "-alTSS.xlsx")))
  data.frame(contagens_unpaired$counts) %>%
    summarize(across(everything(), sum), total_sum = sum(across(everything(), sum)))
  altss <- tsstable |> filter(SubClass == "alTSS")
  
  plotsb[[filename]] <- ggplot(altss, aes(x = gTSSDistance)) + 
    geom_histogram(binwidth = 10) + ggtitle( bquote(italic(.(dados[i,]$label))) ) + 
    xlim(-250, 250) + theme_bw() + xlab("Distância (nt)") + ylab("alTSS")
  
}

# Arrange the plots
plotsbg <- lapply(dados$Filenames, function(name) plotsb[[name]])
graphsb <- grid.arrange(
  grobs = plotsbg,
  nrow = 4,
  top = "Distância entre alTSS para o gTSS de um mesmo gene"
)

ggsave(file.path(start_wd, "./misc/graphsaltssdist.png"), plot = graphsb, units = "px", width=940, height=940, scale=3)




##### Checking proportionality between (alTSS / pTSS) and (5' UTR with >20 nt / pTSS)

# The purpose of this graph is to check if there's an association between the amount of alTSSs/gTSSs with long 5' UTRs (longer than 20nc)
# The reasoning is that one of the proposed purposes of alternative transcription is to generate mRNAs with different sized 5' UTRs
# So if most of your TSSs don't have a significantly sized 5' UTR


df <- data.frame(
  tssCount = as.numeric(),
  ptssCount = as.numeric(),
  altssCount = as.numeric(),
  tssUtr0 = as.numeric(),
  tssUtr5 = as.numeric(),
  tssUtr10 = as.numeric(),
  tssUtr15 = as.numeric(),
  tssUtr20 = as.numeric(),
  tssUtr25 = as.numeric(),
  tssUtr30 = as.numeric(),
  tssUtr40 = as.numeric(),
  tssUtr50 = as.numeric(),
  species = as.factor(as.character())
)
# plotsc <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  # gbff <- read_gbk(file.path(start_wd, genome_folder, paste0(filename, ".gbff")))
  
  everytss <- read_excel(file.path(start_wd, "merged", paste0(filename, ".xlsx")))
  ptss <- read_excel(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx")))
  
  ptss <- ptss |> mutate(UtrLength = as.numeric(gsub("nt upstream.*", "", Location)))
  
  df <- df |> add_row(
    tssCount = nrow(everytss),
    ptssCount = nrow(ptss),
    altssCount = nrow(ptss |> filter(TSSCount > 1)),
    tssUtr0 = nrow(ptss |> filter(UtrLength > 0)),
    tssUtr5 = nrow(ptss |> filter(UtrLength > 5)),
    tssUtr10 = nrow(ptss |> filter(UtrLength > 10)),
    tssUtr15 = nrow(ptss |> filter(UtrLength > 15)),
    tssUtr20 = nrow(ptss |> filter(UtrLength > 20)),
    tssUtr25 = nrow(ptss |> filter(UtrLength > 25)),
    tssUtr30 = nrow(ptss |> filter(UtrLength > 30)),
    tssUtr40 = nrow(ptss |> filter(UtrLength > 40)),
    tssUtr50 = nrow(ptss |> filter(UtrLength > 50)),
    species = filename
  )
  
  # plotsc[[filename]] <- ggplot(ptss, aes(UtrLength)) + geom_histogram(binwidth = 10) + ggtitle(paste0(filename))
}

# Arrange the plots
# plotscf <- lapply(dados$Filenames, function(name) plotsc[[name]])
# grid.arrange(
#   grobs = plotscf,
#   nrow = 4,
#   top = "P-value x TSS Count"
# )

df <- df %>% mutate(alternativeRatio = altssCount / ptssCount,
                    utrRatio = tssUtr20 / ptssCount)


ggplot(df, aes(x = species, y = ptssCount)) +
  geom_col()

model <- lm(utrRatio ~ alternativeRatio, data = df)
summary(model)

utraltssratiorelative <- ggplot(df, aes(x = alternativeRatio, y = utrRatio, color = species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  ylab("TSSs with >20nc 5' UTRs / pTSSs") + xlab("(alTSSs + gTSSs) / pTSSs") + 
  ggtitle(paste0("Ratios, R^2 = ", summary(model)$adj.r.squared))




model <- lm(tssUtr20 ~ altssCount, data = df)
summary(model)

utraltssratioabsolute <- ggplot(df, aes(x = altssCount, y = tssUtr20, color = species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  ylab("TSSs with >20nc 5' UTRs (Absolute Value)") + xlab("alTSSs + gTSSs count") + 
  ggtitle(paste0("Ratios, R^2 = ", summary(model)$adj.r.squared))

ggsave(file.path(start_wd, "./misc/graphsaltssratioabs.png"), plot = utraltssratioabsolute, units = "px", width=940, height=940, scale=3)
ggsave(file.path(start_wd, "./misc/graphsaltssratiorel.png"), plot = utraltssratiorelative, units = "px", width=940, height=940, scale=3)




# Select columns that have featuredistance > 20 ggdf[ggdf$nextFeatureDistance > 19, ]





#tmp

resultsdf <- tibble(
  species = as.factor(as.character()),
  tss = as.numeric(),
  ptss = as.numeric(),
  gtss = as.numeric(),
  altss = as.numeric()
)
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  dfeverytss <- read_excel(file.path(start_wd, "merged", paste0(filename, ".xlsx")))
  dfptss <- read_excel(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx")))
  
  resultsdf <- resultsdf |> add_row(
    species = filename,
    tss = nrow(dfeverytss),
    ptss = nrow(dfptss),
    gtss = nrow({dfptss |> filter(SubClass == "gTSS")}),
    altss = nrow({dfptss |> filter(SubClass == "alTSS")})
  )
}


plotsidk <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  data <- read_xlsx(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx") ))
  
  ggdf <- data.frame(
    nextFeatureDistance = as.integer(sub(".*?(\\d+)nt upstream.*", "\\1", data$Location)),
    intensityTSS = data$LibP + data$LibM, 
    intensityAll = data$LibP_Total + data$LibM_Total
  ) %>% 
    mutate(utr = ifelse(nextFeatureDistance > 5, "5UTR", "Leaderless"))
  
  ggdf2 <- data.frame(
    gene_count = c(
      length(unique({data %>% filter(TSSCount == 1)}$Gene)),
      length(unique({data %>% filter(TSSCount == 2)}$Gene)),
      length(unique({data %>% filter(TSSCount == 3)}$Gene)),
      length(unique({data %>% filter(TSSCount > 3)}$Gene))
    ),
    tss_count = c("1 TSS", "2 TSSs", "3 TSSs", "4+ TSSs")
  )
  
  mytitle <- paste0(dados[i, ]$Shortened, " (pTSS: ", nrow(data), ")")
                        
  # ggdf$range <- cut(ggdf$nextFeatureDistance, breaks = seq(0, max(ggdf$nextFeatureDistance) + 10, by = 10))
  plotsidk[[filename]] <- ggplot(ggdf, aes(nextFeatureDistance)) + 
    geom_histogram(binwidth=10) + 
    ggtitle(LETTERS[i], subtitle = bquote(italic(.(mytitle)))) + # ggtitle(paste0(mytitle, "\n(pTSS = ", nrow(data), ")")) +
    theme_set(theme_bw(base_size = 12, base_family = "Helvetica")) +
    xlim(-12, 250) + xlab("Tamanho 5' UTR") +
    ylab("Contagem pTSS") 
}

# Arrange the plots
plotsidkg <- lapply(dados$Filenames, function(name) plotsidk[[name]])
graphsb <- grid.arrange(
  grobs = plotsidkg,
  nrow = 4 # , top = "Title"
)

ggsave(file.path(start_wd, "./misc/graphsutrfull.png"), plot = graphsb, units = "px", width=940, height=940, scale=3)












# Using dplyr, select only columns "..." and "..." from iris
# Solution: https://tidyselect.r-lib.org/reference/all_of.html

# Using dplyr, get the average petal length for each species. Return a 2x3 tibble
# Solution: iris |> group_by(Species) |> summarise(petalTotal = sum(Petal.Length), n = n())

# Using dplyr, get the 5 biggest petals for each species. Return a tibble.
# Solution: iris |> group_by(Species) |> arrange(desc(Petal.Length)) |> slice_head(n = 5)

# Using dplyr, create a new column called "MaxPetalLength" which is specific for each species.
# iris |> group_by(Species) |> arrange(desc(Petal.Length)) |> mutate(test = max(Petal.Length))

# https://www.njtierney.com/post/2022/08/08/fun-across/
# Using dplyr, multiply the value of all columns that start with "Sepal" by 100
# iris %>% mutate(across(.cols = starts_with("Sepal"), .fns = function(x) x*100))

# Using dplyr, group by species. Then, get the mean and standard deviation for all columns that have numeric values
# iris |> group_by(Species) |> summarize(across(.cols = where(is.numeric), list(mean=mean, sd=sd) ))

# # First, you need to think: Will I collapse this function (like tunrning 150 lines -> 3 lines)? Yes! Then I'll use summarize.
# # When summarizing, you could specify the name of each column, like the one below:
# #  iris |> group_by(Species) |> summarize(FirstMean = mean(Sepal.Length))
# # But instead, you want to use the across(.cols). 
# #   Across applies a function for all rows, and by specifying the ".cols" parameter.. You get the idea.

# Using dplyr, filter only for numeric values. Then, SUM ALL NUMBERS, EVERYTHING
# iris |> select(where(is.numeric)) |> summarize(across(everything(), sum), total_sum = sum(across(everything(), sum)))

# data.frame(contagens_unpaired$counts) %>% summarize(across(everything(), sum), total_sum = sum(across(everything(), sum)))








