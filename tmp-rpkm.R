
library(tidyverse)

# same arg used in altss.R
frtc_f <- path.expand(file.path("~/files/frtc"))
dados <- read_tsv("./data/dados.tsv")

for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  print(paste0("Running for ", filename))
  
  ptss <- read_tsv(paste0("./results/altss/", filename, "-pTSS.tsv"))
  
  gff <- tibble(as.data.frame(
    readGFF(paste0("./data/genomes/", filename, ".gff"), 
            tags=c("ID", "Parent", "Name"))
  ))
  gff$Length <- abs(gff$end - gff$start)
  
  geneLength <- gff %>% 
    mutate(Parent = sub("gene-", "", Parent)) %>%
    select(Parent, Length) %>%
    deframe()

  system2("bash", args = c("./tmp-rpkmfetch.sh", frtc_f, filename), 
          stdout = "./tmp-rpkmnfetch.log",
          stderr = "./tmp-rpkmnfetch.log")
  
  samples <- read_csv(paste0("./data/samples/", filename, ".csv"))
  samples <- samples %>% pivot_longer(1:2, values_to = "library")
  
  mrfn <- paste0("./results/readcount/", filename, "-mappedreads.tsv")
  if (!file.exists(mrfn)) {
    print(paste0("Could not find mapped reads for ", filename))
    next;
  }
  mr <- read_tsv(mrfn)
  
  mappedreads <- mr %>%
    mutate(library = sub("-.*", "", sample)) %>%
    left_join(select(samples, library, OriginalFile = Output), by = "library") %>%
    group_by(OriginalFile) %>% 
    summarize(sum(mapped)) %>% 
    rename(TotalMappedReads = "sum(mapped)") %>%
    drop_na() %>% 
    deframe()

  print(mappedreads)
  ##### HRPM #####

  headreads <- read_tsv(paste0("./results/readcount/", filename,"-window50.bed"), col_names = FALSE)
  colnames(headreads) <- c("Chromossome", "PosStart", "PosEnd", "Strand", "OriginalFile", "Headreads")
  headreads <- headreads %>% 
    mutate(mapreads = mappedreads[OriginalFile], 
           hrpm = 1e6*Headreads/mapreads,
           Position = ifelse(Strand == "+", PosStart, PosEnd)) %>%
    select(-PosStart, -PosEnd)
  
  ptss <- ptss %>% 
    left_join(headreads, by = c("Chromossome", "Position", "Strand", "OriginalFile"), suffix = c(".old", "")) %>%
    select(-ends_with(".old"))
  
  ##### RPKM #####
  # small fix first: let's merge -paired and -unpaired.csv
  unpaired_readcount_fn <- paste0("./results/readcount/", filename, "-unpaired.csv")
  paired_readcount_fn <- paste0("./results/readcount/", filename, "-paired.csv")
  if (file.exists(unpaired_readcount_fn) & file.exists(paired_readcount_fn)) {
    unpr_readcount <- read_csv(unpaired_readcount_fn)
    pear_readcount <- read_csv(paired_readcount_fn)
    merged_readcount <- cbind(unpr_readcount[1], 
                              unpr_readcount[-1] + pear_readcount[-1])
  } else {
    merged_readcount <- read_csv(unpaired_readcount_fn)
  }
  
  libraryName <- samples %>% select(library, Output) %>% deframe()
  # rownames are genes
  merged_readcount <- as.data.frame(merged_readcount)
  rownames(merged_readcount) <- merged_readcount[[1]]
  merged_readcount <- merged_readcount[,-1]
  # filter for only valid library names (discards all samples not listed in samples.csv)
  validLibraryName <- sub("-.*", "", colnames(merged_readcount)) %in% names(libraryName)
  merged_readcount <- merged_readcount[, validLibraryName]
  # replaces colname by repname, then sums
  names(merged_readcount) <- libraryName[sub("-.*", "", names(merged_readcount))]
  merged_readcount <- t(rowsum(t(merged_readcount), group = sub("\\..*", "", colnames(merged_readcount))))
  # makes column for genes
  merged_readcount <- as.data.frame(merged_readcount)
  merged_readcount$Gene <- sub("gene-", "", rownames(merged_readcount))
  merged_readcount <- as_tibble(merged_readcount)
  # merged_readcount
  head(merged_readcount)
  
  ptss$Genereads <- NA
  for (j in 1:nrow(ptss)) {
    matchedGene <- {merged_readcount %>% filter(Gene == ptss[j,]$Gene)}
    if (nrow(matchedGene) > 0) {
      ptss[j,]$Genereads <- matchedGene[[ptss[j,]$OriginalFile]]
    }
  }
  
  ptss <- ptss %>% mutate(rpkm = 1e9 * Genereads / (mapreads * geneLength[Gene]))
  
  ptss <- write_tsv(ptss, paste0("./results/altss/", filename, "-pTSS.tsv"))
  
}
