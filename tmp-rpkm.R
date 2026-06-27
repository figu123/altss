

library(tidyverse)

# same arg used in altss.R
frtc_f <- path.expand(file.path("~/files/frtc"))
dados <- read_tsv("./data/dados.tsv")

for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  
  ptss <- read_tsv(paste0("./results/altss/", filename, "-pTSS.tsv"))
  
  gff <- tibble(as.data.frame(
    readGFF(paste0("./data/genomes/", filename, ".gff"), 
            tags=c("ID", "Parent", "Name"))
  ))
  gff$Length <- abs(gff$end - gff$start)
  
  ##### HRPM #####
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
  
  headreads <- read_tsv("./results/readcount/mjannaschii-window50.bed", col_names = FALSE)
  colnames(headreads) <- c("Chromossome", "PosStart", "PosEnd", "Strand", "OriginalFile", "Headreads")
  headreads <- headreads %>% 
    mutate(mapreads = mappedreads[OriginalFile], 
           hrpm = 1e6*Headreads/mapreads,
           Position = ifelse(Strand == "+", PosStart, PosEnd)) %>%
    select(-PosStart, -PosEnd)
  
  # samtools flagstat file.bam
  ptss <- ptss %>% left_join(headreads, by = c("Chromossome", "Position", "Strand", "OriginalFile"))
  
  
}
