library(rstudioapi)
library(rtracklayer)
library(seqinr)

# Setup Global Variables #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
start_wd <- getwd()

# Gff Info
files <- list.files(pattern = "\\.gff")
for (i in 1:length(files)) {
  print(paste("Reading",files[i]))
  gff_data <- as.data.frame(rtracklayer::readGFF(files[i], tags=c("ID", "Parent", "Name")))
  annotated_genes <- gff_data |> filter(grepl("^gene-", ID)) |> 
    group_by(ID) |> slice_head(n=1) |> ungroup()
  print(paste("Number of genes:",toString(nrow(annotated_genes))))
  size <- gff_data |> filter(grepl("^region", type)) |> summarize(total = sum(end)) |> pull(total)
  print(paste("Genome Size:",toString(size)))
}

# Fasta Info
files <- list.files(pattern = "\\.fasta")
for (i in 1:length(files)) {
  print(paste("Reading",files[i]))
  fasta_seq <- seqinr::read.fasta(files[[i]])[[1]]
  print(paste("Nucleotide Count:",length(fasta_seq)))
}
