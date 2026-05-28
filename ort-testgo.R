library(UniProt.ws)

# reference
# columns(up)
# keytypes(up)

# Protein accession
acc2 <- c("WP_010903180.1", "WP_004060682.1")

up <- UniProt.ws()

res <- select(
  up,
  keys = acc2,
  keytype = "RefSeq_Protein",
  columns = c("protein_name", "cc_function", "go_id", "go_c", "go_f", "go_p",
              "ec", "cc_subcellular_location",
              "protein_existence", "length")
)

res



start_wd <- "~/files/filesaguia-orthofinder/orthofinder/analise"
setwd(start_wd)


filesuffix <- "pTSS" # altss / ptss


for (i in 1:nrow(dados)) {
  filename <- dados$Filenames[i]
  assembly <- dados$Assembly[i]
  print(paste("For file", filename))
  
  
  gbff <- read_gbk(file.path(start_wd, genome_folder, paste0(filename, ".gbff")))
  gbff <- gbff |> filter(type == "CDS")
  
  altss <- read_excel(file.path(start_wd, "altss", paste0(filename, "-", filesuffix, ".xlsx")))
  # allowed_genes <- altss$Gene
  
  altss$go_function <- NA
  for (j in 1:nrow(altss)) {
    currentgene <- altss[j,]$Gene
    aa <- gbff |> filter(grepl(currentgene, feat_id))
    if (nrow(aa) > 0) {
      altss[j,]$go_function <- aa[1,]$go_function
    }
  }
  
  df2 <- altss %>% filter(!is.na(go_function)) %>% 
    dplyr::count(go_function) %>%
    arrange(desc(n))
  
  df2 <- head(df2, 10)
  
  ggplot(df2, aes(x=reorder(go_function, n), y=n)) +
    geom_point(size=8) +
    coord_flip() + x_label()
  
}

