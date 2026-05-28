
start_wd <- "~/files/filesaguia-orthofinder/orthofinder/analise"


resultdf <- data.frame(
  species <- as.factor(),
  ptss <- int()
)
for (i in 1:nrow(dados)) {
  filename <- dados$Filenames[i]
  assembly <- dados$Assembly[i]
  print(paste("For file", filename))
  
  altss <- read_excel(file.path(start_wd, "altss", paste0(filename, "-alTSS.xlsx")))
  altssuq <- altss |> group_by(Gene) |> slice_head() |> ungroup()
  ptssuq <- read_excel(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx")))
  ptssuq <- ptssuq |> group_by(Gene) |> slice_head() |> ungroup()
  
  tot_ptss <- nrow(ptssuq)
  hasortho_ptss <- nrow({ptssuq |> filter(!is.na(Orthogroup))})
  tot <- nrow(altssuq)
  hasortho <- nrow({altssuq |> filter(!is.na(Orthogroup))})
  
  print(paste0(hasortho_ptss, "/", tot_ptss, " --- ", hasortho, "/", tot))
  
  
}
