
library(readxl)

dados <- read.table("./dados.tsv", sep = "\t", header = TRUE)

for (i in 1:nrow(dados)) {
  df <- read_excel(paste0("./altss/", dados[i,]$Filenames, "-pTSS.xlsx"))

  df$gtfinfo <- ""
  for (j in 1:nrow(df)) {
    df[j,]$gtfinfo <- paste0('position \"', toString(df[j,]$Position), '\";',
                              'tss_class \"', df[j,]$Class, '\";',
                              'tss_subclass \"', df[j,]$SubClass, '\";',
                              'tss_location \"', gsub(" ", "_", toString(df[j,]$Location)), '\";',
                              'tss_strand \"', toString(df[j,]$Strand), '\";',
                              'p_value \"', toString(df[j,]$Pvalue), '\";',
                              'run_id \"', toString(df[j,]$OriginalFile), '\";',
                              'product \"', toString(df[j,]$Product), '\";')
  }
  
  
  hsv2 <- data.frame(
    V1 = df$Chromossome,
    V2 = "TSS",
    V3 = "TSS",
    V4 = df$Position,
    V5 = df$Position+1,
    V6 = ".",
    V7 = df$Strand,
    V8 = ".",
    V9 = df$gtfinfo)

  outfile_name <- paste0("./altss/", dados[i,]$Filenames, "-pTSS")
  write_tsv(hsv2, paste0(outfile_name, ".gtf"), col_names = FALSE, escape = "none")
}
