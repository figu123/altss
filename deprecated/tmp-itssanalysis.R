# i had a thought:

# how far inside are iTSS inside the gene? 
# because if it's like in the very first base pairs, who's to say its not a pTSS flagged as iTSS?

start_wd <-"~/files/analises/full"
setwd(start_wd)
genome_folder <- "./genomes"

dados <- read_tsv("./dados.tsv")

library(gridExtra)

plots <- list()
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
  
  plots[[filename]] <- ggplot(itss, aes(x = intraGeneLocationRelative)) + 
    geom_histogram(binwidth = 0.01) + xlim(0, 1) + ggtitle(filename)
}


grid.arrange(
  plots[["mjannaschii"]], 
  plots[["hsalinarum"]], 
  plots[["pfuriosus"]], 
  plots[["hvolcanii"]], 
  plots[["tonnurineus"]],
  plots[["tkodakarensis"]],
  plots[["sacidocaldarius"]],
  plots[["hmediterranei"]],
  nrow = 4,
  top = "P-value x TSS Count"
)


plotsb <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  tsstable <- read_excel(file.path(start_wd, "altss", paste0(filename, "-alTSS.xlsx")))
  
  altss <- tsstable |> filter(SubClass == "alTSS")
  
  plotsb[[filename]] <- ggplot(altss, aes(x = gTSSDistance)) + 
    geom_histogram(binwidth = 10) + ggtitle(filename)
  
}

grid.arrange(
  plotsb[["mjannaschii"]], 
  plotsb[["hsalinarum"]], 
  plotsb[["pfuriosus"]], 
  plotsb[["hvolcanii"]], 
  plotsb[["tonnurineus"]],
  plotsb[["tkodakarensis"]],
  plotsb[["sacidocaldarius"]],
  plotsb[["hmediterranei"]],
  nrow = 4,
  top = "P-value x TSS Count"
)








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
  species = as.character()
)
plotsc <- list()
for (i in 1:nrow(dados)) {
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  # gbff <- read_gbk(file.path(start_wd, genome_folder, paste0(filename, ".gbff")))
  
  alltss <- read_excel(file.path(start_wd, "merged", paste0(filename, ".xlsx")))
  ptss <- read_excel(file.path(start_wd, "altss", paste0(filename, "-pTSS.xlsx")))
  
  ptss <- ptss |> mutate(UtrLength = as.numeric(gsub("nt upstream.*", "", Location)))

  df <- df |> add_row(
    tssCount = nrow(alltss),
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
  
  plotsc[[filename]] <- ggplot(ptss, aes(UtrLength)) + geom_histogram(binwidth = 10) + 
    ggtitle(paste0(filename))
}

# Arrange the plots
plotscf <- lapply(dados$Filenames, function(name) plotsc[[name]])
grid.arrange(
  grobs = plotscf,
  nrow = 4,
  top = "P-value x TSS Count"
)

# Trying to get some relationships

df$species <- as.factor(df$species)

# pairwise plot (all numeric columns)

df <- df %>% mutate(alternativeRatio = altssCount / ptssCount,
                    utrRatio = tssUtr20 / ptssCount)


ggplot(df, aes(x = species, y = ptssCount)) +
  geom_col()

ggplot(df, aes(x = utrRatio, y = alternativeRatio, color = species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  ylab("alternativeRatio") + xlab("utrRatio") + ggtitle("Ratios")

model <- lm(alternativeRatio ~ utrRatio, data = df)
summary(model)



ggplot(df, aes(x = altssCount, y = tssUtr20, color = species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  ylab("alternativeRatio") + xlab("utrRatio") + ggtitle("Ratios")

model <- lm(tssUtr20 ~ altssCount, data = df)
summary(model)

