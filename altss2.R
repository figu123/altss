library(tidyverse)    # Core used for working with datatables
library(rtracklayer)  # Reading genomic annotation files
library(seqinr)       # Extracting upstream regions (promoter regions)
library(Rsubread)     # Used in getting read counts for each gene.


# Used for 5' UTR graphs
# library(ggplot2)


# first, run fetch-genomes.sh

# --------------------- #
##### CONFIGURATION #####

start_wd <- "~/files/analises/altss"
setwd(start_wd)

# how many threads is Rsubread allowed to use. it can run on 1.
rsubread_nthreads <- 1

# How much upstream/downstream of each pTSS is extracted for promoter analysis.
# (For motif searching. The length of the sequence below will be appended to the final spreadsheet.)
upstream_count <- 100
downstream_count <- 10


# minimal DISTANCE in base pairs between alTSSs.
# A distance of 5 is the equivalent of a 10 nucleotide window (5 nt left, 5 nt right)
threshold_mindist <- 5
threshold_pvalue <- 1e-6

# INPUT DATA #

# 1. Assemblies Table (dados.tsv) 
# "dados" is table with 3 columns: Assembly, Filenames and Species
# Where the first line is the header with the name of each column.
#
# Assembly: String. ID of the NCBI RefSeq assembly from NCBI datasets.
# Filenames: String. Will be used for the pathing of your files.
# Species: String. An identifier string. Can be anything you want.
#
### Example
# Assembly	Filenames	Species
# GCF_000006805.1	hsalinarum	Halobacterium salinarum NRC-1
###
dados <- read_tsv("./data/dados.tsv")

# 2. FRTC folder
# It is expected for it to have one folder for each species, each one containing
# the subfolder "coverage" from running frtc.sh, and the subfolder "tssarresult" 
# from running the scripts in ./tssar.
#   https://github.com/alanlorenzetti/frtc
### Example
#  "${frtc_f}/hsalinarum/coverage"
#  "${frtc_f}/hsalinarum/tssarresult"
# must be a valid folders, if you've defined "hsalinarum" in the "Filenames"
# of the dados.tsv table (see variable above)
###
frtc_f <- file.path("~/files/frtc")
# file.path(frtc_f, "hsalinarum", "tssarresult")

# 3. Samples folder
# Inside this folder, there should be .csv files with information about which 
# runs are TEX+, which runs are TEX-, and a string identifier.
# The runs should be the SRR ID (same as the .sra downloaded from NCBI)
# The identifier can be anything you want, but make them unique from each other.
# The first column should be formatted as: LibP,LibM,Output
### Example: (samples_folder_path)/hsalinarum.csv
# LibP,LibM,Output
# SRR6953855,SRR6953869,ref1a
# SRR6953856,SRR6953870,ref1b
# SRR6953857,SRR6953871,ref2a
# SRR6953858,SRR6953872,ref2b
# SRR6953859,SRR6953865,17ha
# SRR6953860,SRR6953866,17hb
# SRR6953861,SRR6953867,37ha
# SRR6953862,SRR6953868,37hb
# SRR6953863,SRR6953873,86ha
# SRR6953864,SRR6953874,86hb
###
samples_f <- file.path("./data/samples")

# 4. Genomes folder
# Should have all genomic data. Given a valid Assemblies Table (dados.tsv),
# you can generate all the contents of this folder by running fetch-genomes.sh
genome_f <- file.path("./data/genomes")

# --------------------- #
##### MAIN #####

merge_tssar_results <- function(filename) {
  # Unifies the results of TSSAR into a single table.
  # This unified result DOES INCLUDE duplicates 
  #   (eg. a TSS that is found in one sample will be found in another)
  # The removal of duplicates happens in the next step.
  files <- list.files(
    path = file.path(frtc_f, filename, "tssarresult"),
    pattern = "\\.table")
  files <- files[!grepl("smear.table", files)]    # removes smear.table
  pathes <- file.path(frtc_f, filename, "tssarresult", files)
  
  df <- tibble(data.frame(File = files, Path = pathes) |>
                 mutate(Id = sub("_([^_]+)$", "", sub("_([^_]+)$", "", File))))
  print("Found the following files | chromossomes:")
  head(df)
  # TODO: this is messy. there's a lot of subs to manipulate the filename.
  
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
  colnames(merged) <- c("Chromossome", "Position", "Strand", "TSSid", "V5", "d", 
                        "Pvalue", "V8", "Class", "Location", "OriginalFile")
  
  # the TSSid column is outdated, since you merged multiple TSSAR results, multiple
  # ids could exist. so these should be reassigned right now.
  merged <- merged |>
    arrange(Chromossome, Position) |>
    mutate(TSSid = paste0("TSS_", sprintf("%07d", row_number())))
  
  dir.create(file.path("./results/merged"), showWarnings = FALSE)
  write_tsv(merged, file.path("./results/merged/", paste0(filename, ".tsv")))
}

gene_readcount <- function(filename) {
  # Amount of reads with coverage in each gene is determined by the Rsubread package.
  # It will read all the .bam files found inside the coverage folder.
  # (Read about github/frtc to understand more about the coverage folder).
  
  # -- Listing Files -- #
  files_paired <- list.files(path = file.path(frtc_f, filename, "coverage"),
                             pattern = "-paired\\.bam$")
  if (length(files_paired > 0)) {
    files_paired <- paste0(file.path(frtc_f, filename, "coverage"), "/", files_paired)
  }
  files_unpaired <- list.files(path = file.path(frtc_f, filename, "coverage"),
                               pattern = "-unpaired\\.bam$")
  files_unpaired <- paste0(file.path(frtc_f, filename, "coverage"), "/", files_unpaired)
  # -- Counting -- #
  contagens_unpaired <- featureCounts(
    files = files_unpaired,
    annot.ext = file.path(genome_f, paste0(filename, "-original.gff1")),
    isGTFAnnotationFile = TRUE,    # vale para GTF e GFF como usamos GFF so trocar
    GTF.featureType = "gene",      # tipo de feição no GFF.
    GTF.attrType = "ID",           # nome do atributo do gene no GFF. Sempre olhe se o gff tem o "ID" do gene alguns nao tem (a grande maioria tem)
    isPairedEnd = FALSE,            # coloque FALSE se for single-end, sempre verificar
    nthreads = rsubread_nthreads
  )
  
  head(contagens_unpaired$counts) #Olhe se deu certo SEMPRE 
  
  dir.create(file.path("./results/readcount"), showWarnings = FALSE)
  # criar um arquivo com as contagens
  write.csv(contagens_unpaired$counts, 
            file = file.path("./results/readcount", paste0(filename, "-unpaired.csv")))
  
  if (length(files_paired) > 0) {
    contagens_paired <- featureCounts(
      files = files_paired,
      annot.ext = file.path(genome_f, paste0(filename, "-original.gff1")),
      isGTFAnnotationFile = TRUE,    # vale para GTF e GFF como usamos GFF so trocar
      GTF.featureType = "gene",      # tipo de feição no GFF.
      GTF.attrType = "ID",           # nome do atributo do gene no GFF. Sempre olhe se o gff tem o "ID" do gene alguns nao tem (a grande maioria tem)
      isPairedEnd = TRUE,             # coloque FALSE se for single-end, sempre verificar
      nthreads = rsubread_nthreads
    )
    
    head(contagens_paired$counts) #Olhe se deu certo SEMPRE 
    
    # criar um arquivo com as contagens
    
    write.csv(contagens_paired$counts, 
              file = file.path("./results/readcount", paste0(filename, "-paired.csv")))
  }
}

merge_altss <- function(df) {
  # This is like a greedy signal-intensity (d column) clustering
  # Essentially, you define a window, and alTSSs within that window are "merged"
  # -
  # You should run this function following this https://dplyr.tidyverse.org/reference/group_map.html
  # .x is the data frame for one gene.
  if (nrow(df) == 1) {
    df$dMerged <- df$d
    return(df)
  }
  
  # "df" is grouped by gene. So expect a short list of TSSs for the same gene.
  # sorts by "d" column, and sequentially (from highest d, to lowest d)
  # merges all TSSs within that nucleotide window.
  df <- df %>% arrange(desc(d)) 
  
  used <- rep(FALSE, nrow(df)) 
  # "used" TSSs either have 
  #   (1) already merged all alTSSs within their window (with the for loop below)
  #   (2) already been merged by another TSSs (that's why used[close] reassigns)
  output <- rep(0, nrow(df)) 
  
  for (j in 1:nrow(df)) {
    if (used[j]) next
    
    center <- df$Position[j]
    close <- which(!used & abs(df$Position - center) < threshold_mindist)
    
    output[j] <- sum(df$d[close])
    used[close] <- TRUE
  }
  
  df$dMerged <- output
  return(df)
}

extract_sequence_from_position <- function(PositionStart, PositionEnd, Chromossome, Strand, genomes) {
  # Useful for extract promoter regions of genes.
  # Requires the variables:
  #   upstream_count -- How many bases upstream to extract from the starting position of the gene?
  #   downstream_count -- How many bases downstream to extract from the starting position of the gene?
  #   genomes -- A list containing the sequence of each chromossome.
  # Genes that are in the positive strand begin at "PositionStart".
  # Genes that are in the negative strand begin at "PositionEnd", and the result is the reverse of the complementary
  Chromossome <- sub("\\..*", "", Chromossome)
  # print(paste("Accessing", chr_name))
  fasta_seq <- genomes[[Chromossome]]
  if (is.null(fasta_seq)) {
    print(paste0("Failed to access genomes[[",Chromossome,"]]"))
    return(FALSE)
  }
  # Defines the start/end position of the sequence, based on if it's an antisense or sense
  if (Strand == "+") {
    # For a gene with start_position 65274
    # getFrag(genomes[[1]], begin=65274, end=65276)
    # [1] "a" "t" "g"
    start_pos <- as.numeric(PositionStart) - upstream_count
    end_pos <- start_pos + upstream_count + downstream_count
  } else {
    # For a gene with end_position 30378
    # rev(comp(getFrag(genomes[[1]], begin=30376, end=30378)))
    # [1] "a" "t" "g"
    # Without rev(comp), the result is [1] "c" "a" "t"
    end_pos <- as.numeric(PositionEnd) + upstream_count
    start_pos <- end_pos - upstream_count - downstream_count
  }
  # clamps to min/max value of the sequence # 
  # print(paste0("endpos", toString(end_pos), ",,,,", "fasta_seq_len", toString(length(fasta_seq))))
  if (is.na(end_pos) || is.na(start_pos)) {
    print(paste("Error: Trying to fetch position for", toString(PositionStart), ",", toString(PositionEnd), " for chr ", toString(Chromossome)," resulted in NAs, somehow."))
    return(FALSE)
  }
  if (end_pos > length(fasta_seq)) {
    print(paste("Warning: Fragment starting at",toString(start_pos),"is at the end of the genome, resulting in a sequence smaller than the one requested."))
    end_pos = length(fasta_seq)
  }
  if (start_pos < 0) {
    print(paste("Warning: Fragment starting at",toString(start_pos),"is at the beginning of the genome, resulting in a sequence smaller than the one requested."))
    start_pos = 0
  }
  # Uses seqinr to extract the sequences
  if (Strand == "+") {
    frag <- paste0(getFrag(fasta_seq, begin=start_pos, end=end_pos),collapse="")
  } else {
    frag <- paste0(rev(comp(getFrag(fasta_seq, begin=start_pos, end=end_pos))),collapse="")
  }
  return(frag)
}

signal_counting <- function(merged, grs, samples) {
  # Information of signal intensity from the .bedgraph files is read, and then
  # saved into the "merged" table. So this function defines five columns:
  #  LibP : Signal intensity of the TEX+ Library of this TSS.
  #  LibM : Signal intensity of the TEX- Library of this TSS.
  #  LibP_Total : Signal intensity of all TEX+ libraries in the position of this TSS.
  #  LibM_Total : Same as above, for TEX- libraries.
  #  Both_Total : Signal intensity of all libraries in the position of this TSS.
  ##
  # Example: Suppose you have the following samples:
  #   LibP       LibM           Output
  # 1 SRR1951998 SRR1951997     b1
  # 2 SRR1952000 SRR1951999     b2
  # 3 SRR1952002 SRR1952001     b3
  # If you have a TSS that is found from the difference of libraries SRR1952000 
  # and SRR1951999, then the intensity of these libraries will be saved in LibP
  # and LibM columsn, respectively.
  # The intensity of all LibP values in that TSS position is saved in LibP_Total.
  # etc.
  ##
  # Input:
  #   merged: Table with all TSSs.
  #   grs: List with GenomicRanges of all BedGraphs.
  #   samples: Table that specifies which samples are LibP, which ones are LibM.
  ##
  # Output:
  #   merged_new: Table with all TSSs and columns with signal intensity of each TSS.
  ##
  # In the case of paired files, there will exist 3 files in $filename/coverage:
  #  Run-paired-fwd.bedgraph
  #  Run-unpaired-fwd.bedgraph
  #  Run-fwd.bedgraph <--- This is the only one that will be considered for counting.
  # The other two aren't considered to avoid redundancy.
  
  
  # defining columns
  merged$LibP <- 0
  merged$LibP_Total <- 0
  merged$LibM <- 0
  merged$LibM_Total <- 0
  merged$Both_Total <- 0
  
  # for each line in the merged table...
  for (j in 1:nrow(merged)) {
    # merged[j,]$LibP <- samples[samples$Output == merged[j,]$OriginalFile, ]$LibP
    
    # .. filters only for the libraries that have the same strand direction
    foo <- ifelse(merged[j,]$Strand == "+", "fwd", "rev")
    valid_names <- names(grs)[grepl(foo, names(grs))]
    
    # .. finds the LibP and LibM that was used to determine the TSS of this line
    LibP <- samples[samples$Output == merged[j,]$OriginalFile, ]$LibP
    LibM <- samples[samples$Output == merged[j,]$OriginalFile, ]$LibM
    
    # .. for each library that is in the same strand as this one...
    for (k in 1:length(valid_names)) {
      # .. determines the score for that TSS location for that library
      score <- grs[[valid_names[k]]][start(grs[[valid_names[k]]]) <= merged[j,]$Position & 
                                       end(grs[[valid_names[k]]]) >= merged[j,]$Position &
                                       grepl(merged[j,]$Chromossome, as.character(seqnames(grs[[valid_names[k]]]))),
                                     "score"]$score
      
      # .. sets that score to column LibP or LibM, if we're analysing the
      # libP or libM library that was used to spot this TSS in particular
      if (grepl(LibP, valid_names[k])) {
        merged[j,]$LibP <- score
      } else if (grepl(LibM, valid_names[k])) {
        merged[j,]$LibM <- score
      }
      
      # .. otherwise just split it on whether its a positive or a neg library
      is_positive_library <- any(grepl(sub("-.*", "", valid_names[k]), samples$LibP) == TRUE)
      if (is_positive_library) {
        merged[j,]$LibP_Total <- merged[j,]$LibP_Total+score
      } else {
        merged[j,]$LibM_Total <- merged[j,]$LibM_Total+score
      }
      
      merged[j,]$Both_Total <- merged[j,]$Both_Total+score
    }
  }
  
  # finished redefining the columns
  return(merged)
}

generate_gtf <- function(df, outfile_name) {
  # converts table in altss to gtf file that can be imported in IGV for visualization.
  
  # assumes that these colnames exist:
  # colnames(df) <- c("Chromossome", "Position", "Strand", "V4", "V5", "d", "Pvalue", "V8", "Class", "Location")
  # "Chromossome" "Position" "Strand" "V4" "V5" "d" "Pvalue" "V8" "Class" "Location" "OriginalFile" 
  # "LibP" "LibP_Total" "LibM" "LibM_Total" "Both_Total" "Gene" "TSSCount" "SubClass" "gTSSDistance" "Sequence" "Product"
  
  hsv <- df
  hsv$gtfinfo <- ""
  
  if (nrow(hsv) == 0) {
    return(FALSE)
  }
  for (j in 1:nrow(hsv)) {
    if (any(colnames(hsv) == "SubClass")) {
      # filtered df -> gtf
      hsv[j,]$gtfinfo <- paste0('position \"', toString(hsv[j,]$Position), '\";',
                                'tss_class \"', hsv[j,]$Class, '\";',
                                'tss_subclass \"', hsv[j,]$SubClass, '\";',
                                'tss_location \"', gsub(" ", "_", toString(hsv[j,]$Location)), '\";',
                                'tss_strand \"', toString(hsv[j,]$Strand), '\";',
                                'p_value \"', toString(hsv[j,]$Pvalue), '\";',
                                'run_id \"', toString(hsv[j,]$OriginalFile), '\";',
                                'product \"', toString(hsv[j,]$Product), '\";')
    } else {
      # merged df -> gtf
      hsv[j,]$gtfinfo <- paste0('position \"', toString(hsv[j,]$Position), '\";',
                                'tss_class \"', hsv[j,]$Class, '\";',
                                'tss_location \"', gsub(" ", "_", toString(hsv[j,]$Location)), '\";',
                                'tss_strand \"', toString(hsv[j,]$Strand), '\";',
                                'p_value \"', toString(hsv[j,]$Pvalue), '\";',
                                'run_id \"', toString(hsv[j,]$OriginalFile), '\";')
    }
  }
  
  hsv2 <- data.frame(
    V1 = hsv$Chromossome,
    V2 = "TSS",
    V3 = "TSS",
    V4 = hsv$Position,
    V5 = hsv$Position+1,
    V6 = ".",
    V7 = hsv$Strand,
    V8 = ".",
    V9 = hsv$gtfinfo
  )
  
  write_tsv(hsv2, paste0(outfile_name, ".gtf"), col_names = FALSE, escape = "none")
}

generate_memeinput <- function(filtered, output_f) {
  # converts table in altss to .fasta files that can run in MEME for mapping motifs
  
  if (nrow(filtered) == 0) {
    return(FALSE)
  }
  for (j in 1:nrow(filtered)) {
    header <- paste0(">", filtered[j,]$Chromossome, 
                     "|", filtered[j,]$Strand,
                     "|", filtered[j,]$Position,
                     "|", filtered[j,]$d,
                     "|", filtered[j,]$Class)
    sequence <- toupper(filtered[j,]$Sequence)
    
    write(header, 
          file=file.path("./", output_f, paste0(filename, "-pTSS-all.fasta")), 
          append=TRUE)
    write(sequence, 
          file=file.path("./", output_f, paste0(filename, "-pTSS-all.fasta")), 
          append=TRUE)
    
    if (filtered[j,]$TSSCount > 1) {
      write(header, 
            file=file.path("./", output_f, paste0(filename, "-pTSS-has-alTSS.fasta")), 
            append=TRUE)
      write(sequence, 
            file=file.path("./", output_f, paste0(filename, "-pTSS-has-alTSS.fasta")), 
            append=TRUE)
    } 
    if (filtered[j,]$SubClass == "gTSS") {
      write(header, 
            file=file.path("./", output_f, paste0(filename, "-pTSS-uniq.fasta")), 
            append=TRUE)
      write(sequence, 
            file=file.path("./", output_f, paste0(filename, "-pTSS-uniq.fasta")), 
            append=TRUE)
    }
    
    
  }
}


##### RUNNING #####

sink("altss2.log", append = TRUE)

print(paste0("Running at ", getwd()))
dir.create(file.path("./results"), showWarnings = FALSE)

for (i in 1:nrow(dados)) {
  tryCatch(
    {
      
      
      filename <- dados[i,]$Filenames
      gcf_id <- dados[i,]$Assembly
      
      print("---------------------------------------")
      print(paste0("Running for species: ", filename))
      
      ##### checking input files #####
      print("1. Checking input files...")
      if (!dir.exists(file.path(frtc_f, filename))) {
        paste0("Error: Folder not found at ", file.path(frtc_f, filename),". Skipping.")
        next
      }
      if (!dir.exists(file.path(frtc_f, filename, "tssarresult"))) {
        print("Error: The FRTC folder doesn't have a tssarresult folder. Skipping.")
        print("Make sure to setup & run tssar.sh first. Check the README.")
        next
      }
      if (!file.exists( file.path(samples_f, paste0(filename, ".csv")) )) {
        print("Error: There's no samples.csv file for this species. Skipping.")
        print("For formatting on how to make this file, please check the README.")
        next
      }
      if (!file.exists( file.path(genome_f, paste0(filename, ".gff")) )) {
        print("Error: Missing genome files. Download them by running fetch-genomes.sh.")
        next
      }
      
      # Assigning it to genomes
      print("1. Checking genomes...")
      genomes <- list()
      chromossomes <- list.files(
        path = file.path(start_wd, genome_f),
        pattern = "\\.fasta")
      # filters for current species
      chromossomes <- chromossomes[grepl(filename, chromossomes)]
      pathes <- file.path(start_wd, genome_f, chromossomes)
      chromossomes <- sub(paste0(filename,"-"), "", sub(".fasta", "", chromossomes))
      for (j in 1:length(chromossomes)) {
        genomes[[sub("\\..*", "", chromossomes[j])]] <- read.fasta(pathes[[j]])[[1]]
      }
      print("1.1. Chromossomes:")
      print(chromossomes)
      print(summary(genomes))
      
      gff_data <- tibble(as.data.frame(
        readGFF(file.path(start_wd, genome_f, paste0(filename, ".gff")), 
                tags=c("ID", "Parent", "Name"))
      ))
      
      print("1.2. GFF_Data")
      print(head(gff_data))
      
      features <- read_tsv(file.path(start_wd, genome_f, paste0(filename, "-features.tsv")))
      names(features)[1] <- "feature"
      
      print("1.3. Features")
      print(head(features))
      
      gaf_data <- read_tsv(file.path(start_wd, genome_f, paste0(filename, ".gaf")),
                           comment = "!", col_names = FALSE)
      print("1.4. Gene Ontology")
      print(head(gaf_data))
      
      # # This is deprecated. It was only used in the signal counting function.
      # print("1.2. Checking bedgraphs")
      # files <- list.files(path = file.path(frtc_f, filename, "coverage"), 
      #                     pattern = "\\.bedgraph")
      # files <- files[!grepl(".*paired-.*", files)]
      # 
      # grs <- List()
      # for (j in 1:length(files)) {
      #   gr <- import(BEDGraphFile(
      #     file.path(frtc_f, filename, "coverage", files[j])
      #   ))
      #   grs[[files[j]]] <- gr
      # }
      
      ##### processing tables #####
      print("2. Merging TSSAR tables")
      merge_tssar_results(filename)
      
      print("3. Processing merged table pTSSs. ")
      ptss <- read_tsv(file.path("./results/merged/", paste0(filename, ".tsv")))
      ptss <- ptss %>% filter(Class == "P")
      
      print("3.0. Fixing Chromossome Names")
      # i don't know *when this started happening*, but some chromsosomes names are truncated.
      # since we have already loaded the Chromossome names from the genome annotations above,
      # it makes sense to replace the truncated chromossomes with fixed ones.
      ptss$Chromossome <- sapply(ptss$Chromossome, function(x) { chromossomes[startsWith(chromossomes, x)] })
      
      print("3.1. Removing duplicate entries")
      # if a TSS is in the same chromossome, same position, same strand
      # then we assume its the same TSS.
      ptss <- ptss %>%
        group_by(Chromossome, Position, Strand) %>%
        slice_max(d, n = 1, with_ties = FALSE) %>% 
        ungroup()
      
      print("3.2. Creating column with Gene, Protein IDs, 5' UTR Length")
      ptss$Gene <- NA
      ptss$UTRLength <- NA
      ptss <- ptss %>% 
        mutate(
          Gene = ifelse(str_detect(Location, "upstream of gene"), 
                        sub(";.*$", "", sub("^.*upstream of gene ", "", Location)), NA),
          UTRLength = ifelse(str_detect(Location, "upstream of gene"), 
                             as.numeric(sub(";.*$", "", sub("nt upstream of gene .*", "", Location))), NA)
        )
      
      ptss$Protein <- NA
      ptss$Product <- NA
      for (j in 1:nrow(features)) {
        cds <- features %>% 
          filter(feature == "CDS") %>% 
          filter(str_detect(locus_tag, ptss[j,]$Gene))
        
        if (nrow(cds) > 0) {
          ptss[j,]$Product = cds[1,]$name
          ptss[j,]$Protein = cds[1,]$product_accession
        }
      }
      
      print("3.3. Creating column with upstream/downstream sequence")
      ptss <- ptss %>%
        mutate(Sequence = mapply(extract_sequence_from_position, 
                                 Position, Position, Chromossome, Strand, 
                                 MoreArgs = list(genomes = genomes)))
      
      print("3.4. Columns TSSCount, SubClass, gTSSDistance")
      # Assigns gTSS for one with lowest value in "d" (equivalent of p-value)
      ptss <- ptss %>%
        group_by(Gene) %>%
        mutate(TSSCount = n()) %>%
        # here, if you do
        # filter(n()>1) |>  # it filters only for alTSSs
        # gTSS = Genuine TSS, alTSS = Alternative TSS.
        arrange(desc(d)) %>%
        mutate(SubClass = ifelse(row_number() == 1, "gTSS", "alTSS")) %>%
        # The distance is from the TSS relative to the gTSS. For example, an alTSS
        # with distance -50, indicates that this alTSS is 50bp upstream of the gTSS.
        mutate(gTSSDistance = first(Position) - Position) %>%
        ungroup() %>%
        arrange(TSSid)
      
      paste0("  Saving results to altss/", filename,"-unfiltered-pTSS.tsv")
      dir.create(file.path("./results/altss"), showWarnings = FALSE)
      write_tsv(ptss, paste0("./results/altss/", filename, "-unfiltered-pTSS.tsv"))
      # These unfiltered-pTSSs are useful for analysing information (see graphs.R)
      #   by analysing how information is before applying filters, you determine
      #   what is the best p-value, or merging window (for the filters below)
      
      
      print(paste0("3.5. Applying filters: p-value and minimum alTSS distance"))
      print(paste0("   p-value threshold: ", threshold_pvalue))
      print(paste0("   merging pTSS within a distance of: ", threshold_mindist, "nt from each other"))
      
      ptss <- ptss %>% 
        # pvalue
        filter(Pvalue < threshold_pvalue) %>% 
        # merging by distance
        group_by(Gene) %>%
        group_modify(~ merge_altss(.x)) %>%
        ungroup() %>%
        filter(dMerged != 0) %>%
        # now that this is done, you need to reapply step 3.4., since the 
        # TSS count per gene will have channged
        # Assigns gTSS for one with lowest value in "d" (equivalent of p-value)
        group_by(Gene) %>%
        mutate(TSSCount = n()) %>%
        # here, if you do
        # filter(n()>1) |>  # it filters only for alTSSs
        # gTSS = Genuine TSS, alTSS = Alternative TSS.
        arrange(desc(d)) %>%
        mutate(SubClass = ifelse(row_number() == 1, "gTSS", "alTSS")) %>%
        # The distance is from the TSS relative to the gTSS. For example, an alTSS
        # with distance -50, indicates that this alTSS is 50bp upstream of the gTSS.
        mutate(gTSSDistance = first(Position) - Position) %>%
        ungroup() %>%
        arrange(TSSid)
      
      paste0("  Saving results to altss/", filename,"-pTSS.tsv")
      write_tsv(ptss, paste0("./results/altss/", filename, "-pTSS.tsv"))
      
      ##### analysing signal intensity #####
      print("4. Reading information of signal intensity")
      
      print("4.1. Using RSubRead to get read counts for each gene")
      gene_readcount(filename)
      
      # print("4.2. Finding intensity for each TSS (This may take a while)")
      # Deprecated.
      # The intensity of each TSS could be used to define gTSS/alTSS,
      #   but because gene-expression is regulated based on an organisms conditions, 
      #   it's not like this would give me any comparable data.
      # For example: I can define the alTSS/gTSS based on intensity for
      #   Halobacterium salinarum's dRNA-Seq experiment, which was run in
      #   different points of the curve under optimal conditions.
      #   And I can define the alTSS/gTSS based on intensity for dRNA-Seq of
      #   Halobacterium mediterranei's dRNA-Seq experiment, which was run in
      #   different mediums.
      # That still doesn't mean the alTSS of hsalinarum is comparable to hmediterranei.
      #   Since those intensities were only defined by the conditions these organisms were in.
      #   And these conditions aren't comparable, since they were designed for different studies.
      # So really, it doesn't contribute for my analysis to define who is gTSS who is alTSS.
      # What contributes is seeing the number of alTSSs and comparing orthogroup-ed genes between species.
      #  
      # samples <- read.csv(file.path(samples_f, paste0(filename, ".csv")))
      # ptss <- signal_counting(ptss, grs, samples)
      # paste0("  Saving results to altss/", filename,"-pTSS.tsv")
      # write_tsv(ptss, paste0("./results/altss/", filename, "-pTSS.tsv"))
      
      ##### generating input/summary files #####
      print("5. Generating input files")
      dir.create(file.path("./results/memeinput"), showWarnings = FALSE)
      generate_memeinput(ptss, file.path("./results/memeinput/"))
      
      generate_gtf(ptss, file.path("./results/altss/", paste0(filename, "-pTSS")))
      
      print("6. Writing summary.tsv")
      if (file.exists("./results/summary.tsv")) {
        summarydf <- read_tsv("./results/summary.tsv")
      } else {
        summarydf <- tibble(
          DateRun = as.character(),
          Filenames = as.character(),
          Genes_total = as.numeric(),
          Genes_ptss = as.numeric(),
          Genes_multiple_ptss = as.numeric(),
          Transcripts_5utr = as.numeric(),
          Transcripts_leaderless = as.numeric()
        )
      }
      summarydf <- summarydf %>% add_row(
        DateRun = paste0("TIME", as.character(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
        Filenames = filename,
        Genes_total = nrow(gff_data %>% filter(type == "CDS")),
        Genes_ptss = length(unique(ptss$Gene)),
        Genes_multiple_ptss = length(unique({ptss %>% filter(TSSCount > 1)}$Gene)),
        Transcripts_5utr = nrow(ptss %>% filter(UTRLength > 10)),
        Transcripts_leaderless = nrow(ptss %>% filter(UTRLength <= 10))
      )
      write_tsv(summarydf, "./results/summary.tsv")
      
      
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
    }
  )
}

print("Done.")

sink()




