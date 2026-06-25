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

# Used for parsing gbff file
library(gggenomes)

# Setup:
# Set your current folder name. If using rstudio, it can be done programatically.
# library(rstudioapi)
# start_wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd("~/files/analises/full")
start_wd <- getwd()
print(paste0("Running at ", getwd()))

# dados table.
# "dados" is table with 3 columns: Assembly, Species and Filenames.
# The first line show have these colnames, and the following lines should have
# Assembly: String. ID of the NCBI RefSeq assembly from NCBI datasets.
# Species: String. An identifier string. Can be anything you want.
# Filenames: String. Will be used for the pathing of your files.
### Example
# Assembly	Species	Filenames
# GCF_000006805.1	Halobacterium salinarum NRC-1	hsalinarum
###
dados <- read_tsv("./dados.tsv")

# frtc folder.
# Inside this folder, the following should exist:
#  frtc/$filename/tssarinput     --- Contains .bam files used as input for TSSAR tool.
#  frtc/$filename/tssarresult    --- Generated after running tssar.sh.
#  frtc/$filename/coverage       --- Contains coverage .bam files.
# Where "$filename" is defined by the 'dados.tsv' table of your species, eg. "hsalinarum"
frtc_folder <- "~/files/frtc/"

# samples folder.
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
samples_folder <- "~/files/samples/"

# Will attempet to use the AssemblyID specified in dados to download the genomes.
# If you prefer, you can provide the following genomes in the ./genomes folder:
# $filename.gff               --- Genes annotation file
# $filename.fna               --- Full genome
# $filename-NucCoreID.fasta   --- .fasta files of
download_genomes <- TRUE

# How much upstream/downstream of each pTSS is extracted for promoter analysis.
# (For motif searching. The length of the sequence below will be appended to the final spreadsheet.)
upstream_count <- 100
downstream_count <- 10

### 
# hardcoded variables
###

# folders found in this script
genome_folder <- "genomes"
merged_folder <- "merged"
filtered_folder <- "altss"
memeinput_folder <- "memeinput"
misc_folder <- file.path(start_wd, "misc")
readcount_folder <- file.path(start_wd, "readcount")
dir.create(readcount_folder, showWarnings = FALSE)
dir.create(misc_folder, showWarnings = FALSE)


# folders found inside the frtc folder
coverage_folder <- "coverage" # "coverage" or "coverage_uniq"
# Coverage files that will be used in readcount / gene and calculating TSS intensity.
# warning: fiveprime, threeprime and TSSAR input will not be uniq
# When quantifying TSS intensity (such as in defining alternative TSS) it is interesting
# to use "coverage" rather than "coverage_uniq", because
tssarout_folder <- "tssarresult" # generated from running tssar.sh

# Download all genomic information needed for the scripts. (.gff, .fna, .fasta)
# Uses the GCF_ID provided in 'dados' in order to fetch these files.
fetch_genome <- function(gcf_id, filename) {
  # Searches in the entrez database for this GCF_ID.
  summary <- entrez_search(db = "assembly", term = gcf_id, retmax = 1)
  if (length(summary$ids) == 0) { 
    print(paste0("ERR: Could not find assembly: ", gcf_id))
    return(FALSE) 
  }
  # Gets summary
  doc <- entrez_summary(db = "assembly", id = summary$ids[[1]])
  ftp <- doc$ftppath_refseq
  if (is.null(ftp)) { 
    print(paste0("ERR: Could not find ftp for assembly: ", gcf_id))
    return(FALSE)
  }
  # 1. Fetching GFF (List of Genes)
  
  file_url <- paste0(ftp, "/", basename(ftp), "_genomic.gff.gz")
  download.file(file_url, 
                destfile = file.path(start_wd, genome_folder, paste0(filename,".gff.gz")), 
                mode = "wb")
  R.utils::gunzip(file.path(start_wd, genome_folder, paste0("./",filename,".gff.gz")), overwrite = TRUE)
  # 2. Fetching FNA (Visualize in IGV)
  file_url <- paste0(ftp, "/", basename(ftp), "_genomic.fna.gz")
  download.file(file_url, 
                destfile = file.path(start_wd, genome_folder, paste0(filename,".fna.gz")), 
                mode = "wb")
  R.utils::gunzip(file.path(start_wd, genome_folder, paste0("./",filename,".fna.gz")), overwrite = TRUE)
  # 3. Fetching GBFF (Genomic Information)
  file_url <- paste0(ftp, "/", basename(ftp), "_genomic.gbff.gz")
  download.file(file_url, 
                destfile = file.path(start_wd, genome_folder, paste0(filename,".gbff.gz")), 
                mode = "wb")
  R.utils::gunzip(file.path(start_wd, genome_folder, paste0("./",filename,".gbff.gz")), overwrite = TRUE)
  # 4. Fetching AssemblyReport, in order to find the chromossomes.
  file_url <- paste0(ftp, "/", basename(ftp), "_assembly_report.txt")
  download.file(file_url, 
                destfile = file.path(start_wd, genome_folder, paste0(filename,"-asreport.txt")), 
                mode = "wb")
  # lazy way to extract the nucleotide IDs
  report <- readLines(file.path(start_wd, genome_folder, paste0(filename,"-asreport.txt")))
  report <- report[grepl("^[^#]", report)]
  writeLines(report, file.path(start_wd, genome_folder, paste0(filename,"-asreport2.tmp.txt")))
  report <- read.table(file.path(start_wd, genome_folder, paste0(filename,"-asreport2.tmp.txt")), header = FALSE, sep = "\t")
  # return(TRUE)
  for (j in 1:nrow(report)) {
    nucid <- report[j,]$`V7`
    print(paste("Fetching fasta for",nucid))
    summary <- entrez_search(db = "nuccore", term = nucid, retmax = 1, use_history = TRUE)
    if (length(summary$ids) == 0) { 
      print(paste0("Error: Could not fetch genome ",nucid))
      return(FALSE)
    }
    fasta_sequence <- entrez_fetch(db = "nuccore", web_history = summary$web_history, rettype = "fasta")
    write(fasta_sequence, 
          file = file.path(start_wd, genome_folder, paste0(filename,"-",nucid,".fasta")))
    
  }
}
extract_sequence_from_position <- function(PositionStart, PositionEnd, Chromossome, Strand) {
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

# Amount of reads with coverage in each gene is determined by the Rsubread package.
# It will read all the .bam files found inside the coverage folder.
# (Read about github/frtc to understand more about the coverage folder).
readcount <- function(filename) {
  #### Listing Files ####
  files_paired <- list.files(path = file.path(frtc_folder, filename,  coverage_folder),
                             pattern = "-paired\\.bam$")
  if (length(files_paired > 0)) {
    files_paired <- paste0(file.path(frtc_folder, filename, coverage_folder), "/", files_paired)
  }
  files_unpaired <- list.files(path = file.path(frtc_folder, filename, coverage_folder),
                               pattern = "-unpaired\\.bam$")
  files_unpaired <- paste0(file.path(frtc_folder, filename, coverage_folder), "/", files_unpaired)
  #### Counting ####
  contagens_unpaired <- featureCounts(
    files = files_unpaired,
    annot.ext = file.path(start_wd, genome_folder, paste0(filename, ".gff")),
    isGTFAnnotationFile = TRUE,    # vale para GTF e GFF como usamos GFF so trocar
    GTF.featureType = "gene",      # tipo de feição no GFF.
    GTF.attrType = "ID",           # nome do atributo do gene no GFF. Sempre olhe se o gff tem o "ID" do gene alguns nao tem (a grande maioria tem)
    isPairedEnd = FALSE            # coloque FALSE se for single-end, sempre verificar
  )
  
  head(contagens_unpaired$counts) #Olhe se deu certo SEMPRE 
  
  # criar um arquivo com as contagens
  write.csv(contagens_unpaired$counts, 
            file = file.path(readcount_folder, paste0(filename, "-unpaired.csv")))
  
  if (length(files_paired) > 0) {
    contagens_paired <- featureCounts(
      files = files_paired,
      annot.ext = file.path(start_wd, genome_folder, paste0(filename, ".gff")),
      isGTFAnnotationFile = TRUE,    # vale para GTF e GFF como usamos GFF so trocar
      GTF.featureType = "gene",      # tipo de feição no GFF.
      GTF.attrType = "ID",           # nome do atributo do gene no GFF. Sempre olhe se o gff tem o "ID" do gene alguns nao tem (a grande maioria tem)
      isPairedEnd = TRUE             # coloque FALSE se for single-end, sempre verificar
    )
    
    head(contagens_paired$counts) #Olhe se deu certo SEMPRE 
    
    # criar um arquivo com as contagens
    
    write.csv(contagens_paired$counts, 
              file = file.path(readcount_folder, paste0(filename, "-paired.csv")))
  }
}

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

TSS_signal_counting <- function(merged, grs, samples) {
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
      score <- grs[[valid_names[k]]][seqnames(grs[[valid_names[k]]]) == merged[j,]$Chromossome & 
                                       start(grs[[valid_names[k]]]) <= merged[j,]$Position & 
                                       end(grs[[valid_names[k]]]) >= merged[j,]$Position, 
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

# converts into gtf file that can be imported in IGV for visualization.
generate_gtf <- function(df, outfile_name) {
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


generate_memeinput <- function(filtered) {
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
          file=file.path(start_wd, memeinput_folder, paste0(filename, "-pTSS-all.fasta")), 
          append=TRUE)
    write(sequence, 
          file=file.path(start_wd, memeinput_folder, paste0(filename, "-pTSS-all.fasta")), 
          append=TRUE)
    
    if (filtered[j,]$TSSCount == 1) {
      write(header, 
            file=file.path(start_wd, memeinput_folder, paste0(filename, "-pTSS-uniq.fasta")), 
            append=TRUE)
      write(sequence, 
            file=file.path(start_wd, memeinput_folder, paste0(filename, "-pTSS-uniq.fasta")), 
            append=TRUE)
    } else {
      if (filtered[j,]$SubClass == "gTSS") {
        
        write(header, 
              file=file.path(start_wd, memeinput_folder, paste0(filename, "-gTSS.fasta")), 
              append=TRUE)
        write(sequence, 
              file=file.path(start_wd, memeinput_folder, paste0(filename, "-gTSS.fasta")), 
              append=TRUE)
        
      } else {
        
        write(header, 
              file=file.path(start_wd, memeinput_folder, paste0(filename, "-alTSS.fasta")), 
              append=TRUE)
        write(sequence, 
              file=file.path(start_wd, memeinput_folder, paste0(filename, "-alTSS.fasta")), 
              append=TRUE)
        
      }
    }
  }
  
}

print(dados)
results <- list()
genomes <- list()

for (i in 1:nrow(dados)) {
  tryCatch(
    {
      
  filename <- dados[i,]$Filenames
  gcf_id <- dados[i,]$Assembly
  
  print("---------------------------------------")
  print(paste0("Running for species: ", filename))

  ### TSSAR ###
  print("Checking for TSSAR/Sample files...")
  if (!dir.exists(file.path(frtc_folder, filename))) {
    print("Error: There's no frtc folder for this species. Skipping.")
    next
  }
  if (!dir.exists(file.path(frtc_folder, filename, tssarout_folder))) {
    print("Error: There's no tssar results for this species. Skipping. Make sure to setup & run tssar.sh first.")
    next
  }
  ### Samples ###
  if (!file.exists( file.path(samples_folder, paste0(filename, ".csv")) )) {
    print("Error: There's no samples.csv file for this species. Skipping.")
    print("For formatting on how to make this file, please check information in samples folder variable.")
    next
  }
  
  #if (!dir.exists(file.path(frtc_folder, filename, "tssarresult"))) {
  #  print("2.1. Running TSSAR")
  #  # Copy everything from ./tssar to frtc_folder/tssarinput
  #  files <- list.files("/home/jfigueiredo/files/teste/rtracklayer/tssar", full.names = TRUE)
  #  file.copy(files, file.path(frtc_folder, filename, "tssarinput"))
  #  # Then runs the tss_tudo.sh script, which runs tssar alongside other scripts.
  #  system(paste0("bash ", 
  #                file.path(frtc_folder, filename, "tssarinput"), 
  #                "/TSS_tudo.sh >> ", 
  #                file.path(frtc_folder, filename, "tssarinput"), 
  #                "/log-tssar.log"))
  #}
  
  
  print("1. Setting up Genomes")
  #### Genomes ####
  print("1.1. Fetching Genomes")
  if (download_genomes) {
    fetch_genome(gcf_id, filename)
  }
  
  # Part 1.2. Reads genome data, saves into the following variables:
  #   chromossomes  -- c() with chromossomes IDs
  #   genomes       -- list() with sequence of these chromossomes
  #   gff_data      -- datatable with GFF information
  print("1.2. Reading Genomes")
  chromossomes <- list.files(
    path = file.path(start_wd, genome_folder),
    pattern = "\\.fasta")
  # filters for current species
  chromossomes <- chromossomes[grepl(filename, chromossomes)]
  pathes <- file.path(start_wd, genome_folder, chromossomes)
  chromossomes <- sub(paste0(filename,"-"), "", sub(".fasta", "", chromossomes))
  genomes <- list()
  
  for (j in 1:length(chromossomes)) {
    genomes[[sub("\\..*", "", chromossomes[j])]] <- read.fasta(pathes[[j]])[[1]]
  }
  print("Read the following Chromossomes:")
  print(chromossomes)
  print(summary(genomes))
  
  gff_data <- as.data.frame(readGFF(file.path(start_wd, genome_folder, paste0(filename,".gff")), 
                                    tags=c("ID", "Parent", "Name")))
  print("GFF_Data")
  print(head(gff_data))
  
  #### Merging ####
  # Merges the TSSAR result files
  #   merged      -- dataframe
  print("2. Merging TSSAR tables")
  files <- list.files(
    path = file.path(frtc_folder, filename, tssarout_folder),
    pattern = "\\.table")
  files <- files[!grepl("smear.table", files)]    # removes smear.table
  pathes <- file.path(frtc_folder, filename, tssarout_folder, files)
  
  df <- data.frame(File = files, Path = pathes) |>
    mutate(Id = sub("_([^_]+)$", "", sub("_([^_]+)$", "", File)))
  print("Found the following files | chromossomes:")
  head(df)
  
  # TODO: this is messy. there's a lot of subs to manipulate the filename.
  # im sure there's an easier way to do this. but i dont want to do a regexp.
  
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
    mutate(File = { chr <- sub("_([^_]+)$", "", sub("_([^_]+)$", "", File))
                    chromossomes[grepl(chr, chromossomes)]
                  }) |>
    ungroup()
  colnames(merged) <- c("Chromossome", "Position", "Strand", "V4", "V5", "d", 
                        "Pvalue", "V8", "Class", "Location", "OriginalFile")
  
  # removing duplicate entries.
  #   if a TSS is in the same chromossome, same position, same strand
  #   then we assume its the same TSS.
  merged <- merged |>
    group_by(Chromossome, Position, Strand) |>
    slice_max(d, n = 1, with_ties = FALSE) |> ungroup()
  
  write_xlsx(merged, file.path(start_wd, merged_folder, paste0(filename, ".xlsx")))
  
  
  #### Coverage ####
  # Loads the GR for your species
  #    grs      -- list() containing genomic ranges of each run
  print("3. Loading coverage information")
  
  print("3.1. Checking gene counts")
  readcount(filename)
  
  print("3.2. Reading bedgraph files...")
  files <- list.files(path = file.path(frtc_folder, filename, coverage_folder), 
                      pattern = "\\.bedgraph")
  files <- files[!grepl(".*paired-.*", files)]
  
  grs <- List()
  for (j in 1:length(files)) {
    gr <- import(BEDGraphFile(
      file.path(frtc_folder, filename, coverage_folder, files[j])
    ))
    grs[[files[j]]] <- gr
  }
  
  print("3.3. Finding intensity for each TSS (This may take a while)")
  samples <- read.csv(file.path(samples_folder, paste0(filename, ".csv")))
  merged <- TSS_signal_counting(merged, grs, samples)
  write_xlsx(merged, file.path(start_wd, merged_folder, paste0(filename, ".xlsx")))
  generate_gtf(merged, file.path(start_wd, merged_folder, paste0(filename, "-TSS")))
  
  print("3.4. Creating column with Gene IDs only of pTSSs")
  merged$Gene <- NA
  merged <- merged %>% 
    mutate(
      Gene = ifelse(str_detect(Location, "upstream of gene"), 
                    sub(";.*$", "", sub("^.*upstream of gene ", "", Location)), NA)
    )
    
  print("3.5. Reading gene products for pTSSs (proteins from the .gbff)")
  # To annotate information about the genes, we will be reading the .gbff of this species.
  gbff <- read_gbk(file.path(start_wd, genome_folder, paste0(filename, ".gbff")))
  gbff_cds <- gbff |> filter(type == "CDS")
  merged$Product <- NA
  if (nrow(merged) == 0) {
    print("WARN WARNING: Merged table is empty...?")
  } else {
    for (j in 1:nrow(merged)) {
      if (!str_detect(merged[j,]$Class, "P")) {
        # ignores if there's no primary class.
        next
      }
      cds <- gbff_cds %>% 
        filter(str_detect(feat_id, merged[j,]$Gene))
      
      if (nrow(cds) > 0) {
        merged[j,]$Product = cds[1,]$product
      } else {
        print(paste0("WARN: Could not find product for gene ", merged[j,]$Gene))
      }
    }
  }
  write_xlsx(merged, file.path(start_wd, merged_folder, paste0(filename, ".xlsx")))
  
  print("4. Filtering for multiple pTSSs per gene and extracting upstream sequence (for promoter region identification)")
  filtered <- merged %>%
    # Filtering for pTSS and p-value
    filter(Class == "P") %>% 
    # filter(Pvalue <= pvalue_filter) %>%
    # Makes a column with "Genes", and further analysis will be grouped by them.
    # update: i moved this into the step prior to this one.
    # mutate(Gene = sub(";.*$", "", sub("^.*upstream of gene ", "", Location))) |>
    group_by(Gene) |> 
    mutate(TSSCount = n()) |>
    # filter(n()>1) |> 
    # Both_Total is intensity of the peaks of all samples. 
    # Whichever sample has the highest peak is called "gTSS", all others are "alTSS".
    # gTSS = Genuine TSS, alTSS = Alternative TSS.
    arrange(desc(Both_Total)) |>
    mutate(SubClass = ifelse(row_number() == 1, "gTSS", "alTSS")) |>
    # The distance is from the TSS relative to the gTSS. For example, an alTSS
    # with distance -50, indicates that this alTSS is 50bp upstream of the gTSS.
    mutate(gTSSDistance = first(Position) - Position) |>
    # extract_sequence_from_position is a function that relies in the 
    #   global variable "genomes" and package "seqinr"
    # It will read the genome file found in genomes/$filename.gff
    mutate(Sequence = mapply(extract_sequence_from_position, Position, Position, Chromossome, Strand)) |>
    ungroup()
  
  
  write_xlsx(filtered, file.path(start_wd, filtered_folder, paste0(filename, "-pTSS.xlsx")))
  write_xlsx({filtered %>% filter(TSSCount > 1)}, file.path(start_wd, filtered_folder, paste0(filename, "-alTSS.xlsx")))
  generate_gtf(filtered, file.path(start_wd, filtered_folder, paste0(filename, "-pTSS")))
  generate_gtf({filtered %>% filter(TSSCount > 1)}, file.path(start_wd, filtered_folder, paste0(filename, "-alTSS")))
  
  # Generates the .fasta files to run MEME Motif later.
  print("4.2. Generating meme input")
  generate_memeinput(filtered)
  
  results[[filename]] <- c(
    Filenames = filename,
    Species = dados[i,]$Species,
    Assembly = gcf_id,
    Genes_total = nrow(gbff %>% filter(type == "gene")),
    Genes_ptss = length(unique(filtered$Gene)),
    Genes_multiple_ptss = length(unique({filtered %>% filter(TSSCount > 1)}$Gene)),
    Transcripts_5utr = nrow(ggdf %>% filter(utr == "5UTR")),
    Transcripts_leaderless = nrow(ggdf %>% filter(utr == "Leaderless"))
  )
  print(paste0("Finished ", filename))
  print("--------------")
  
   
  
  },
  error=function(e) {
    message('An Error Occurred')
    print(e)
  }
  )
}

r2 <- as.data.frame(t(as.data.frame(results)))
write_tsv(r2, file.path(start_wd, "results.tsv"))

print("Done.")

# nohup Rscript altss.R >> log-altss.txt 2>&1 &


