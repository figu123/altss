# Script for parsing data post running FRTC for multiple files,
# which gives tables that can be used later for plotting boxplots, etc.

setwd("~/files")

library("tidyverse")

frtc_folder <- "./frtc"

directories <- list.dirs(path = frtc_folder, recursive = FALSE)

readscounts_results <- data.frame(
  species = as.character(),
  library = as.character(),
  paired = as.character(),             # TRUE/FALSE
  reads_raw = as.numeric(),            # trimmed.log
  reads_trimmed = as.numeric(),        # trimmed.log
  reads_aligned = as.numeric(),        # readCounts.txt     <-- generated in FRTC step 6.1 : Trimmed + Aligned + MMR
  reads_aligned_uniq = as.numeric()    # readCountsUniq.txt <-- generated in FRTC step 6.2 : Trimmed + Aligned Only Once
)
alignment_results_full <- data.frame(
  species = as.character(),
  library = as.character(),
  type = as.character(),               # unpaired_R1, unpaired_R2, paired
  overall_align_rate = as.numeric()    # sam/align.log
)


# This script parses data from three locations:
#  1. misc/readCounts.txt misc/readCountsUniq.txt --- Contains information of total aligned reads (Paired + Unpaired Libraries)
#  2. trimmed/*.log   --- Information about trimming
#  3. sam/*.log       --- Information about aligning

# 1 and 2 is used in generating readscounts_results dataframe
# 3 is used in generating alignment_results_full dataframe

# In the occasion one of these locations is missing, the species will be skipped.

sink(file = "./check-alignments.log")

print(paste0("Running at ", getwd()))
print(paste0("FRTC folder is ", frtc_folder))
print(paste0("Found ", str(length(directories)), " directories"))

for (d in directories) {
  # lets go through step by step. starting by trimming information.
  # it should have a single library.log file, for each library of the experiment
  # if FRTC was run in paired-end=y, it should look like:
  # # Input Read Pairs: 1570241 Both Surviving: 889164 (56.63%) Forward Only Surviving: 47459 (3.02%) Reverse Only Surviving: 18003 (1.15%) Dropped: 615615 (39.21%)
  # otherwise, it looks like:
  # # Input Reads: 54118564 Surviving: 41919963 (77,46%) Dropped: 12198601 (22,54%)

  species = sub(".*/", "", d)
  print(paste0("Running for species: ", species))
  
  if (!file.exists(paste0(d, "/misc/readCounts.txt"))) {
    print("Couldn't find readCounts file. Skipping.")
  } else if (!dir.exists(paste0(d, "/trimmed"))) {
    print("Couldn't find trimmed folder. Skipping.")
  } else {
    read_counts = read_tsv(paste0(d, "/misc/readCounts.txt"), col_names = FALSE, show_col_types = FALSE)
    colnames(read_counts) = c("library", "counts")
    read_counts_uniq = read_tsv(paste0(d, "/misc/readCountsUniq.txt"), col_names = FALSE, show_col_types = FALSE)
    colnames(read_counts_uniq) = c("library", "counts")
    
    logfiles = list.files(paste0(d, "/trimmed"), pattern = "\\.log")
    for (lf in logfiles) {
      logcontents = readr::read_file(paste0(d, "/trimmed/", lf))
      logcontents = str_split(logcontents, "\n")[[1]]
      
      paired_end = any(grepl("Input Read Pairs", logcontents))
      reads_raw = -1
      reads_trimmed = -1
      library = sub(".log", "", lf)
      
      if (paired_end) {
        logcontents = logcontents[grep("Input Read", logcontents)]
        # "Input Read Pairs: 1570241 Both Surviving: 889164 (56.63%) Forward Only Surviving: 47459 (3.02%) Reverse Only Surviving: 18003 (1.15%) Dropped: 615615 (39.21%)"
        reads_raw = as.numeric(sub(".*Pairs: ", "", sub(" Both.*", "", logcontents)))      # 1570241
        reads_dropped = as.numeric(sub(" \\(.*", "", sub(".*Dropped: ", "", logcontents))) # 615615
        reads_trimmed = reads_raw - reads_dropped # 954626 = Both Surviving + Fwd Surviving + Rev Surviving
      } else {
        logcontents = logcontents[grep("Input Read", logcontents)]
        # "Input Reads: 58212099 Surviving: 45513643 (78,19%) Dropped: 12698456 (21,81%)"
        reads_raw = as.numeric(sub(".*eads: ", "", sub(" Surv.*", "", logcontents)))      # 58212099
        reads_trimmed = as.numeric(sub(".*viving: ", "", sub(" \\(.*", "", logcontents))) # 45513643
      }paste0(d, "/misc/readCounts.txt")
      
      reads_aligned = read_counts[grep(library, read_counts$library),]$counts
      reads_aligned_uniq = read_counts_uniq[grep(library, read_counts_uniq$library),]$counts
      
      readscounts_results <- rbind(readscounts_results,
                                   data.frame(
                                     species = species,
                                     library = library,
                                     paired = paired_end,
                                     reads_raw = reads_raw,
                                     reads_trimmed = reads_trimmed,       
                                     reads_aligned = reads_aligned,
                                     reads_aligned_uniq = reads_aligned_uniq
                                   ))
    }
  }
  
  # final information we will be parsing is alignment information
  # found in the sam folder, paired-end results should have:
  # # library-paired.log      containing aligned paired reads
  # # library-unpaired_R1.log containing aligned reads R1 that had its R2 eliminated by the trimming (unpaired or orphan R1)
  # # library-unpaired_R2.log containgin aligned reads R2 that had its R1 eliminated by the trimming (unpaired or orphan R2)

  if (!dir.exists(paste0(d, "/sam"))) {
    print("Couldn't find sam folder. Skipping.")
  } else {
    logfiles2 <- list.files(paste0(d, "/sam"), pattern = "\\.log")
    for (lf in logfiles2) {
      logcontents = readr::read_file(paste0(d, "/sam/", lf))
      logcontents = str_split(logcontents, "\n")[[1]]
      
      type = sub(".*-", "", sub(".log", "", lf))
      library = sub("-.*", "", lf)
      
      if (type == "paired") {
        # [1] "3199120 reads; of these:"
        # [2] "  3199120 (100.00%) were paired; of these:"              
        # [3] "    378928 (11.84%) aligned concordantly 0 times"
        # [4] "    2615303 (81.75%) aligned concordantly exactly 1 time"
        # [5] "    204889 (6.40%) aligned concordantly >1 times"         
        # [6] "88.16% overall alignment rate"
        overall_alignment = logcontents[grep("overall alignment", logcontents)]
        overall_alignment = as.numeric(sub("\\% overall.*", "", overall_alignment)) # 88.16
      } else {
        # [1] "75203 reads; of these:"                     
        # [2] "  75203 (100.00%) were unpaired; of these:" 
        # [3] "    10336 (13.74%) aligned 0 times"        
        # [4] "    51064 (67.90%) aligned exactly 1 time"  
        # [5] "    13803 (18.35%) aligned >1 times"        
        # [6] "86.26% overall alignment rate"             
        overall_alignment = logcontents[grep("overall alignment", logcontents)]
        overall_alignment = as.numeric(sub("\\% overall.*", "", overall_alignment)) # 86.26
      }
      
      
      alignment_results_full <- rbind(alignment_results_full,
                                      data.frame(
                                        species = species,
                                        library = library,
                                        type = type,
                                        overall_align_rate = overall_alignment
                                      ))
      
      
    }
  }
  
  print("----------------------")
}


print("Saving results.")
write_tsv(alignment_results_full, "check-align-results-align.tsv")
write_tsv(readscounts_results, "check-align-results-reads.tsv")

print("Done.")
sink()


