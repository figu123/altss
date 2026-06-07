
# https://xml2.r-lib.org/
library("xml2")
library("tidyverse")
library("readxl")
library("writexl")

dados <- read_tsv("./dados.tsv")

# insert here the same values that you used for your memeinput.
upstream_count <- 100
downstream_count <- 10

for (i in 1:nrow(dados)) {
  fn <- dados[i,]$Filenames
  print(paste0("Parsing meme results for ", fn))
  memefile <- file.path(paste0("./memeoutput/", fn, "-pTSS-all/meme.xml"))
  
  if (!file.exists(memefile)) {
    print("Could not find meme results.")
    next
  }
  mfxml <- read_xml(memefile)
  
  # 1. reading motif data
  mdatainfo <- tibble(
    id = as.character(),
    seq = as.character(),
    width = as.numeric(),
    sites = as.numeric(),
    pval = as.numeric()
  )
  
  mdataxml <- xml_find_all(mfxml, ".//motif")
  for (j in 1:length(mdataxml)) {
    mdataattrs <- xml_attrs(mxml[j])[[1]]
    mdatainfo <- mdatainfo %>% add_row(
      id = mdataattrs[["id"]],
      seq = mdataattrs[["name"]],
      width = as.numeric(mdataattrs[["width"]]),
      sites = as.numeric(mdataattrs[["sites"]]),
      pval = as.numeric(mdataattrs[["p_value"]])
    ) 
  }
  mdatainfo
  
  # 2. appending to TSS what motifs they have
  memeresults <- tibble(
    chr = as.character(),
    pos = as.numeric(),
    strand = as.character(),
    seqid = as.character(),
    seqlen = as.numeric(),
    combinedmatchpvalue = as.numeric(),
    numsites = as.numeric(),
    sites = as.character()
  )
  
  # 3. information for each motif (will be used to make a .gtf)
  motifstable <- tibble(
    chr = as.character(),
    start = as.numeric(),
    end = as.numeric(),
    strand = as.character(),
    id = as.character(),
    seq = as.character(),
    pvalue = as.numeric(),
    parentid = as.character(),
    parentstrand = as.character(),
    parentposition = as.numeric()
  )
  
  seqs <- xml_find_all(mfxml, ".//sequence")
  for (j in 1:length(seqs)) {
    seqattrs <- xml_attrs(seqs[j])[[1]]
    
    pchr <- str_split_i(seqattrs[["name"]], "\\|", 1)
    pstr <- str_split_i(seqattrs[["name"]], "\\|", 2)
    ppos <- as.numeric( str_split_i(seqattrs[["name"]], "\\|", 3))
    
    memeresults <- memeresults %>% add_row(
      chr = pchr,
      pos = ppos,
      strand = pstr,
      seqid = seqattrs[["id"]],
      seqlen = as.numeric(seqattrs[["length"]])
    )
  }
  
  scannedsites <- xml_find_all(mfxml, ".//scanned_sites")
  for (j in 1:length(scannedsites)) {
    seqattrs <- xml_attrs(scannedsites[j])[[1]]
    
    motifs <- xml_children(scannedsites[j])
    motifsstr <- ""
    if (length(motifs) > 0) {
      for (k in 1:length(motifs)) {
        motifattrs <- xml_attrs(motifs[k])[[1]]
        motifsstr <- paste0(motifsstr,
                            motifattrs[["motif_id"]], "#",
                            motifattrs[["strand"]], "#",
                            motifattrs[["position"]], "#",
                            motifattrs[["pvalue"]], "%"
        )
        
        motifinfo <- mdatainfo %>% filter(id == motifattrs[["motif_id"]])
        pseqid <- seqattrs[["sequence_id"]]
        parentinfo <- xml_attrs(seqs[xml_attr(seqs, "id") == pseqid])[[1]][["name"]]
        pchr <- str_split_i(parentinfo, "\\|", 1)
        pstr <- str_split_i(parentinfo, "\\|", 2)
        ppos <- as.numeric(str_split_i(parentinfo, "\\|", 3))
        
        
        motifstable <- motifstable %>% add_row(
          chr = pchr,
          start = ifelse(motifattrs[["strand"]] == "plus", 
                         ppos - upstream_count + as.numeric(motifattrs[["position"]]),
                         ppos - upstream_count + as.numeric(motifattrs[["position"]]) - motifinfo$width
                        ),
          end = ifelse(motifattrs[["strand"]] == "plus", 
                         ppos - upstream_count + as.numeric(motifattrs[["position"]]) + motifinfo$width,
                         ppos - upstream_count + as.numeric(motifattrs[["position"]])
                        ),
          strand = ifelse(motifattrs[["strand"]] == "plus", "+", "-"),
          id = motifinfo$id,
          seq = motifinfo$seq,
          pvalue = as.numeric(motifattrs[["pvalue"]]),
          parentid = pseqid,
          parentstrand = pstr,
          parentposition = ppos
        )
        
      }
    }
    
    seqmotifs <- tibble(
      seqid = seqattrs[["sequence_id"]],
      combinedmatchpvalue = as.numeric(seqattrs[["pvalue"]]),
      numsites = as.numeric(seqattrs[["num_sites"]]),
      sites = motifsstr
    )
  
    memeresults <- rows_update(memeresults, seqmotifs, by = "seqid", unmatched = "ignore")
  }
  
  memeresults
  motifstable
  
  # 3. gtf with all the motifs for this species 
  motifstableextrainfo <- paste0('motifid \"', as.character(motifstable$id), '\";',
                     'seq \"', as.character(motifstable$seq), '\";',
                     'pval \"', as.character(motifstable$pvalue), '\";',
                     'parentstrand \"', as.character(motifstable$parentstrand), '\";',
                     'parentpos \"', as.character(motifstable$parentposition), '\";')
  memegtf <- data.frame(
    V1 = motifstable$chr,
    V2 = "MOTIF",
    V3 = "MOTIF",
    V4 = motifstable$start,
    V5 = motifstable$end,
    V6 = ".",
    V7 = motifstable$strand,
    V8 = ".",
    V9 = motifstableextrainfo)
  
  outfile_name <- paste0("./altss/", fn, "-meme-pTSS.gtf")
  write_tsv(memegtf, outfile_name, col_names = FALSE, escape = "none")
  
  # 4. applying to tss table the motif stuff.
  ptss <- read_xlsx(paste0("./altss/", fn, "-pTSS.xlsx"))
  
  ptss$memecombinedpval <- NA
  ptss$memenumsites <- NA
  ptss$memesites <- NA
  
  for (j in 1:nrow(ptss)) {
    rowresult <- memeresults %>% filter(pos == ptss[j,]$Position &
                                        chr == ptss[j,]$Chromossome &
                                        strand == ptss[j,]$Strand)
    
    if (nrow(rowresult) > 0) {
      ptss[j,]$memecombinedpval <- rowresult$combinedmatchpvalue
      ptss[j,]$memenumsites <- rowresult$numsites
      ptss[j,]$memesites <- rowresult$sites
    }
  }
  
  write_xlsx(ptss, paste0("./altss/", fn, "-pTSS.xlsx"))
  print("Saved on altss table.")
}



