
# https://xml2.r-lib.org/
library("xml2")
library("tidyverse")

dados <- read_tsv("./dados.tsv")

# insert here the same values that you used for your memeinput.
upstream_count <- 100
downstream_count <- 10

for (i in 1:nrow(dados)) {
  fn <- dados[i,]$Filenames
  print(paste0("Parsing meme results for ", fn))
  xmlfile <- file.path(paste0("./results/memeoutput/", fn, "-pTSS-all/meme.xml"))
  
  if (!file.exists(xmlfile)) {
    print("Could not find meme results.")
    next
  }
  myxml <- read_xml(xmlfile)
  
  datasequences <- tibble(
    id = as.character(),
    chr = as.character(),
    strand = as.character(),
    position = as.numeric()
  )
  datamotifs <- tibble(
    id = as.character(),
    seq = as.character(),
    width = as.numeric(),
    sites = as.numeric(),
    pval = as.numeric()
  )
  motifsites <- tibble(
    motif_id = as.character(),
    motif_seq = as.character(),
    
    sequence = as.character(),
    flankedsequence = as.character(),
    pvalue = as.numeric(),
    position = as.numeric(),
    parent_seqid = as.character(),
    
    chr = as.character(),
    start = as.numeric(),
    end = as.numeric(),
    strand = as.character(),
    parent_position = as.numeric()
  )
  
  
  myxml_sequences <- xml_find_all(myxml, ".//sequence")
  for (j in 1:length(myxml_sequences)) {
    sequences_attrs <- xml_attrs(myxml_sequences[j])[[1]]
    pchr <- str_split_i(sequences_attrs[["name"]], "\\|", 1)
    pstr <- str_split_i(sequences_attrs[["name"]], "\\|", 2)
    ppos <- as.numeric( str_split_i(sequences_attrs[["name"]], "\\|", 3))
    
    datasequences <- datasequences %>% add_row(
      id = as.character(sequences_attrs[["id"]]),
      chr = pchr,
      strand = pstr,
      position = ppos
    )
  }
  datasequences
  
  myxml_motifs <- xml_find_all(myxml, ".//motif")
  for (j in 1:length(myxml_motifs)) {
    motifs_attrs <- xml_attrs(myxml_motifs[j])[[1]]
    datamotifs <- datamotifs %>% add_row(
      id = motifs_attrs[["id"]],
      seq = motifs_attrs[["name"]],
      width = as.numeric(motifs_attrs[["width"]]),
      sites = as.numeric(motifs_attrs[["sites"]]),
      pval = as.numeric(motifs_attrs[["p_value"]])
    )
    
    # here i learned that "contributing_site" is what describes the motifs.
    sites <- xml_children(xml_find_all(myxml_motifs[j], ".//contributing_sites"))
    for (k in 1:length(sites)) {
      sitesequence <- paste0(xml_attr(xml_find_all(sites[k], ".//letter_ref"), "letter_id"), collapse = "")
      fullsiteseq <- paste0(xml_text(xml_find_first(sites[k], "left_flank")),
                            sitesequence,
                            xml_text(xml_find_first(sites[k], "right_flank")))
      
      pseqid <- xml_attr(sites[k], "sequence_id")
      pdata <- datasequences[datasequences$id == pseqid, ]
      pchr <- pdata$chr
      pstr <- pdata$strand
      ppos <- pdata$position
      
      motifsites <- motifsites %>% add_row(
        motif_id = as.character(motifs_attrs[["id"]]),
        motif_seq = as.character(motifs_attrs[["name"]]),
        sequence = as.character(sitesequence),
        flankedsequence = as.character(fullsiteseq),
        pvalue = as.numeric(xml_attr(sites[k], "pvalue")),
        position = as.numeric(xml_attr(sites[k], "position")),
        parent_seqid = as.character(pseqid),
        chr = pchr,
        start = ifelse(pstr == "+", 
                       ppos - upstream_count + position,
                       ppos + upstream_count - position - nchar(sequence)),
        end = start + nchar(sequence),
        strand = pstr,
        parent_position = ppos
      )
    }
  }
  # datamotifs
  # motifsites
  # write_tsv(motifsites, paste0("./debug/meme-", fn, ".tsv"))

  motifsitesextrainfo <- paste0('motifid \"', as.character(motifsites$motif_id), '\";',
                     'motifseq \"', as.character(motifsites$motif_seq), '\";',
                     'pval \"', as.character(motifsites$pvalue), '\";',
                     'sequence \"', as.character(motifsites$sequence), '\";',
                     'fullseq \"', as.character(motifsites$flankedsequence), '\";',
                     'parentposition \"', as.character(motifsites$parent_position), '\";')
  
  memegtf <- data.frame(
    V1 = motifsites$chr,
    V2 = "MOTIF",
    V3 = "MOTIF",
    V4 = motifsites$start,
    V5 = motifsites$end,
    V6 = ".",
    V7 = motifsites$strand,
    V8 = ".",
    V9 = motifsitesextrainfo)

  outfile_name <- paste0("./results/altss/", fn, "-meme-pTSS.gtf")
  write_tsv(memegtf, outfile_name, col_names = FALSE, escape = "none")
  
  ptss <- read_tsv(paste0("./results/altss/", fn, "-pTSS.tsv"))
  ptss$memenumsites <- 0
  
  for (j in 1:nrow(motifsites)) {
    motifexists <- ptss$Position == motifsites$parent_position[j] &
                   ptss$Chromossome == motifsites$chr[j] &
                   ptss$Strand == motifsites$strand[j]
    
    ptss$memenumsites[motifexists] <- ptss$memenumsites[motifexists] + 1
  }

  write_tsv(ptss, paste0("./altss/", fn, "-pTSS.tsv"))
  print("Generated meme.gtf and added numsites column for altss table.")

  # 2. appending to TSS what motifs they have
  # memeresults <- tibble(
  #   chr = as.character(),
  #   pos = as.numeric(),
  #   strand = as.character(),
  #   seqid = as.character(),
  #   seqlen = as.numeric(),
  #   combinedmatchpvalue = as.numeric(),
  #   numsites = as.numeric(),
  #   sites = as.character()
  # )
  
  # 3. information for each motif (will be used to make a .gtf)
  # motifstable <- tibble(
  #   chr = as.character(),
  #   start = as.numeric(),
  #   end = as.numeric(),
  #   strand = as.character(),
  #   id = as.character(),
  #   seq = as.character(),
  #   pvalue = as.numeric(),
  #   parentid = as.character(),
  #   parentstrand = as.character(),
  #   parentposition = as.numeric()
  # )
  # 
  # seqs <- xml_find_all(mfxml, ".//sequence")
  # for (j in 1:length(seqs)) {
  #   seqattrs <- xml_attrs(seqs[j])[[1]]
  #   
  #   pchr <- str_split_i(seqattrs[["name"]], "\\|", 1)
  #   pstr <- str_split_i(seqattrs[["name"]], "\\|", 2)
  #   ppos <- as.numeric( str_split_i(seqattrs[["name"]], "\\|", 3))
  #   
  #   memeresults <- memeresults %>% add_row(
  #     chr = pchr,
  #     pos = ppos,
  #     strand = pstr,
  #     seqid = seqattrs[["id"]],
  #     seqlen = as.numeric(seqattrs[["length"]])
  #   )
  # }
  # 
  # scannedsites <- xml_find_all(mfxml, ".//scanned_sites")
  # for (j in 1:length(scannedsites)) {
  #   seqattrs <- xml_attrs(scannedsites[j])[[1]]
  #   
  #   motifs <- xml_children(scannedsites[j])
  #   motifsstr <- ""
  #   if (length(motifs) > 0) {
  #     for (k in 1:length(motifs)) {
  #       motifattrs <- xml_attrs(motifs[k])[[1]]
  #       motifsstr <- paste0(motifsstr,
  #                           motifattrs[["motif_id"]], "#",
  #                           motifattrs[["strand"]], "#",
  #                           motifattrs[["position"]], "#",
  #                           motifattrs[["pvalue"]], "%"
  #       )
  #       
  #       motifinfo <- mdatainfo %>% filter(id == motifattrs[["motif_id"]])
  #       pseqid <- seqattrs[["sequence_id"]]
  #       parentinfo <- xml_attrs(seqs[xml_attr(seqs, "id") == pseqid])[[1]][["name"]]
  #       pchr <- str_split_i(parentinfo, "\\|", 1)
  #       pstr <- str_split_i(parentinfo, "\\|", 2)
  #       ppos <- as.numeric(str_split_i(parentinfo, "\\|", 3))
  #       
  #       
  #       motifstable <- motifstable %>% add_row(
  #         chr = pchr,
  #         start = ifelse(motifattrs[["strand"]] == "plus", 
  #                        ppos - upstream_count + as.numeric(motifattrs[["position"]]),
  #                        ppos - upstream_count + as.numeric(motifattrs[["position"]])
  #                       ),
  #         end = ifelse(motifattrs[["strand"]] == "plus", 
  #                        ppos - upstream_count + as.numeric(motifattrs[["position"]]) + motifinfo$width,
  #                        ppos - upstream_count + as.numeric(motifattrs[["position"]]) + motifinfo$width
  #                       ),
  #         strand = ifelse(motifattrs[["strand"]] == "plus", "+", "-"),
  #         id = motifinfo$id,
  #         seq = motifinfo$seq,
  #         pvalue = as.numeric(motifattrs[["pvalue"]]),
  #         parentid = pseqid,
  #         parentstrand = pstr,
  #         parentposition = ppos
  #       )
  #       
  #     }
  #   }
  #   
  #   seqmotifs <- tibble(
  #     seqid = seqattrs[["sequence_id"]],
  #     combinedmatchpvalue = as.numeric(seqattrs[["pvalue"]]),
  #     numsites = as.numeric(seqattrs[["num_sites"]]),
  #     sites = motifsstr
  #   )
  # 
  #   memeresults <- rows_update(memeresults, seqmotifs, by = "seqid", unmatched = "ignore")
  # }
  # 
  # memeresults
  # motifstable
  # 
  # # 3. gtf with all the motifs for this species 
  # motifstableextrainfo <- paste0('motifid \"', as.character(motifstable$id), '\";',
  #                    'seq \"', as.character(motifstable$seq), '\";',
  #                    'pval \"', as.character(motifstable$pvalue), '\";',
  #                    'parentstrand \"', as.character(motifstable$parentstrand), '\";',
  #                    'parentpos \"', as.character(motifstable$parentposition), '\";')
  # memegtf <- data.frame(
  #   V1 = motifstable$chr,
  #   V2 = "MOTIF",
  #   V3 = "MOTIF",
  #   V4 = motifstable$start,
  #   V5 = motifstable$end,
  #   V6 = ".",
  #   V7 = motifstable$strand,
  #   V8 = ".",
  #   V9 = motifstableextrainfo)
  # 
  # outfile_name <- paste0("./altss/", fn, "-meme-pTSS.gtf")
  # write_tsv(memegtf, outfile_name, col_names = FALSE, escape = "none")
  # 
  # # 4. applying to tss table the motif stuff.
  # ptss <- read_xlsx(paste0("./altss/", fn, "-pTSS.xlsx"))
  # 
  # ptss$memecombinedpval <- NA
  # ptss$memenumsites <- NA
  # ptss$memesites <- NA
  # 
  # for (j in 1:nrow(ptss)) {
  #   rowresult <- memeresults %>% filter(pos == ptss[j,]$Position &
  #                                       chr == ptss[j,]$Chromossome &
  #                                       strand == ptss[j,]$Strand)
  #   
  #   if (nrow(rowresult) > 0) {
  #     ptss[j,]$memecombinedpval <- rowresult$combinedmatchpvalue
  #     ptss[j,]$memenumsites <- rowresult$numsites
  #     ptss[j,]$memesites <- rowresult$sites
  #   }
  # }
  # 
  # write_xlsx(ptss, paste0("./altss/", fn, "-pTSS.xlsx"))
  # print("Saved on altss table.")
}



