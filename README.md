# altss

Pipeline for mapping and analysing alternative promoters in archea based on dRNA-seq data.

This is a repository organizing all the tools that I've been using throughout my scientific initiation.
Initial sequencing processing was done thanks to the pipeline of https://github.com/alanlorenzetti/frtc/ .
This Pipeline expects that you will have the same results as obtained by running FRTC.
It is also integrated with the TSSAR commandline tool, http://rna.tbi.univie.ac.at/TSSAR , having some modifications that were done in its clusterization as described in (IBRAHIM Amr, 2023) https://doi.org/10.11606/T.95.2023.tde-17082023-215744

-----

### Setup

1. Clone the repository. 

2. Inside the `data` folder, provide the following files:

2. 1. Create a tab separated value (.tsv) file called "dados.tsv". Use the exact columns as described below, substituting each row name by your sample. 

| Assembly | Species | Filenames | Shortened | 
| -------- | ------- | --------- | --------- | 
| GCF_000006805.1 | Halobacterium salinarum NRC-1 | hsalinarum | H. salinarum | 
| GCF_000025685.1 | Haloferax volcanii DS2 | hvolcanii | H. volcanii | 
| GCF_000009965.1 | Thermococcus kodakarensis KOD1 | tkodakarensis | T. kodakarensis |
| GCF_000012285.1 | Sulfolobus acidocaldarius DSM 639 | sacidocaldarius | S. acidocaldarius |

2. 2. Create a folder called `samples/`. Inside this folder, place the exact same term used in the column "`Filenames`" above, naming it `{filename}`.csv. This comma separated value should have information about each run, specifying the TEX+ and TEX- library, as well as a descriptor for identification. 

For example, given the species `Halobacterium salinarum NRC-1` of the table above, supposing it has 10 different conditions. Each condition has two libraries, one TEX+ and one TEX-. You can give any descriptor for the condition, but it is recommended keeping it simple, such as simplifying "17 hours, sample A" into "17ha". Here's how that can look like:

```csv
LibP,LibM,Output
SRR6953855,SRR6953869,ref1a
SRR6953856,SRR6953870,ref1b
SRR6953857,SRR6953871,ref2a
SRR6953858,SRR6953872,ref2b
SRR6953859,SRR6953865,17ha
SRR6953860,SRR6953866,17hb
SRR6953861,SRR6953867,37ha
SRR6953862,SRR6953868,37hb
SRR6953863,SRR6953873,86ha
SRR6953864,SRR6953874,86hb
```

2. 3. Change the pathes inside the `altss.R` script to your own.

Make one folder for each `{filename}` specified in dados.tsv (2.1.), each one should have the files from running the frtc pipeline ( https://github.com/alanlorenzetti/frtc/ ). 

For example, if you've used FRTC on the 4 archeal organisms described in the table `dados.tsv` above, it is assumed that the folder path of frtc will have the following files:

```
./frtc/
 |- ./hsalinarum
 |  |- ./coverage
 |  |- ./tssarinput
 |- ./hvolcanii
 |  |- ./coverage
 |  |- ./tssarinput
 |- ./tkodakarensis
 |  |- ./coverage
 |  |- ./tssarinput
 |- ./sacidocaldarius
 |  |- ./coverage
 |  |- ./tssarinput
```

-----

After all of this is done, this is how the final folder structure should look like

```
# input folders
(frtcfolder)/
 |- {filenames}/
 |  |- coverage
 |  |- tssarinput
./data/
 |- dados.tsv
 |- samples/
 |  |- {filenames}.csv
# scripts (WIP -- Prone to be changed)
./tssar/              # folder containing scripts used in tssar.sh
./1downloading-sras.sh   # <-- IF YOU NEED TO RUN FRTC, THIS SCRIPT HELPS in fetching raw sequencing files (sra)
./1fetch-genomes.sh   # run this first
./2tssar.sh           # then run this
./3altss.R            # then run this
./4orthofinder.sh     # then run this
./5meme.sh            # then run this
./meme-parsingresults # puts meme results into a table
./ortho-parsingresults# puts orthofinder results into a table
./quarto-candidates   # filters the table for visualizing alTSS candidates in Quarto
./6graphs.R           # then run this [[ WIP ]]
./7tree.sh (mafft.sh?)# then run this [[ WIP ]]
# output folders (automatically generated)
./results/
 |  |- altss - contains table information for each TSS
 |  |- memeinput - input data for meme for motif identification
 |  |- memeoutput - output data for meme
 |  |- merged - TSSAR results placed into one single spreadsheet
 |  |- misc - some graphs
 |  |- readcount - read count per gene (RSubread)
```

## Overview / Usage

Once setup is done, just run the scripts in order. Dependencies are the same as the ones described in FRTC, TSSAR, including the packages that can be found at the top of each .R script.

Things likely are broken, since this work is still in progress.

Other things you need to be able to run:
- basically all dependencies of FRTC and TSSAR (for running them in the first place)
- SRA toolkit (for prefetch + fasterq-dump in 1downloading-sras.sh)
- wget  (for 1fetch-genomes.sh)
- bedtools/samtools  (for 1fetch-genomes.sh to generate the .fai files for IGV)
- MEME (i used version 5.0.5)
- Orthofinder (i used version 2.5.4)
- IQTree
- Quarto (for visualizing coverage data in browser)
- python + package RangeHTTPServer (for igvserver.py)
- some R packages (tidyverse, rtracklayer, seqinr, Rsubread, ggplot2)

