# altss

Pipeline for identification of alternative promoters of multiple prokaryotic genomes based on dRNA-seq data.

This is a repository organizing all the tools that I've been using throughout my scientific initiation.
It has the same requirements described on https://github.com/alanlorenzetti/frtc/ , and relies on the output files of this program for running it.
It is also integrated with the TSSAR commandline tool, http://rna.tbi.univie.ac.at/TSSAR , having some modifications that were done in its clusterization as described in (IBRAHIM Amr, 2023) https://doi.org/10.11606/T.95.2023.tde-17082023-215744

TODO: Maybe reduce the amount of packages that are requirements in R.

### Setup

1. Clone the repository. 

2. Create a tab separated value (.tsv) file called "dados.tsv". Use the exact columns as described below, substituting each row name by your sample. 

| Assembly | Species | Filenames | Shortened | 
| -------- | ------- | --------- | --------- | 
| GCF_000006805.1 | Halobacterium salinarum NRC-1 | hsalinarum | H. salinarum | 
| GCF_000025685.1 | Haloferax volcanii DS2 | hvolcanii | H. volcanii | 
| GCF_000009965.1 | Thermococcus kodakarensis KOD1 | tkodakarensis | T. kodakarensis |
| GCF_000012285.1 | Sulfolobus acidocaldarius DSM 639 | sacidocaldarius | S. acidocaldarius |

3. Create a folder called `samples/`. Inside this folder, place the exact same term used in the column "`Filenames`" above, naming it {filename}.csv. This comma separated value should have information about each run, specifying the TEX+ and TEX- library, as well as a descriptor for identification. For example, given the species `Halobacterium salinarum NRC-1` of the table above, supposing it has 10 different conditions. Each condition has two libraries, one TEX+ and one TEX-. You can give any descriptor for the condition, such as simplifying "17 hours, sample A" into "17ha". Here's how that can look like:

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

4. Finally, specify the folder that contains results of the pipeline frtc ( https://github.com/alanlorenzetti/frtc/ ) for each filename. frtc is a pipeline that produces various outputs, only two are of interest: `tssarinput` and `coverage`. It is assumed that this folder will also follow the {filename} naming convention established above. So for example, if you've used FRTC on the 4 archeal organisms described in the table `dados.tsv` above, it is assumed that the folder path of frtc will have the following files:

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

You must specify this folder inside the script altss.R, in the frtc_folder variable.
Other variables in there can also be customized, at least until you reach the "hardcoded variables" spot.

After all of this is done, this is how the final folder structure should look like

```
# input folders
./dados.tsv
./samples/{filename}.csv
(frtcfolder)/{filename}/coverage and tssarinput
# scripts
./2tssar.sh
./3altss.R
./4meme.sh
./5graphs.R
./6ortho.sh
./7tree.sh
# output folders (automatically generated)
./altss/
./genomes/
./memeinput/
./merged/
./misc/
./readcount/
./tssar/
```

## Overview / Usage

Once setup is done, just run the scripts in order. Dependencies are the same as the ones described in FRTC, TSSAR, including the packages that can be found at the top of each .R script.

