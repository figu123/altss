#!/bin/bash

cd "./iqtreeinput"

# https://mafft.cbrc.jp/alignment/software/

for f in *.fa; do
  mafft --auto "$f" > "${f%.fa}.aln.fa"
done

#for f in *.aln.fa; do
#  iqtree3 -s "$f" -m MFP
#done;

# no boostrap cus most orthogroups are just 3 filenames

#for f in *.aln.fa; do
#  ./iqtree/bin/iqtree3.exe -s "$f" -m MFP
#done;

# cat *.treefile > all_trees.tre
# java -jar astral.5.7.8.jar -i all_trees.tre -o species_tree.tre