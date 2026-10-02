#!/bin/bash

# run iqtreeinput.R first
# conda activate mafft

cd "./results/iqtreeinput"

# https://mafft.cbrc.jp/alignment/software/

rm ./*.aln.fa ./*.mafft.log
find . -empty -delete

for f in *.fa; do
  out="${f%.fa}.aln.fa"

  # avoids rerunning if its already been run
  [[ -f "$out" ]] && continue
  [[ "$f" == *.aln.fa ]] && continue
  
  mafft "$f" > "$out" 2> "${f%.fa}.mafft.log" || {
      echo "MAFFT failed for $f" >&2
      rm "$out"
      continue
  }
done

# check
for f in *.fa; do
  out="${f%.fa}.aln.fa"
  [[ -f "$out" ]] && continue
  [[ "$f" == *.aln.fa ]] && continue
  echo "Did not run MAFFT for $f"
done

# conda activate iqtree
cd "/home/jfigueiredo/files/analises/altss/results"
mkdir iqtreeinput2
mkdir iqtreeoutput
cd ./iqtreeoutput
cp ../iqtreeinput/*.aln.fa ../iqtreeinput2
iqtree -s ../iqtreeinput2 -T 2
mv ../iqtreeinput2.* .

