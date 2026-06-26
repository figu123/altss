#!/bin/bash

# Fetches coverage files from frtc folder, places them inside ./data

frtc_folder="~/files/frtc"
dados="./data/dados.tsv"

mkdir -p ./data/genomes/frtc >/dev/null 2>&1

# https://unix.stackexchange.com/questions/435814/compact-way-to-get-tab-separated-fields-into-variables
# tail skips the header. the rest is for reading a tab separated .tsv file.
tail -n +2 "$dados" | while IFS=$'\t' read -r col1 col2 col3; do 
  mkdir ./data/genomes/frtc/coverage >/dev/null 2>&1
  echo "$frtc_folder/$col2/coverage/*.bedgraph -> ./data/frtc/$col2/coverage/"
  cp "$frtc_folder/$col2"/coverage/*.bedgraph "./data/frtc/$col2/coverage/"
done

