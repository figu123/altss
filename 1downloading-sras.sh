#!/bin/bash

f=$(basename "$(pwd)")
echo "filename: $f"

sf="../samples/$f.csv" # make sure to specify the samples folder path here


if [ ! -f "$sf" ]; then
    echo "err: no samples.csv file for this filename."
    exit 1;
fi

mkdir sras 2>/dev/null
cp "$sf" ./sras/samples.csv
cd ./sras

awk -F, 'NR > 1 && $1 != "" && $2 != "" {print $1 "," $2}' samples.csv | tr ',' '\n' > ./acclist.txt

echo "acclist.txt is"
cat acclist.txt

prefetch --option-file acclist.txt &&
for file in ./**/*.sra; do mv "$file" .; done &&
for file in ./*.sra; do fasterq-dump --split-files "$file"; done &&
for file in ./*.fastq; do gzip "$file"; done; >> log-download.txt

