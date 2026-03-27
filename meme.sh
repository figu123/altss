#!/bin/bash

memeinput_folder="./memeinput"
memeoutput_folder="./memeoutput"

mkdir "$memeoutput_folder"

if [ -z "$1" ]; then
  echo "Please specify the number of threads. Usage: bash meme.sh <number_of_threads>"
  exit 1
fi

# Sets directory location to where the script is found in.
cd "$(dirname "$0")"

# meme sequence.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 5 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0

find . -mindepth 2 -maxdepth 2 -iwholename "$memeinput_folder/*.fasta" | while read file; do
  echo "For file $file.."
  newfoldername="${file/#*\/}"
  newfoldername="${newfoldername/%\.*}"
  meme "$file" -dna -oc "$memeoutput_folder/$newfoldername" -time 18000 \
    -mod zoops -nmotifs 5 -minw 6 -maxw 50 -p $1 \
    -objfun classic -revcomp -markov_order 0 >> "$memeoutput_folder/$newfoldername.log.txt"
done;
