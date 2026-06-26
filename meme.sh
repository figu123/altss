#!/bin/bash


# perosnal note: remember to activate your enviroment
# source ~/miniconda3/bin/activate; conda activate meme-env

memeinput_folder="./results/memeinput"
memeoutput_folder="./results/memeoutput"

mkdir "$memeoutput_folder" >/dev/null 2>&1

if [ -z "$1" ]; then
  echo "Please specify the number of threads. Usage: bash meme.sh <number_of_threads>"
  exit 1
fi

# Sets directory location to where the script is found in.
cd "$(dirname "$0")"

# meme sequence.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 5 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0

#find . -mindepth 2 -maxdepth 2 -iwholename "$memeinput_folder/*.fasta" | while read file; do
#  echo "For file $file.."
#  newfoldername="${file/#*\/}"
#  newfoldername="${newfoldername/%\.*}"
#  meme "$file" -dna -oc "$memeoutput_folder/$newfoldername" -time 18000 \
#    -mod zoops -nmotifs 5 -minw 5 -maxw 10 -p $1 \
#    -objfun classic -revcomp -markov_order 0
#done;

find . -iwholename "$memeinput_folder/*.fasta" | \
while read -r file; do
  newfoldername="${file/#*\/}"
  newfoldername="${newfoldername/%\.*}"
  echo "For file $file.. (Saving results at $memeoutput_folder/$newfoldername"
  # -revcomp should not be used, because we already give the reverse-complementary of reverse strands.
  # doing so would be... finding motifs in the reverse strand of the TSS, which is counterintuitive.
  meme "$file" -dna -oc "$memeoutput_folder/$newfoldername" -brief 100000 -time 18000 \
    -mod zoops -nmotifs 4 -minw 5 -maxw 10 -minsites 40 -p "$1" \
    -objfun classic -markov_order 0 < /dev/null
done
