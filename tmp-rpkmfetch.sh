#!/bin/bash

# fetches the mapped reads from flagstat command.

frtcf="$1"
fn="$2"
logfile="./tmp-rpkmnfetch.log"

# echo "Running command: $frtcf $fn"
[ -z "$1" ] && echo "Please specify the FRTC folder path as the first argument" >> "$logfile" && exit 1;
[ -z "$2" ] && echo "Please specify the filename of your species as the second argument" >> "$logfile" && exit 1;

echo "$(date) Running $frtcf $fn" >> "$logfile"
echo -e "sample\tmapped" > "./results/readcount/$fn-mappedreads.tsv"
find "$frtcf/$fn/coverage" -iname "*.bam" | while read file; do
  # fs="${file%.*}.flagstat.txt"
  # [ ! -f "$fs" ] && echo "Running flagstat for $file" >> "$logfile" && samtools flagstat "$file" > "$fs"
  # mapped=`head -n 1 "$fs" | awk '{print $1}'`
  mapped=`samtools view -c "$file"`
  bn=`basename "$file"`
  echo -e "${bn}\t${mapped}" >> "./results/readcount/$fn-mappedreads.tsv"
done;


# does the same protocol as described in https://www.nature.com/articles/srep09209#Sec2

awk -v FS="\t" -v OFS="\t" '
NR > 1 {
    chr=$2; pos=$3; strand=$4; rep=$12
    
    if(strand=="+"){
        print chr, pos, pos+50, strand, rep
    } else {
        print chr, pos-50, pos, strand, rep
    }
}' "./results/altss/$fn-pTSS.tsv" > "./results/readcount/$fn-window50.bed"

# https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
# https://bedtools.readthedocs.io/en/latest/content/tools/map.html
tail -n +2 "./data/samples/$fn.csv" | while IFS=',' read texp texm repid; do 
  echo "Running for $repid";
  # output | c1: chromossome; c2: start; c3: end; c4: rep; c5: count
  # awk collapses all columns that have identical chr, start, end and repid.
                   
  bedtools intersect -a "./results/readcount/$fn-window50.bed" \
                   -b "$frtcf/$fn/fiveprime/$texp-5primeprofile-fwd.bedgraph" \
                   "$frtcf/$fn/fiveprime/$texp-5primeprofile-rev.bedgraph" \
                   "$frtcf/$fn/fiveprime/$texm-5primeprofile-fwd.bedgraph" \
                   "$frtcf/$fn/fiveprime/$texm-5primeprofile-rev.bedgraph" \
                   -wa -wb | awk -v FS="\t" -v OFS="\t" '{
                        key=$1"\t"$2"\t"$3"\t"$4"\t"$5;
                        sum[key]+=$9
                    }
                    END{
                        for (k in sum) print k, sum[k]
                    }' | awk -v repid="$repid" '$5 == repid' | sort -k1,1 -k2,2n > "./results/readcount/$fn-win50-$repid.bed"

done

rm "./results/readcount/$fn-window50.bed"
cat "./results/readcount/$fn-win50-"*.bed | sort -k1,1 -k2,2n > "./results/readcount/$fn-window50.bed"
rm "./results/readcount/$fn-win50-"*.bed

