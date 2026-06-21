#!/bin/bash

dados="./data/dados.tsv"

# expected folder structure
mkdir data/genomes >/dev/null 2>&1

# https://unix.stackexchange.com/questions/435814/compact-way-to-get-tab-separated-fields-into-variables
# tail skips the header. the rest is for reading a tab separated .tsv file.
tail -n +2 "$dados" | while IFS=$'\t' read -r col1 col2 col3; do 
	echo "Running for $col2 ($col3).."; 
	ncbi_accession="$col1"
	fn="$col2"
	
    if [ -f "./data/genomes/$fn.fna" ] && \
       # [ -f "./data/genomes/$fn.gbff" ] && \
       [ -f "./data/genomes/$fn.gaf" ] && \
       [ -f "./data/genomes/$fn.gff" ] && \
       [ -f "./data/genomes/$fn.faa" ] && \
       [ -f "./data/genomes/$fn-features.tsv" ]; then
        echo "Genome files for $ncbi_accession already exist"
        return 0
    fi
	
	# assembling the URL for NCBI path based on the GCF/GCA Assembly ID
    prefix=$(echo "$ncbi_accession" | cut -d_ -f1) # GCF or GCA
    digits=$(echo "$ncbi_accession" | cut -d_ -f2 | cut -d. -f1) # 9 digits before .

    p1=${digits:0:3}
    p2=${digits:3:3}
    p3=${digits:6:3}
	
    base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/$prefix/$p1/$p2/$p3"
    full_name=$(curl -s "$base_url/" | grep -o "${ncbi_accession}_[^\"/]*" | head -n 1)
    if [ -z "$full_name" ]; then
        echo "ERROR: Folder not found for $ncbi_accession in NCBI"
        exit 1
    fi

    ftp_path="$base_url/$full_name"

    fetch() {
        local file_suffix="$1"
        local output_name="$2"
        local remote_file="${full_name}${file_suffix}"
        local url="${ftp_path}/${remote_file}"

        echo "Downloading $remote_file -> $output_name ..."
        wget -qO- "$url" | gunzip -c > "./data/genomes/${output_name}"
    }

    fetch "_genomic.fna.gz" "$fn.fna"
    fetch "_gene_ontology.gaf.gz" "$fn.gaf"
    fetch "_genomic.gff.gz" "$fn-original.gff1"
    # fetch "_genomic.gbff.gz" "$fn.gbff"
    fetch "_feature_table.txt.gz" "$fn-features.tsv"
    fetch "_protein.faa.gz" "$fn.faa"  

    echo "Cleaning annotation .gff"
    awk -v OFS="\t" -v FS="\t" \
        '{if(/^#/){print}else{if($3=="CDS" || $3=="tRNA" || $3=="rRNA" || $3=="ncRNA"){print}}}' \
        "./data/genomes/$fn-original.gff1" > "./data/genomes/$fn.gff"


    echo "FTP download finished for $fn ($col3)"
    
	# extracting FASTAs from .fna
	
	# adapted from https://www.biostars.org/p/14067/
	awk -v fpath="./data/genomes/$fn" '
	$0 ~ /^>/ {                    # regexp match current line ("$0") for starting with with ">" symbol (FASTA header)
	  if (out) close(out)          # close previous file
	  name=$1                      # first word of header
	  sub(/^>/,"",name)            # removes ">"
	  out=fpath "-" name ".fasta"  # prefix with species name (fn)
	}
	{print >> out}                 # if regexp doesnt match, write sequence lines to file
	' "./data/genomes/$fn.fna"     # input file, complete genome
	
	
	echo "-----------------------------------"
	
done


# samtools faidx 

exit 0;
