

tail -n +2 dados.tsv | while IFS=$'\t' read -r assembly species filename shortened project
do
    echo "$filename"
done


# samtools faidx 