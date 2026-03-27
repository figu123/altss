#!/bin/bash
# Requer todos os scripts de classificação, clusterização de TSS
# Requer um arquivo csv contendo informações de "TEX+, TEX-, output" nomeada 'samples.csv'

# Usage: bash TSS_tudo.sh

# Sets directory location to where the script is found in.
cd "$(dirname "$0")"

# Procura por NC_ para extrair.
ids=($(ls *.bam 2>/dev/null | grep -o 'NC_[0-9]\+\.[0-9]\+' | sort -u))
ids+=($(ls *.bam 2>/dev/null | grep -o 'NZ_CP[0-9]\+\.[0-9]\+' | sort -u))


echo "IDs encontrados: ${ids[@]}"
# alternativamente, da pra passar os ids assim:  'ids=(NC_007488.2 NC_007489.1 NC_007490.2 NC_007493.2 NC_007494.2 NC_009007.1 NC_009008.1)'
# eval "$1"

for id in "${ids[@]}"; do
  echo "Baixando genoma $id"
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$id&rettype=fasta&retmode=text" -o "$id.fasta"
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$id&rettype=gbwithparts&retmode=text" -o "$id.gbff"
done
# Gera tabelas de proteína/RNA
for f in *.gbff; do mv "$f" "${f%.gbff}.gb"; done
python3 gbkToPttAndRnt.py ./*.gb

mkdir files;
while IFS=, read -r libP libM output; do
  [[ "$libP" == "LibP" ]] && continue
  echo "Rodando TSSAR para $libP e $libM"
  for id in "${ids[@]}"; do
    idn="${id%%.*}";
    libpbam="./${libP}-tssar-input-${id}.bam"
    libmbam="./${libM}-tssar-input-${id}.bam"
    gsize=$(grep -v "^>" "./${id}.fasta" | tr -d '\n' | wc -c)
    [[ -f "$libpbam" && -f "$libmbam" ]] || echo "[!] Arquivos de input estão faltando."
    
    samtools view -h -o "./${libP}-tssar-input-${id}.sam" "$libpbam"
    samtools view -h -o "./${libM}-tssar-input-${id}.sam" "$libmbam"

    perl ./TSSAR_12345GaV2_gff \
      --libP "./${libP}-tssar-input-$id.sam" \
      --libM "./${libM}-tssar-input-$id.sam" \
      --fasta "./$id.fasta" \
      --g_size $gsize \
      --verbose --noclean --nocluster --tmpdir '.'
    
    mv TSS_unclustered.bed "${idn}_${output}_TSS_u".bed

    perl ./gff_cluster.pl --score 'd' --i "${idn}_${output}_TSS_u".bed \
      --o "${idn}_${output}_TSS_c".bed --workdir '.' --tmpdir '.' --range 2

    cp ./smear.table ./files

    perl ./TSS_classification.pl --ptt ./"$idn".ptt --rnt ./"$idn".rnt \
      -i ./"${idn}_${output}_TSS_c".bed -o ./"${idn}_${output}_TSS".table --workdir '.'

  done
done < samples.csv


mkdir ../tssarresult
cp ./*.table ../tssarresult



