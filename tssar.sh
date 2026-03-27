#!/bin/bash

# This script generates the "tssarresult" folder that is used by the main script
# Full explanation of what the folders below are can be found in the main script
frtc_folder="/home/jfigueiredo/files/frtc/"
samples_folder="/home/jfigueiredo/files/samples/"
tssar_folder="./tssar"

# Sets directory location to where the script is found in.
cd "$(dirname "$0")"

while IFS=$'\t' read -r assembly species filename; do
  # Skip the first line.
  [[ "${assembly,,}" == "assembly"* ]] && continue;
  fp="$frtc_folder""/$filename"
  fs="$samples_folder""/$filename.csv"
  [[ -d "$fp""/tssarresult" ]] && echo "TSSAR was already run for $filename. Skipping." && continue;
  [[ ! -f "$fs" ]] && echo "Could not find samplefolder/$filename.csv. Stopping." && continue;
  # Copying tssar program inside tssarinput, and running it.
  cp "$tssar_folder/"* "$fp/tssarinput"
  cp "$fs" "$fp/tssarinput/samples.csv"
  echo "Running TSSAR for $filename..."
  echo "call: bash $fp/tssarinput/TSS_tudo.sh"
  bash "$fp/tssarinput/TSS_tudo.sh" >> "$fp/tssarinput/tssar-log.txt"
done < "./dados.tsv"
