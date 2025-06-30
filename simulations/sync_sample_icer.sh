#!/bin/bash

DEST_HOST="hegazyab@rsync.hpcc.msu.edu"
DEST_BASE="/mnt/scratch/hegazyab/dock_comp/simulations/complexes"
SAMPLE_FILE="sample"

if [[ ! -f "$SAMPLE_FILE" ]]; then
    echo "Error: sample file '$SAMPLE_FILE' not found."
    exit 1
fi

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" ]] && continue

    IFS='/' read -ra parts <<< "$line"
    if [[ ${#parts[@]} -lt 5 ]]; then
        echo "Skipping malformed path: $line"
        continue
    fi

    species="${parts[2]}"
    protein="${parts[3]}"
    ligand="${parts[4]}"

    remote_dir="${DEST_BASE}/${species}/${protein}/${ligand}"

    # Create remote directory (using SSH)
    ssh "$DEST_HOST" "mkdir -p '$remote_dir'" < /dev/null

    # Paths for rsync
    md_src="${line}/complex.pdb"
    md_dest="${remote_dir}/MD_complex.pdb"

    echo "Syncing MD: $md_src -> $md_dest"
    rsync -auvP "$md_src" "${DEST_HOST}:${md_dest}" < /dev/null

    af_line="${line/MD_complexes_c3/AF_complexes}"
    af_src="${af_line}/complex.pdb"
    af_dest="${remote_dir}/AF3_complex.pdb"

    echo "Syncing AF: $af_src -> $af_dest"
    rsync -auvP "$af_src" "${DEST_HOST}:${af_dest}" < /dev/null

done < "$SAMPLE_FILE"
