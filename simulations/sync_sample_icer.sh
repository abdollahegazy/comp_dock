#!/bin/bash

DEST_HOST="hegazyab@rsync.hpcc.msu.edu"
DEST_BASE="/mnt/scratch/hegazyab/dock_comp/simulations/complexes"
SAMPLE_FILE="k_lowest_rmsd_sample"

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

    complex_type="${parts[1]}"  # MD_complexes_c3 or AF_complexes
    species="${parts[2]}"
    protein="${parts[3]}"
    ligand="${parts[4]}"

    # Determine file name based on complex type
    if [[ "$complex_type" == "MD_complexes_c3" ]]; then
        dest_filename="MD_complex.pdb"
    elif [[ "$complex_type" == "AF_complexes" ]]; then
        dest_filename="AF3_complex.pdb"
    else
        echo "Skipping erroneous line type in path: $line"
        continue
    fi

    remote_dir="${DEST_BASE}/${species}/${protein}/${ligand}"
    src_pdb="${line}/complex.pdb"
    dest_pdb="${remote_dir}/${dest_filename}"

    # Create remote directory
    ssh "$DEST_HOST" "mkdir -p '$remote_dir'" < /dev/null

    echo "Syncing: $src_pdb -> $dest_pdb"
    rsync -auvP "$src_pdb" "${DEST_HOST}:${dest_pdb}" < /dev/null

done < "$SAMPLE_FILE"
