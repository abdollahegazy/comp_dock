#!/bin/bash -l
species="Arabidopsis DouglasFir Eucalyptus Human" #list of species

for s in $species
do
    for d in $(ls -f ../Simulations/$s/fasta/*.fasta)
    do
        id=${d##../Simulations/$s/fasta/}
        uid=${id%%.fasta}
        if [ ! -f ../Simulations/$s/alphafold/$uid/ranked_0.pdb ]; then
            echo $uid
            mkdir ../Simulations/$s/alphafold/$uid
            rsync -auvrP borendun@rsync.hpcc.msu.edu:/mnt/gs18/scratch/users/borendun/IMPACTS2022/$s/alphafold/$uid/* ../Simulations/$s/alphafold/$uid/ 
        fi
    done
done