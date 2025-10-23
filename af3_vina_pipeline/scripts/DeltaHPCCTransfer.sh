#!/bin/bash -l
species="Arabidopsis DouglasFir Eucalyptus Human" #list of species

for s in $species
do
    for d in $(ls -f /mnt/gs21/scratch/borendun/IMPACTS2022/$s/equilibration/)
    do
        id=${d##../Simulations/$s/equilibration/}
        #uid=${id%%.fasta}
            echo $id
            #mkdir ../Simulations/$s/alphafold/$id
            echo rsync -auvrP /mnt/gs21/scratch/borendun/IMPACTS2022/$s/equilibration/$id/* boren1@login.delta.ncsa.illinois.edu:/scratch/bbft/boren1/IMPACTS/$s/equilibration/$id/
        fi
    done
done