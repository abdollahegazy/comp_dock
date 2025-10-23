#!/bin/bash

species_list=("DouglasFir" "Human" "Arabidopsis" "Eucalyptus")

for species in "${species_list[@]}"
do
    
    nohup python -u 03_make_prepare_ligand.py "$species" > 03log/03${species}.out 2>&1 &
done

