#!/bin/bash

species_list=("DouglasFir" "Human" "Arabidopsis" "Eucalyptus")

for species in "${species_list[@]}"
do
    
    nohup python -u 04_run_docking.py "$species" > 04log/04${species}.out 2>&1 &
done

