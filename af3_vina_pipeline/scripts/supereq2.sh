#!/bin/bash
#this program runs eq2.sh for all the folders in NSFIMPACTS/[SpeciesName]/Equilibration, most recently set to a runtime of 8 hours

species="Arabidopsis DouglasFir Eucalyptus Human" #list of species
for s in $species
do
	for d in $(ls /mnt/gs21/scratch/borendun/IMPACTS2022/$s/equilibration)
	do
		id=${d##/mnt/gs21/scratch/borendun/IMPACTS2022/$s/equilibration/}
                n=`squeue --me | grep -c eq`
       		echo $n
       		if [ $n -gt 399 ]
        	then
                	exit 0
        	fi
		if [ ! -f /mnt/gs21/scratch/borendun/IMPACTS2022/$s/equilibration/$id/eq2.coor ];then
			echo $id
			cd /mnt/gs21/scratch/borendun/IMPACTS2022/$s/equilibration/$id
			sbatch -J eq$id eq2.sh
		fi
    	done
done
