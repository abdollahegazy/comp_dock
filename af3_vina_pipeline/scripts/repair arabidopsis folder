# !/bin/bash

for d in $(ls -f /mnt/gs18/scratch/users/borendun/IMPACTS2022/Arabidopsis/fasta/*.fasta)
do
	id=${d##/mnt/gs18/scratch/users/borendun/IMPACTS2022/$s/fasta/}
	uid=${id%%.fasta}
	if [ -f /mnt/gs18/scratch/users/borendun/IMPACTS2022/Arabidopsis/alphafold/$uid/$uid/ranked_0.pdb ]; then
        echo $uid
        echo $id
        mv /mnt/gs18/scratch/users/borendun/IMPACTS2022/Arabidopsis/alphafold/$uid/$uid/* /mnt/gs18/scratch/users/borendun/IMPACTS2022/Arabidopsis/alphafold/$uid
        rm -d /mnt/gs18/scratch/users/borendun/IMPACTS2022/Arabidopsis/alphafold/$uid/$uid/
        fi
done