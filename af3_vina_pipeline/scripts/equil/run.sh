#!/bin/bash
#SBATCH -C v100
#SBATCH --array=0-10%1
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:0:0
#SBATCH --job-name=runDUMMY_NAME

cd $SLURM_SUBMIT_DIR
#modules are loaded automatically by the NAMD module.
module use /mnt/home/vermaasj/modules
module load NAMD/3.0a9-gpu
NUM=`ls *dcd | wc -l` 
PRINTNUM=`printf "%03d" $(expr $NUM - 1)`
srun namd3 ++ppn 4 +devices 0 run.namd > run${PRINTNUM}.out
