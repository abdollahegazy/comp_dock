#!/bin/bash
#SBATCH -C v100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4:0:0 
#SBATCH --job-name=eqDUMMY_NAME

cd $SLURM_SUBMIT_DIR
#modules are loaded automatically by the NAMD module.
module use /mnt/home/vermaasj/modules
module load NAMD/3.0a9-gpu
srun namd3 +ppn 11 +ignoresharing eq.namd > eq.out
