#!/bin/bash --login

#SBATCH -C v100
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:0:0

cd $SLURM_SUBMIT_DIR
#modules are loaded automatically by the NAMD module.
module use /mnt/home/vermaasj/modules ;
module load NAMD/3.0.1-gpu
NUM=`ls system_eq[0-9][0-9][0-9].dcd 2>/dev/null | wc -l`
PRINTNUM=`printf "%03d" $NUM`
srun namd3 ++ppn 1 +devices 0 system_eq.namd > system_eq${PRINTNUM}.out 2> system_eq${PRINTNUM}.e
