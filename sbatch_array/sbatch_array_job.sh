#!/bin/bash 
#SBATCH -p short 
#SBATCH -n 1 
#SBATCH --mem=6g
sed -n ${SLURM_ARRAY_TASK_ID}p array_tasks.list | bash
