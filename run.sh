#!/bin/bash
#SBATCH -p medium 
#SBATCH -n 20
#SBATCH -N 1 
#SBATCH --mem=60g
#SBATCH -o log
cat commands | parallel -j20 
