#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --array=0-6

# IMPORTS
module purge
source myconda
module load snakemake/7.7.0

# RUN SCRIPT
snakemake \
    --cores all \
    --profile profile/snakemake_profile \
    --use-conda -f annotate #--unlock

#! WARNING - if the slurm-*.txt files say that there is a locked file error, then:
# - uncomment the --unlock flag above
# - run this script once
# - it will just unlock the files
# - comment it back out 
# - and run the script again   