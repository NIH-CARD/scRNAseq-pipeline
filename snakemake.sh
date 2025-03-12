#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:2
#SBATCH --array=0-6

module purge
module load apptainer
module load snakemake/7.7.0
module load cuda/


# Decision tree for analysis type
# 1. Single cell
# 2. Multiome -> this pipeline
# 3. Bulk

# Results for 2
# QC metric values
# Data directory
# Results directory
# Location

# Pull profile, this will only run once, and is required for running on Biowulf
git clone https://github.com/NIH-HPC/snakemake_profile.git

# Pull the containers
apptainer pull envs/snapATAC2.sif oras://quay.io/adamcatchingdti/snapatac2
apptainer pull envs/single_cell_cpu.sif oras://quay.io/adamcatchingdti/single_cell_cpu:0.2
apptainer pull envs/single_cell_cpu.sif oras://quay.io/adamcatchingdti/single_cell_gpu
apptainer pull envs/scenicplus.sif oras://hub.docker.com/litd/docker-scenicplus:latest

# Load singularity
module load singularity/4.1.5

# Bind external directories on Biowulf
. /usr/local/current/singularity/app_conf/sing_binds

# RUN SCRIPT
snakemake --cores all --profile snakemake_profile --use-singularity 

#! WARNING - if the slurm-*.txt files say that there is a locked file error, then:
# - uncomment the --unlock flag above
# - run this script once
# - it will just unlock the files
# - comment it back out 
# - and run the script again   
