#!/bin/bash -l
#SBATCH -t 96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=dimml002@umn.edu 
#SBATCH --job-name Phase
#SBATCH -o ./jobs/%j.Phase.out
#SBATCH -e ./jobs/%j.Phase.err
#SBATCH -p amdsmall

set -e

cd $SLURM_SUBMIT_DIR

micromamba activate phase 

snakemake -s phase.smk \
    --profile slurm.Phase \
    --rerun-triggers mtime \
    --configfile config.yaml \
    --keep-going \
    --use-conda

