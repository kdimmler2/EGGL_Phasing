#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dimml002@umn.edu
#SBATCH --job-name subset.slurm
#SBATCH -o %j.subset.out
#SBATCH -e %j.subset.err
#SBATCH -p agsmall,aglarge

bcftools view -Oz -o reference_panel.vcf.gz -R /panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/100bp_UpDown_targets.bed goldenPath.Phased.vcf.gz 
