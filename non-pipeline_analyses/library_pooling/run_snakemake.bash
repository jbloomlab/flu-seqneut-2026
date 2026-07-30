#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#
# Submit the library_pooling notebooks to SLURM via the Snakemake slurm profile.
# Run from this directory with: sbatch run_snakemake.bash

snakemake --profile profiles/ -s workflow/Snakefile
