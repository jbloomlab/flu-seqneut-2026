#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=1:00:00

snakemake --profile profiles/ -s workflow/Snakefile
