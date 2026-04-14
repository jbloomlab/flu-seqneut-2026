#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1G

snakemake --profile profiles/ -s Snakefile

