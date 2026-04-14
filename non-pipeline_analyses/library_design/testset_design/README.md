# Testset library design
This subdirectory contains a simplified Snakemake pipeline for generating the input nucleotide 
sequences needed for construct generation for a small specified set of influenza HA strains. 
This set is much smaller than the typical library and is being used as a test set for various 
vector backbone and process changes and will likely be folded into the larger library downstream.

The main input data are TSVs of the test haplotypes (for H3N2 and H1N1) and the output is a TSV with 
the selected test strains, corresponding genbank accessions, amino acid sequences, and nucleotide 
sequences.

## How to use and look at results
The input and configuration are all specified in [config.yaml](config.yaml).

The starting point are the *selected_haplotypes* in the files defined in [config.yaml](config.yaml).
Those TSV files manually select which haplotypes you want to include in the test set. There is no
additional selection logic beyond their inclusion here in these files. 

The pipeline will then validate all the selected haplotypes to identify haplotypes with ambiguous or 
missing amino acids (ideally none). The results of this validation, which you should inspect to 
check everything makes sense, are in 
[results/selected_haplotype_validation/](results/selected_haplotype_validation/).
If any of the selected haplotypes have ambiguous or missing amino acids, you will get an error!

Finally, a final TSV with all of the selected strains, the HA protein sequence, and the closest 
GenBank accession and nucleotide sequence are saved in 
[results/aggregated_library_strains/library_strains.tsv](results/aggregated_library_strains/library_strains.tsv).
This file includes columns for the representative strain name, derived haplotype name, 
HA protein sequence, closest GenBank nucleotide accession, and the matched nucleotide sequence. Use 
this as the final file for next steps.

## Workflow structure

### Configuration
All pipeline configuration is in [config.yaml](config.yaml), which specifies input data and various 
configurations and options, and should be largely self explanatory.

### Input data
All input data are in [./data/](data):

 - Selected HA haplotypes (created by Andrew Butler, April-13-2026):
   + [2026-testset-strain-selection-H1N1.tsv](data/2026-testset-strain-selection-H1N1.tsv)
   + [2026-testset-strain-selection-H3N2.tsv](data/2026-testset-strain-selection-H3N2.tsv)

### Workflow
The workflow is run by the `snakemake` pipeline in [Snakefile](Snakefile), and all results are 
placed in [./results/](results).

To run the pipeline, build and activate the `conda` environment in 
[environment.yml](environment.yml), then run the pipeline with:

    snakemake --use-conda -j <ncpus>
    
or on the Hutch cluster you can also just do:

    sbatch run_Hutch_cluster.bash

Note that the step that matches protein haplotypes with Genbank nucleotide accessions can take a 
little while to run.
