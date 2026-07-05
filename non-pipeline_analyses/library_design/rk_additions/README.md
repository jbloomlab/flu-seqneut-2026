# Library additions 
This subdirectory contains a simplified Snakemake pipeline for generating the input nucleotide 
sequences needed for construct generation for a small specified set of influenza HA strains. 
This set is was recommended by Rebecca Kondor.

The main input data are TSVs of current haplotypes, and manual annotations of which haplotypes to
include.

## How to use and look at results

 The input and configuration are all specified in [config.yaml](./config.yaml).

 The starting point are the `recent_haplotypes` in the files defined in [config.yaml](./config.yaml).
 Those TSV files are used in manually selecting which haplotypes to initially consider for inclusion
 in the library. To include your choices of haplotypes, add files to the [./data][./data]
 subdirectory containing your list of haplotypes (one haplotype per line in a plain text file). Then
 update [config.yaml](./config.yaml) with these new files under the `select_recent_haplotype_files`, 
 with separate files for H1N1 and H3N2 keyed under the respective strain.

 The pipeline will then validate all the recent haplotypes to identify haplotypes with ambiguous or 
 missing amino acids. The results of this validation, which you should inspect to check that 
 everything makes sense, are in 
 [./results/recent_haplotype_validation/](./results/recent_haplotype_validation/). If any of the 
 selected haplotypes you selected for inclusion have ambiguous or missing amino acids, you will get
 an error and this will be noted in the log file! You need to either not select those haplotypes, or 
 specify them for exclusion manually under `override_select_recent_haplotypes` in 
 [config.yaml](./config.yaml). Note that you can manually fix the haplotypes so they can be added
 back, if you do not want to alter the `recent_haplotypes` files, instead create new "fixed 
 haplotypes" files and specify them for inclusion under `additional_haplotypes` in 
 [config.yaml](./config.yaml). If you have other haplotypes you definitely want to include, you can 
 also put them under `additional_haplotypes`.

 The pipeline then curates the haplotypes with results placed in 
 [./results/curated_library](./results/curated_library). This subdirectory contains TSV files listing 
 the haplotypes selected for the library and those not selected for the library. These files include
 the distances in Hamming distance and at epitope sites to strains selected in the library. The 
 Hamming distance is the number of amino acid differences across the HA1 and HA2 sequences, while the 
 epitope distances count differences only at known epitope sites (Wolf and Koel sites for H3N2;
 Caton sites for H1N1). For the selected strains (which are sorted with closest distances first),
 you want to make sure you do not have important non-selected haplotypes that are too distant from
 selected library strains. So look carefully to assess this and potentially adjust your selections 
 of library strains.

 Finally, a final TSV with all of the selected strains is saved in 
 [./results/aggregated_library_strains/library_strains.tsv](./results/aggregated_library_strains/library_strains.tsv).
 This file includes columns for the representative strain name, derived haplotype name, HA1+HA2
 protein sequence, closest GenBank accession and nucleotide sequence, and additional metadata. Use
 this as the final file for next steps.

 ## Workflow structure

 ### Configuration

 All pipeline configuration is in [config.yaml](./config.yaml), which specifies input data and 
 various configurations and options. This should be largely self explanatory.

 ### Input data

 All input data are in [./data/](./data/):

 * Recent HA haplotypes observed since December-1-2025 with metadata (created by John Huddleston, 
 May-19-2026):
    * [2026-Sept-VCM-seqneut-library-strain-selection-H1N1.tsv](./data/2026-Sept-VCM-seqneut-library-strain-selection-H1N1.tsv)
    * [2026-Sept-VCM-seqneut-library-strain-selection-H3N2.tsv](./data/2026-Sept-VCM-seqneut-library-strain-selection-H3N2.tsv)
* Manually specified haplotype selection file

### Submodules

One submodule is included; note this is included via `git submodule` at the top level:

* [match_prot_to_genbank_nt](https://github.com/jbloom/match_prot_to_genbank_nt/): 
Matches protein sequences to closest GenBank nucleotide sequences


### Workflow

The workflow is run by the `snakemake` pipeline in [Snakefile](./Snakefile), and all results are 
placed in [results/](./results/).

To run the pipeline, build and activate the `conda` environment in 
[environment.yml](./environment.yml), then run the pipeline with:

```
snakemake --use-conda -j <ncpus>
```

or run on the Hutch cluster with 

```
sbatch run_Hutch_cluster.bash
```

Note that the step that matches protein haplotypes with GenBank nucleotide accessions can take a 
little while to run.
