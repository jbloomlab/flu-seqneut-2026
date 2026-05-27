# README
This repository designs HA insert sequences that can be submitted to 
[GenScript](https://www.genscript.com) for gene fragment synthesis and cloning. These inserts are 
designed to fit between the BsmBI-v2 cut sites in the 
[Bloom lab vector 5764](./data/5764_pHW_SeqNeutCloningVector_WSNHAflankGFP.gb). They should begin 
after the 19th codon for H3N2 strains or 20th codon for H1N1 strains of WSN upstream signal peptide, 
and continue all the way through the end of the HA coding region, followed by a double stop codon 
and a 16-nucleotide barcode. **Note:** these scripts are specifically hardcoded to work with the 
[Bloom lab vector 5764](./data/5764_pHW_SeqNeutCloningVector_WSNHAflankGFP.gb) and thus make many
assumptions that likely will not generalize to other backbones. 

The key output files are located in `results`. This directory will contain the following: 

* `results/orders/{date}_{order}_order.csv` - This is a minimal CSV with just the shortname 
identifier and the insert sequence that can be easily uploaded to the ordering portal.
* `results/plasmid_log/{date}_{order}.csv` - This CSV contains all the information needed for 
recording the newly ordered plasmids in the seqneut plasmid log. 
* `results/genbank/` - This directory will contain individual genbank files for each plasmid.

## Configuration

The main configuration is handled through the `config.yaml` file. To add new orders, add items under
the `orders` heading. See the `testset` order as an example for required fields.  

Other important fields:

* `past_barcodes_to_avoid` - add any files specifying existing barcodes to prevent duplication
* `vaccine_annotation_file` - specify vaccine annotation info for the log. These are keyed off of 
the protein_sequence_HA_ectodomain sequence. This currently points to the file
`data/vaccine_annotations.csv` which can be added to if any new vaccine strains are added to the 
library. 

## Insert design details for H3N2 and H1N1 subtypes
In our design, our H3 HA constructs have an upstream signal peptide from WSN, an HA ectodomain 
matching a currently circulating strain, an endodomain matching a recent H3 consensus sequence, 
and a C-terminal domain from WSN

Our H1 HA constructs are slightly different, with an upstream signal peptide from WSN, an HA 
ectodomain matching a currently circulating strain, and an endodomain and a C-terminal domain also 
from WSN

## Barcode design
We specifically design our barcodes to avoid barcodes used by prior libraries and barcodes starting 
with `GG` due to sequencing issues. Prior libraries are identified in the `config.yml`. 

## Generating the output files 
To run the pipeline, build and activate the environment with:

```
conda env create -f environment.yml
conda activate constructs
```

Then this can be run using Snakemake. To submit to the SLURM scheduler at the Hutch, run with

```
sbatch run_snakemake.bash
```
