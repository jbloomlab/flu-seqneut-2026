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
identifier and the insert sequence that can be easily uploaded to the ordering portal. For 
constructs being ordered in plate format, this will also include a platemap. 
* `results/plasmid_log/{date}_{order}.csv` - This CSV contains all the information needed for 
recording the newly ordered plasmids in the seqneut plasmid log. 
* `results/genbank/` - This directory will contain individual genbank files for each plasmid.
* `results/excluded_from_order/{date}_{order}_excluded.csv` - This CSV contains all the entries that
were excluded from the ordering and logging CSVs because the same HA ectodomain sequence already 
exists either in an already obtained plasmid or within the order set. 
* `results/final_library/final_library.csv` - This CSV contains the final set of sequences along 
with their associated metadata that are to be included in the final library.

## Configuration

The main configuration is handled through the `config.yaml` file. To add new orders, add items under
the `orders` heading. See the `testset` order as an example for required fields.  

Other important fields:

* `past_barcodes_to_avoid` - add any files specifying existing barcodes to prevent duplication
* `annotation_file` - specify vaccine annotation info for the log. These are keyed off of 
the protein_sequence_HA_ectodomain sequence. This currently points to the file
`data/vaccine_annotations.csv` which can be added to if any new vaccine strains are added to the 
library. 
* `past_protein_sequences_to_avoid` - add any files specifying existing HA ectodomain sequences to 
prevent ordering duplicate/identical plasmids. 

## Insert design details for H3N2 and H1N1 subtypes
In our design, our H3 HA constructs have an upstream signal peptide from WSN, an HA ectodomain 
matching a currently circulating strain, an endodomain matching a recent H3 consensus sequence, 
and a C-terminal domain from WSN

Our H1 HA constructs are slightly different, with an upstream signal peptide from WSN, an HA 
ectodomain matching a currently circulating strain, and an endodomain and a C-terminal domain also 
from WSN

## Barcode design
We specifically design our barcodes to avoid barcodes used by prior libraries and barcodes starting 
with `GG` due to sequencing issues. Prior libraries are identified in the `config.yaml`. 

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
