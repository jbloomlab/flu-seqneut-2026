# How well do library strains represent circulating flu?

This analysis compares the HA strains in a viral library CSV to recent publicly
available (GISAID-derived) HA sequences. For each public sequence, it finds the
closest-matching library strain and classifies the sequence by how many
mutations (substitutions, insertions, or deletions) separate it from that match:
identical, some number of mutations different, or more than that many mutations
different. This is done separately for each configured HA region and subtype,
and counts are tracked over time, binned at a configurable frequency.

The mutation-count resolution used for matching and the number of categories
shown individually in the plots (before collapsing the rest into a coarser ">N
mutations different" bucket) are separate config keys, so the plots can be
regrouped without rerunning the (slow) matching step.

## Input data

The library strain file and the public-sequence file (including its provenance
and why it isn't tracked in this repository) are both specified in
[config.yml](config.yml). See the top-level [README.md](../../README.md) for how
the library CSV itself is built.

The following input formatting requirements are hardcoded in
`scripts/count_matches.py` rather than configurable:

- The viral library CSV must have `strain`, `subtype`, and
  `protein_sequence_HA_ectodomain` columns. `subtype` values must match the
  subtypes configured in `config.yml`. All rows sharing the same `strain` (e.g.
  barcode replicates) must have an identical `protein_sequence_HA_ectodomain`.
- The public-sequence tarball must contain one FASTA file per configured
  subtype. Each FASTA header must be exactly four fields delimited by `|`:
  accession, strain name, location, and collection date, in that order, with the
  date fully specified (`YYYY-MM-DD`).

Sequences with non-canonical amino acids (e.g. `X`) are dropped before matching,
in both the library and the public sequences (this raises an error instead if
found in the library, since a library strain's sequence should never be
ambiguous).

## Configuration

Input file locations and provenance, the residue ranges defining each HA region,
the mutation-count resolution used for matching vs. plotting, the GISAID
snapshot date, and the time-binning frequency are all set in
[config.yml](config.yml). Updating to a new sequence set or viral library should
only ever require editing that file.

## Workflow

Matching every public sequence against the library is fuzzy (allows for
mismatches) and CPU-intensive, so `count_matches` is parallelized with
`multiprocessing` across cores and should be submitted to the cluster rather
than run on the login/head node:

    sbatch run_Hutch_cluster.bash

Run this from within this directory. It uses its own conda environment
([environment.yml](environment.yml)).

## Output

- [results/binned_counts.csv](results/binned_counts.csv): time-binned counts of
  public sequences per subtype/region/mutation-count category, at the full
  matching resolution.
- `results/plots/{subtype}_{region}_counts.{html,svg}`: a plot of those counts
  over time for each subtype/region combination, regrouped to the plotting
  resolution.
