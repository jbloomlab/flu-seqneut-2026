# Instructions for Claude Code
  
## Bloom lab coding standards

@bloomlab-coding-standards/CLAUDE.md

The [standards](https://github.com/jbloomlab/bloomlab-coding-standards) are pinned at a
commit of that submodule; update periodically with
`git submodule update --remote bloomlab-coding-standards`.

## `seqneut-pipeline` coding standards

Also read [seqneut-pipeline/CLAUDE.md](seqneut-pipeline/CLAUDE.md), which covers the
pipeline included here as a submodule, and follow its principles as well.

## Overrides of the lab coding standards

The `age` column in `data/sera_metadata/*.csv` holds either an exact age or an age
range, depending on what the cohort released. Both forms are accepted and documented in
`README.md`, and `scripts/aggregate_sera_metadata.py` parses them into `age_numeric`.
This is a deliberate exception to the standards' "one column, one job" rule: do not flag
it, and do not propose splitting it into separate bound columns.

## Testing rules without re-running the whole pipeline

A full run takes many hours, so when testing a new or changed rule, target only the
output files that rule makes and pass `--rerun-triggers mtime`. Check the job count with
`-n` first and stop if the dry run queues upstream rules such as `count_barcodes`.

Never use `-F` or `--forceall`: they force the entire upstream DAG, and a run killed
partway through leaves outputs deleted and Snakemake metadata marked incomplete.
