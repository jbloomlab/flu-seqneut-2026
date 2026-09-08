# Text and figures for the paper about this study

## Manuscript
The current manuscript is [in this Google doc](https://docs.google.com/document/d/114EQK8n18Iq5AU1QisaOMNsDnG3JNO41A1PvDwvefDE/edit?tab=t.0).

## Figures
Paper figures are in [./figures/](figures), and are built by [Snakefile](Snakefile)
from the plots produced by the other analyses in [../](..). Each figure has its own rule
and its own script in [scripts/](scripts), named for what the figure shows; see
[CLAUDE.md](CLAUDE.md) for the conventions those follow.

Every figure is written into [./figures/](figures) twice: under its descriptive name, and
as a `Figure_N.svg` copy for citing in the manuscript. The mapping between the two is at
the top of [Snakefile](Snakefile), and is the only place a figure's number appears.

Building the figures needs only the `seqneut-pipeline` conda environment used by the
main pipeline (see the top-level [README.md](../../README.md)). Run from this directory:

    snakemake -j 1
