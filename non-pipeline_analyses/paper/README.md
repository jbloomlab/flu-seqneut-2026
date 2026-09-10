# Text, figures, and tables for the paper about this study

## Manuscript
The current manuscript is [in this Google doc](https://docs.google.com/document/d/114EQK8n18Iq5AU1QisaOMNsDnG3JNO41A1PvDwvefDE/edit?tab=t.0).

## Figures
Paper figures are in [./figures/](figures), and are built by [Snakefile](Snakefile) from
the plots the main pipeline and the other analyses in [../](..) already produce. Each
figure has its own rule and its own script in [scripts/](scripts), named for what the
figure shows; see [CLAUDE.md](CLAUDE.md) for the conventions those follow.

Every figure is written into [./figures/](figures) twice: under its descriptive name, and
as a `Figure_N.svg` copy for citing in the manuscript. The mapping between the two is at
the top of [Snakefile](Snakefile), and is the only place a figure's number appears.

A figure built from one of the pipeline's interactive charts is frozen at a chosen
selection state and stripped to the panels that belong in print, rendered offline with
`vl_convert` rather than exported by hand from the chart's menu. A figure built from a
`prot-struct-viz` structure page works the same way, except that the page is a Mol\*
viewer and so has to be run to be rendered: a headless browser loads it, changes the
surface, and exports the image Mol\* renders offscreen, and the labels are then drawn
over that image as vector text. Those intermediate images land in this directory's own
`results/`, which is not tracked. Those pages and charts under `results/` are not tracked
in this repository, so rebuilding such a figure requires having run the main pipeline;
the figure's own SVG is tracked, so reading the paper's figures does not.

## Tables
Paper tables are in [./tables/](tables), and are built by [Snakefile](Snakefile) from the
summaries the main pipeline already writes. They are HTML so they can be pasted into the
manuscript with their formatting intact, and are numbered exactly as the figures are:
each is written under its descriptive name and as a `Table_N.html` copy, from a mapping
at the top of [Snakefile](Snakefile).

## Building
Building the figures and tables needs the `seqneut-pipeline` conda environment used by
the main pipeline (see the top-level [README.md](../../README.md)). The structure figures
additionally need a `chromium-browser` on `PATH` and network access, since the pages
they render load Mol\* from a CDN. Run from this directory:

    snakemake -j 1
