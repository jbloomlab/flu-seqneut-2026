# Text, figures, and tables for the paper about this study

## Manuscript
The manuscript source is the LaTeX in [manuscript/](manuscript), which is where edits go.
[Snakefile](Snakefile) typesets it to `results/manuscript/preprint.pdf`, pulling in the
figures and table below. Its formatting follows the lab template at
[jbloomlab/manuscript_formatting](https://github.com/jbloomlab/manuscript_formatting);
`preprint.tex` holds the preamble and title block, `body.tex` the prose and the figure
and table floats, and `references.bib` the bibliography.

It was converted once from the Word version of the Google doc it used to live in; see
[data/manuscript_docx/](data/manuscript_docx), which keeps that document and the
conversion as a record and is not part of the workflow.

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
summaries the main pipeline already writes. Each is written as LaTeX, which the manuscript inputs, and as HTML, which can be pasted
somewhere that wants the formatting inline. They are numbered exactly as the figures are:
each is written under its descriptive name and as a `Table_N` copy, from a mapping at the
top of [Snakefile](Snakefile).

## Building
Building the figures and tables needs the `seqneut-pipeline` conda environment used by
the main pipeline (see the top-level [README.md](../../README.md)). The structure figures
additionally need a `chromium-browser` on `PATH` and network access, since the pages
they render load Mol\* from a CDN. Converting the figures for the manuscript needs
`qpdf`, which makes them byte-reproducible. Typesetting the manuscript needs a TeX
distribution with XeLaTeX; that is a system-level dependency rather than a conda one, and
the rule declares the cluster's `texlive` module, so the build needs `--use-envmodules` to
load it. Run from this directory:

    snakemake -j 1 --use-envmodules
