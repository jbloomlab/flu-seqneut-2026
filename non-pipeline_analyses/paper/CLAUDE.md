# Instructions for Claude Code

Conventions for this directory, which builds the paper's figures and tables. The lab
standards imported by the top-level [CLAUDE.md](../../CLAUDE.md) otherwise apply. A table
is built like a figure and follows every convention below, reading `FIGURE_NUMBERS` as
`TABLE_NUMBERS` and `figures/` as `tables/`.

## Name figures for what they show, not for their number

Each figure is a heavily hand-tuned one-off, so it gets its own rule in `Snakefile` and
its own script in `scripts/`, both named for what the figure shows (not `figure_1`).

Where the same figure is drawn for more than one unit — a subtype, a serum group — that
is one rule with a wildcard and one script, not near-duplicates of either. The figure
files still carry the per-unit descriptive name, so the rule and script are named for
the figure with that part left out.

Mechanics that more than one figure needs — pulling a spec out of a chart's HTML,
rendering one — go in a module in `scripts/` that the figure scripts import, declared as
a rule `input:` so editing it triggers a rerun. Only mechanics: what a given figure
keeps, drops, or relabels stays in that figure's own script.

Nothing here is named or labelled for a figure's number. `manuscript/body.tex` includes
each figure by its descriptive name and labels it `fig:<that name>`, and LaTeX assigns the
number from the order the floats fall; the text refers to them with `\Cref`. Renumbering
is therefore just moving a float, and a number written out by hand anywhere -- in a rule,
a script, a filename, a label, or the prose -- is a bug.

## A figure's configuration lives in its script

For this directory only, a figure's parameters — which panels it draws, how they are
laid out, which legends are dropped — are constants at the top of that figure's script,
commented, rather than keys in a `config.yml`. This is a deliberate exception to the
standards' "no experimental parameters in code": these are layout choices for a single
figure, not scientific parameters, and a config would just be indirection. Anything
scientific stays in the analysis that produced the panels.

## Figures must be byte-reproducible

Never write a date, timestamp, hostname, version string, or any other build metadata
into a generated figure. Unchanged panels must produce a byte-identical figure file
whenever it is rebuilt, so a rerun does not show up as a diff.

## Generated figures and the typeset PDF live outside `results/`

A figure's composed SVG, the PDF converted from it, and the table LaTeX are tracked in
`figures/` and `tables/`, and `latexmk` typesets in place to `manuscript/preprint.pdf`
with its aux files beside it. This is a deliberate exception to the standards' rule that
code writes only to `results/`: it is what lets the manuscript be typeset by hand in an
editor from a clone, with no build directory and nothing to rebuild first. Keep the aux
files out of the rule's `output:` -- `latexmk` reads them to decide how many passes a
rerun needs -- and keep the typeset PDF gitignored, since it is a large regenerable
binary. Only the structure renders' intermediates stay in `results/`.

## Compose panels, do not re-plot them

A figure reads what the other analyses already produce -- an SVG, a chart's embedded
spec, a structure viewer's page -- and composes or re-renders that. Do not re-plot from
their underlying data: that duplicates their plotting code and lets the figure drift
from the analysis it is supposed to show. Where a figure needs something the analysis
did not draw, change what a figure keeps, restyles, or overlays; do not rebuild the
thing itself.

## The manuscript is hand-edited LaTeX

`manuscript/` is source, not output: `body.tex` holds the prose and the figure and table
floats, `frontmatter.tex` the title block, `preprint.tex` the preamble. Edit them
directly. `data/manuscript_docx/` records the one-time conversion from the Word version
of the Google doc the manuscript used to live in; it is history, and re-running it would
overwrite the current LaTeX with a stale document.

Prose is the authors' to write. Change wording only when asked, and never as a side
effect of fixing layout or making something compile.

Tables and supplementary files follow the same rule as the figures: referred to with
`\Cref` against a descriptive label, never by a number written out. Supplementary files
are not floats, so they use the `\supplementaryfile` command that `preprint.tex` defines,
which keeps their counter and their label together.

Citations are `\citep` keys into `references.bib`. Add a reference by adding its BibTeX
entry, not by writing the citation into the prose.
