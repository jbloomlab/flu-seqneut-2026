# Instructions for Claude Code

Conventions for this directory, which builds the paper's figures. The lab standards
imported by the top-level [CLAUDE.md](../../CLAUDE.md) otherwise apply.

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

A figure's number in the manuscript comes only from the `FIGURE_NUMBERS` mapping at the
top of `Snakefile`, and the `numbered_figure` rule copies each figure to a `Figure_N.svg`
alongside it. Never name a rule, a script, or a figure's own output for its number, and
never write a number anywhere else: renumbering a figure, which happens repeatedly while
a paper is in flight, must stay a one-line edit to that mapping.

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

## Compose panels, do not re-plot them

A figure reads what the other analyses already produce -- an SVG, a chart's embedded
spec, a structure viewer's page -- and composes or re-renders that. Do not re-plot from
their underlying data: that duplicates their plotting code and lets the figure drift
from the analysis it is supposed to show. Where a figure needs something the analysis
did not draw, change what a figure keeps, restyles, or overlays; do not rebuild the
thing itself.
