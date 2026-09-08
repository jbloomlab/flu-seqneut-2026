# Instructions for Claude Code

Conventions for this directory, which builds the paper's figures. The lab standards
imported by the top-level [CLAUDE.md](../../CLAUDE.md) otherwise apply.

## Name figures for what they show, not for their number

Each figure is a heavily hand-tuned one-off, so it gets its own rule in `Snakefile` and
its own script in `scripts/`, both named for what the figure shows (not `figure_1`).

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

A figure reads the SVGs the other analyses already produce and composes them. Do not
re-plot from their underlying data: that duplicates their plotting code and lets the
figure drift from the analysis it is supposed to show.
