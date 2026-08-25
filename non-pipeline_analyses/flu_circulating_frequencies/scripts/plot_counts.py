"""Plot weekly counts of public sequences matching a library subtype/region."""

import sys

import altair as alt
import matplotlib.colors
import pandas as pd
from mutation_categories import (
    all_categories,
    category_label,
    mutation_count,
    pretty_label,
)

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

REGION_DISPLAY = {"HA1": "HA1", "ectodomain": "HA ectodomain"}

subtype = snakemake.wildcards.subtype
region = snakemake.wildcards.region
max_mutations_computed = snakemake.params.max_mutations_computed
plot_mutation_threshold = snakemake.params.plot_mutation_threshold
gisaid_snapshot_date = snakemake.params.gisaid_snapshot_date

assert 0 <= plot_mutation_threshold <= max_mutations_computed, (
    f"plot_mutation_threshold ({plot_mutation_threshold}) must be between 0 and "
    f"max_mutations_computed ({max_mutations_computed})"
)

df = pd.read_csv(snakemake.input.counts_csv, parse_dates=["week_start"])
df = df.loc[(df["subtype"] == subtype) & (df["region"] == region)].copy()

# collapse the full-resolution categories onto the coarser plotting threshold,
# e.g. with a threshold of 2: identical/1/2 stay distinct, 3/4/>4 all become ">2
# mutations different" -- this never requires rerunning count_matches.
# `category_label` treats any n above `plot_mutation_threshold` as the catch-all.
n_mutations = df["category"].apply(lambda c: mutation_count(c, max_mutations_computed))
df["category"] = n_mutations.apply(lambda n: category_label(n, plot_mutation_threshold))
df["category_rank"] = df["category"].apply(
    lambda c: mutation_count(c, plot_mutation_threshold)
)
df = df.groupby(["week_start", "category", "category_rank"], as_index=False)[
    "n_sequences"
].sum()
df["category_label"] = df["category"].map(pretty_label)

label_order = [pretty_label(c) for c in all_categories(plot_mutation_threshold)]

# Vega-Lite's built-in `scale.scheme` for a discrete domain samples bin *centers*
# of the colormap, never its true endpoints. Sample the colormap ourselves at
# evenly spaced points, stopping just short of the very end (the last stop of
# most colormaps, e.g. inferno's pale yellow, reads as too light/washed out).
colormap = matplotlib.colormaps["inferno"]
n_shown = len(label_order)
max_extent = 0.93
category_colors = [
    matplotlib.colors.rgb2hex(colormap(i / (n_shown - 1) * max_extent))
    for i in range(n_shown)
]

chart = (
    alt.Chart(df)
    .mark_area(interpolate="step-after", opacity=1)
    .encode(
        x=alt.X(
            "week_start:T",
            title="week",
            axis=alt.Axis(
                format="%b-%Y",
                labelAngle=-90,
                tickCount={"interval": "month", "step": 1},
            ),
        ),
        y=alt.Y(
            "n_sequences:Q",
            title="number of sequences (by week)",
            stack="zero",
            axis=alt.Axis(tickCount=5, titlePadding=10),
        ),
        color=alt.Color(
            "category_label:N",
            sort=label_order,
            scale=alt.Scale(domain=label_order, range=category_colors),
            legend=alt.Legend(
                title="difference from closest library strain",
                titleLimit=300,
                symbolType="square",
            ),
        ),
        order=alt.Order("category_rank:Q"),
        tooltip=[
            alt.Tooltip("week_start:T", title="week"),
            alt.Tooltip("category_label:N", title="category"),
            alt.Tooltip("n_sequences:Q", format=".3g"),
        ],
    )
    .properties(
        width=504,
        height=240,
        title=(
            f"all human {subtype} {REGION_DISPLAY[region]} protein sequences "
            f"as of {gisaid_snapshot_date}"
        ),
    )
    .configure_title(fontSize=17)
    .configure_axis(
        titleFontSize=17, titleFontWeight="normal", labelFontSize=13, grid=False
    )
    .configure_legend(titleFontSize=16, labelFontSize=15)
)

chart.save(snakemake.output.chart_html)
chart.save(snakemake.output.chart_svg)
