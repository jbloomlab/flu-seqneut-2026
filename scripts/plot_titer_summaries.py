"""Interactive Altair plots summarizing the titers of one serum group.

Writes a standalone HTML per subtype, strain set, tree coloring, and chart type. The
recent-strain charts are drawn alongside the
`nextstrain-prot-titers-tree` tree for their subtype, which supplies the strain order,
the strain labels, and the coloring. The vaccine-strain charts have no tree: they are
ordered by collection date and their labels name the vaccine type.

Chart and tree are joined on `derived_haplotype`, which also names the strain a dashed
reference line is optionally drawn at.

Each chart facets by serum cohort, and a serum appears in every cohort it belongs to.

"""

import itertools
import math
import sys

import altair as alt
import pandas as pd
import titer_charts

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

_ = alt.data_transformers.disable_max_rows()

recent_vaccine_strains = snakemake.params.recent_vaccine_strains
circulating_strain_type = snakemake.params.circulating_strain_type
plot_titer_summaries_params = snakemake.params.plot_titer_summaries_params
subtypes = snakemake.params.subtypes

subtype_params = plot_titer_summaries_params["subtype_params"]

# the rule lists one tree per subtype, in `subtypes` order
tree_jsons = dict(zip(subtypes, snakemake.input.trees, strict=True))

# the rule declares its outputs in the same order as `charts`, so each record's chart is
# saved to the output path it is paired with
charts_to_make = [
    dict(record, path=path)
    for record, path in zip(
        snakemake.params.charts, snakemake.output.chart_htmls, strict=True
    )
]

STRAIN_SETS = titer_charts.strain_sets(circulating_strain_type)

titers, metadata, sera_multicohort, viruses = titer_charts.load_and_validate(
    snakemake.input.titers_csv,
    snakemake.input.sera_csv,
    snakemake.input.sera_multicohort_csv,
    snakemake.input.viruses_csv,
    recent_vaccine_strains=recent_vaccine_strains,
    circulating_strain_type=circulating_strain_type,
    subtypes=subtypes,
    subtype_params=subtype_params,
)

if "All" not in sera_multicohort["cohort"].values:
    raise ValueError(
        "Expected 'All' cohort in sera_multicohort but not found. "
        f"Available cohorts: {sera_multicohort['cohort'].unique().tolist()}"
    )

# Get list of all cohorts (for legend), with "All" first
all_cohorts = sera_multicohort["cohort"].unique().tolist()
all_cohorts = ["All"] + sorted([c for c in all_cohorts if c != "All"])
print(f"Cohorts: {all_cohorts}")

# a name that is not a real cohort would select nothing and silently render an empty chart
initial_cohorts = plot_titer_summaries_params["initial_cohorts"]
if initial_cohorts is not None:
    if not initial_cohorts:
        raise ValueError("`initial_cohorts` is empty; use null to show every cohort")
    unknown = [c for c in initial_cohorts if c not in all_cohorts]
    if unknown:
        raise ValueError(
            f"`initial_cohorts` names cohorts not in the data: {unknown}; "
            f"cohorts present are {all_cohorts}"
        )
print(f"Cohorts shown initially: {initial_cohorts or all_cohorts}")

# Aggregate cohorts into a list per serum and add to metadata
cohorts_per_serum = sera_multicohort.groupby("serum")["cohort"].apply(list)
metadata = metadata.merge(
    cohorts_per_serum.rename("cohorts").reset_index(),
    on="serum",
    how="left",
)
sera_missing_cohorts = metadata.loc[metadata["cohorts"].isna(), "serum"].tolist()
if sera_missing_cohorts:
    raise ValueError(
        f"Some sera missing cohort assignments: {sera_missing_cohorts[:10]}"
    )

# ---- selections and sliders, shared by every chart --------------------------------

cohort_selection = alt.selection_point(
    fields=["cohort"],
    bind="legend",
    empty="all",
    toggle="true",
    clear=False,
    # no `value` when `initial_cohorts` is null: an empty selection with `empty="all"`
    # is what already shows every cohort
    **(
        {}
        if initial_cohorts is None
        else {"value": [{"cohort": cohort} for cohort in initial_cohorts]}
    ),
)

min_age_slider, max_age_slider = titer_charts.age_sliders(metadata)

titer_cutoff = plot_titer_summaries_params["titer_cutoff"]

# the slider's floor is the lowest titer measured, put on the same grid as its step
cutoff_slider_step = 5
cutoff_slider_min = cutoff_slider_step * round(
    titers["titer"].min() / cutoff_slider_step
)
print(f"Using {cutoff_slider_min=} and initial {titer_cutoff=}")

titer_cutoff_slider = alt.param(
    value=titer_cutoff,
    bind=alt.binding_range(
        min=cutoff_slider_min,
        max=1000,
        step=cutoff_slider_step,
        name="fraction sera below this cutoff",
    ),
)

# `log10` units for the median-titer sliders; one step is ~12% in titer, fine enough that
# the slider never steps over a serum
MEDIAN_SLIDER_STEP = 0.05

# the slider ends reach 1% past the medians they bound, so `pow(10, slider)` at an end
# cannot land a hair inside the data and drop the lowest- or highest-median serum
MEDIAN_SLIDER_PAD = 1.01

# significant digits the slider ends are rounded to, so the widget reads out a short
# number rather than the full `log10` of a titer
MEDIAN_SLIDER_DIGITS = 3


def _round_outward(x, round_):
    """Return `x` rounded to `MEDIAN_SLIDER_DIGITS` significant digits by `round_`.

    `round_` is `math.floor` for a lower bound and `math.ceil` for an upper one, so the
    result never moves inward and so can never exclude a value `x` includes.

    """
    if x == 0:
        return 0.0
    scale = 10 ** (MEDIAN_SLIDER_DIGITS - 1 - math.floor(math.log10(abs(x))))
    return round_(x * scale) / scale


def median_titer_sliders(serum_medians):
    """Return the sliders bounding which sera are drawn, by their median titer.

    Both are in `log10` titer units, as the medians span orders of magnitude, and each
    starts at an end of the observed range, so nothing is filtered until one is moved.

    """
    medians = serum_medians["median_titer_serum"]
    lo = _round_outward(math.log10(medians.min() / MEDIAN_SLIDER_PAD), math.floor)
    hi = _round_outward(math.log10(medians.max() * MEDIAN_SLIDER_PAD), math.ceil)
    return tuple(
        alt.param(
            value=value,
            bind=alt.binding_range(min=lo, max=hi, step=MEDIAN_SLIDER_STEP, name=name),
        )
        for value, name in [
            (lo, "minimum serum median titer (log10)"),
            (hi, "maximum serum median titer (log10)"),
        ]
    )


VALUE_TITER = titer_charts.VALUE_TITER
VALUE_FOLD_CHANGE = {
    "field": "fold_change",
    "title": ["titer fold change", "from serum's median"],
    "format": ".3g",
}

# dummy chart to bind the selectable legend for serum cohort
cohort_legend = (
    alt.Chart(pd.DataFrame({"cohort": all_cohorts}))
    .add_params(cohort_selection)
    .mark_point(opacity=0)
    .encode(
        fill=alt.Fill(
            "cohort",
            title="serum cohort (click to select)",
            scale=alt.Scale(domain=all_cohorts, range=["gray"]),
            legend=alt.Legend(symbolStrokeColor="black", symbolOpacity=1, columns=6),
        )
    )
    .properties(width=1, height=1)  # tiny plot; legend renders outside
)


# width of the readout view, which `alt.vconcat(center=True)` then centers under the
# chart title; the text is centered within it
READOUT_WIDTH = 300


def median_titer_readout(min_median_slider, max_median_slider):
    """Text naming the titers the median-titer sliders are currently set to.

    The sliders read out in `log10` units, so the titers those correspond to are drawn
    here instead. Its own one-row frame keeps this mark out of the faceted chart's
    filters, so the line still reads when the sliders exclude every serum, and declaring
    the sliders here rather than on the plot lifts them to the top level of the
    concatenated chart, where they are in scope for the filters.

    """
    return (
        alt.Chart(pd.DataFrame({"row": [0]}))
        .add_params(min_median_slider, max_median_slider)
        .transform_calculate(
            # 4 significant digits, and `~` to trim the trailing zeros that leaves; a
            # plain `g` would read out a titer over 1000 in scientific notation
            readout="showing sera with median titer "
            + alt.expr.format(alt.expr.pow(10, min_median_slider), ".4~r")
            + " to "
            + alt.expr.format(alt.expr.pow(10, max_median_slider), ".4~r")
        )
        .mark_text(align="center", fontSize=12)
        .encode(text="readout:N", x=alt.value(READOUT_WIDTH / 2), y=alt.value(7))
        # no stroke: this is a line of text, not a panel, and `configure_view` boxes
        # every view in the chart
        .properties(
            width=READOUT_WIDTH, height=14, view=alt.ViewBackground(stroke=None)
        )
    )


# sera fields looked up per serum; the lookup frame is cut to these so no chart embeds
# columns it never draws
METADATA_LOOKUP_FIELDS = [
    "cohorts",
    "serum_collection_date",
    "age",
    "age_numeric",
    "sex",
]

# the aggregates split by cohort as well, since the facet is what separates the cohorts
virus_groupby = [*titer_charts.VIRUS_GROUPBY, "cohort"]

# the titer stays in the tooltip even when the chart plots the fold change
FOLD_CHANGE_AGGREGATE_EXTRAS = {"median_titer": "median(titer)"}
FOLD_CHANGE_TOOLTIP_EXTRAS = [
    alt.Tooltip("median_titer:Q", format=VALUE_TITER["format"])
]

# the serum's own annotations, tooltipped on its line
SERUM_TOOLTIPS = [
    alt.Tooltip(
        "median_titer_serum:Q",
        title="serum median titer",
        format=VALUE_TITER["format"],
    ),
    alt.Tooltip("serum_collection_date:N", title="serum date"),
    alt.Tooltip("age:N", title="age"),
    alt.Tooltip("sex:N"),
]


def median_points(base, value):
    """Median of the plotted value per strain and cohort, as points."""
    fold_change = value is VALUE_FOLD_CHANGE
    return titer_charts.median_points(
        base,
        value,
        groupby=virus_groupby,
        aggregate_extras=FOLD_CHANGE_AGGREGATE_EXTRAS if fold_change else None,
        tooltip_extras=FOLD_CHANGE_TOOLTIP_EXTRAS if fold_change else (),
    )


def serum_lines(base, value):
    """One line per serum across the strains."""
    titer_tooltip = (
        [alt.Tooltip("titer:Q", format=VALUE_TITER["format"])]
        if value is VALUE_FOLD_CHANGE
        else []
    )
    return titer_charts.serum_lines(
        base, value, tooltip_extras=[*titer_tooltip, *SERUM_TOOLTIPS]
    )


def interquartile_range(base, value):
    """Shaded interquartile range of the plotted value for each strain."""
    fold_change = value is VALUE_FOLD_CHANGE
    return titer_charts.interquartile_range(
        base,
        value,
        groupby=["axis_label"],
        aggregate_extras=FOLD_CHANGE_AGGREGATE_EXTRAS if fold_change else None,
        tooltip_extras=FOLD_CHANGE_TOOLTIP_EXTRAS if fold_change else (),
    )


def frac_below_cutoff(base):
    """Fraction of sera with a titer below the slider's cutoff, as bars."""
    return (
        base.add_params(titer_cutoff_slider)
        .transform_calculate(below_cutoff=alt.datum["titer"] < titer_cutoff_slider)
        .transform_aggregate(
            n_below_cutoff="sum(below_cutoff)",
            n_total="distinct(serum)",
            groupby=virus_groupby,
        )
        .transform_calculate(
            frac_below_cutoff=alt.datum["n_below_cutoff"] / alt.datum["n_total"]
        )
        .encode(
            y=alt.Y("frac_below_cutoff:Q", title="fraction below cutoff"),
            tooltip=[
                *titer_charts.virus_tooltips,
                alt.Tooltip("frac_below_cutoff:Q", format=".2f"),
            ],
            color=alt.condition(
                titer_charts.virus_selection, alt.value("red"), alt.value("black")
            ),
        )
        .mark_bar(opacity=0.8)
    )


def reference_line(chart_titers, ref_axis_label, chart_type, value):
    """Dashed rule at `ref_axis_label`'s value within each cohort facet.

    Built from the same frame object as the rest of the layer so `altair` hoists the data
    to the layer, which puts this mark downstream of the facet's cohort and age filters:
    the line is the reference strain's own value in that panel and follows the sliders.
    `titer_cutoff_slider` is referenced but not re-added, as a param may be declared only
    once.

    """
    line = alt.Chart(chart_titers).transform_filter(
        alt.datum["axis_label"] == ref_axis_label
    )
    if chart_type == "frac_below_cutoff":
        line = (
            line.transform_calculate(
                below_cutoff=alt.datum["titer"] < titer_cutoff_slider
            )
            .transform_aggregate(
                n_below_cutoff="sum(below_cutoff)", n_total="distinct(serum)"
            )
            .transform_calculate(
                ref_value=alt.datum["n_below_cutoff"] / alt.datum["n_total"]
            )
        )
        encoding = alt.Y("ref_value:Q", title="fraction below cutoff")
    else:
        line = line.transform_aggregate(ref_value=f"median({value['field']})")
        encoding = alt.Y(
            "ref_value:Q", title=value["title"], scale=titer_charts.titer_scale
        )
    # a burnt orange dark enough to read over the interquartile band, and far enough
    # from the red of a hovered point not to be mistaken for one
    return line.encode(y=encoding).mark_rule(
        color="#D95F02", strokeWidth=2, strokeDash=[4, 3]
    )


# each chart type names its title, how it is built from the mark builders, and the value
# it plots on the titer axis
CHART_TYPES = {
    "individual_sera": {
        "title": "median (points) and per-serum (lines) titers",
        "build": lambda base, value: serum_lines(base, value)
        + median_points(base, value),
        "value": VALUE_TITER,
    },
    "interquartile_range": {
        "title": "median (points) and interquartile range titers",
        "build": lambda base, value: interquartile_range(base, value)
        + median_points(base, value),
        "value": VALUE_TITER,
    },
    "frac_below_cutoff": {
        "title": "fraction sera below titer cutoff",
        "build": lambda base, value: frac_below_cutoff(base),
        "value": VALUE_TITER,
    },
    "individual_sera_fold_change": {
        "title": (
            "median (points) and per-serum (lines) fold change from the serum's median"
        ),
        "build": lambda base, value: serum_lines(base, value)
        + median_points(base, value),
        "value": VALUE_FOLD_CHANGE,
    },
    "interquartile_range_fold_change": {
        "title": "median (points) and interquartile range fold change from the serum's median",
        "build": lambda base, value: interquartile_range(base, value)
        + median_points(base, value),
        "value": VALUE_FOLD_CHANGE,
    },
}

# the rule names the charts to make, and only these are implemented here
unknown_chart_types = {r["chart_type"] for r in charts_to_make} - set(CHART_TYPES)
if unknown_chart_types:
    raise ValueError(
        f"unknown chart_type(s) {sorted(unknown_chart_types)}; "
        f"implemented chart types are {sorted(CHART_TYPES)}"
    )
unknown_strain_sets = {r["strain_set"] for r in charts_to_make} - set(STRAIN_SETS)
if unknown_strain_sets:
    raise ValueError(
        f"unknown strain_set(s) {sorted(unknown_strain_sets)}; "
        f"defined strain sets are {sorted(STRAIN_SETS)}"
    )

# a fold change is relative to the serum's median over the strains the chart draws, which
# is only a baseline worth plotting against for the recent strains
misplaced_fold_change = sorted(
    {
        (r["chart_type"], r["strain_set"])
        for r in charts_to_make
        if CHART_TYPES[r["chart_type"]]["value"] is VALUE_FOLD_CHANGE
        and r["strain_set"] != "recent"
    }
)
if misplaced_fold_change:
    raise ValueError(
        f"fold-change charts are only made for the recent strains: {misplaced_fold_change}"
    )


def facet_and_add_lookups(
    chart, chart_viruses, serum_medians, min_median_slider, max_median_slider
):
    """Facet `chart` by cohort and look up the serum and virus annotations.

    Scoping when layering and faceting charts with `transform_lookup` requires the
    faceting to be done before the lookups, so both happen here. The fold change and the
    median-titer filters follow the lookup that brings in the serum's median, which also
    puts them upstream of the per-cohort counts in the facet labels: those fall as the
    sliders exclude sera.

    """
    return (
        chart.facet(row=alt.Row("cohort_n:N", title=None))
        .transform_lookup(
            lookup="serum",
            from_=alt.LookupData(
                data=metadata[["serum", *METADATA_LOOKUP_FIELDS]],
                key="serum",
                fields=METADATA_LOOKUP_FIELDS,
            ),
        )
        .transform_lookup(
            lookup="axis_label",
            from_=alt.LookupData(
                data=chart_viruses,
                key="axis_label",
                fields=["virus", "strain_type", "subclade"],
            ),
        )
        .transform_lookup(
            lookup="serum",
            from_=alt.LookupData(
                data=serum_medians, key="serum", fields=["median_titer_serum"]
            ),
        )
        .transform_calculate(
            fold_change=alt.datum["titer"] / alt.datum["median_titer_serum"]
        )
        # flatten cohorts list (from sera_multicohort) to one row per cohort
        .transform_flatten(["cohorts"], as_=["cohort"])
        # filter by cohort, age, and the serum's median titer
        .transform_filter(cohort_selection)
        .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
        .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
        .transform_filter(
            alt.datum["median_titer_serum"] >= alt.expr.pow(10, min_median_slider)
        )
        .transform_filter(
            alt.datum["median_titer_serum"] <= alt.expr.pow(10, max_median_slider)
        )
        # make facet labels w n per cohort
        .transform_joinaggregate(n_per_cohort="distinct(serum)", groupby=["cohort"])
        .transform_calculate(
            cohort_n="datum.cohort + ' (n=' + datum.n_per_cohort + ')'"
        )
    )


# ---- build every chart, slicing the data once per subtype and strain set -----------

for (subtype, strain_set), records in itertools.groupby(
    sorted(charts_to_make, key=lambda r: (r["subtype"], r["strain_set"])),
    key=lambda r: (r["subtype"], r["strain_set"]),
):
    chart_titers, chart_viruses, strain_order, ref_axis_label, serum_medians = (
        titer_charts.slice_chart_data(
            titers,
            viruses,
            subtype,
            strain_set,
            STRAIN_SETS[strain_set],
            subtype_params,
        )
    )
    medians = serum_medians["median_titer_serum"]
    print(
        f"{subtype} {strain_set}: {len(chart_viruses)} strains, {len(chart_titers)} "
        f"titers, serum median titers {medians.min()} to {medians.max()}"
    )
    base = titer_charts.base_chart(
        chart_titers,
        strain_order,
        plot_titer_summaries_params["facet_size"],
        [
            titer_charts.virus_selection,
            titer_charts.serum_selection,
            cohort_selection,
            min_age_slider,
            max_age_slider,
        ],
    )
    min_median_slider, max_median_slider = median_titer_sliders(serum_medians)

    for record in records:
        chart_type = CHART_TYPES[record["chart_type"]]
        value = chart_type["value"]
        layer = chart_type["build"](base, value)
        subtitle = ""
        if ref_axis_label is not None:
            # layered last so the thin line draws over the interquartile band
            layer += reference_line(
                chart_titers, ref_axis_label, record["chart_type"], value
            )
            subtitle = f"dashed orange line marks {ref_axis_label}"
        chart = facet_and_add_lookups(
            layer, chart_viruses, serum_medians, min_median_slider, max_median_slider
        )
        title = f"{chart_type['title']} for {subtype} {strain_set} strains"
        if strain_set == "recent":
            chart = titer_charts.add_tree(
                chart,
                tree_jsons[subtype],
                subtype_params[subtype],
                record["color_label"],
            )
            title += f", tree colored by {record['color_label']}"
        chart = titer_charts.finalize(
            chart,
            title,
            subtitle,
            above=[median_titer_readout(min_median_slider, max_median_slider)],
            below=[cohort_legend],
        )

        print(f"Saving to {record['path']!r}")
        chart.save(record["path"])
