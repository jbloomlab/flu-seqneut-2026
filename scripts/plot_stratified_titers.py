"""Interactive Altair plots splitting the sera by their relative titer to two strains.

Writes a standalone HTML per subtype, strain set, tree coloring, and chart type. Dropdowns
choose a reference and a comparator strain, and a slider sets a fold-change threshold; each
serum is then assigned to one of two groups by whether its comparator-to-reference titer
ratio exceeds that threshold, and the two groups are drawn across every strain in
contrasting colors. Beside the titers is the age distribution of each group.

The split is computed by `vega` transforms rather than in `pandas`, since the strains and
the threshold are chosen in the browser. Sera are pooled rather than faceted by cohort: a
faceted chart cannot be aligned with the density panel beside it, and a cohort dropdown
serves the same purpose.

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
stratified_config = snakemake.params.stratified_config

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

# the groups the sera are split into, in the order they are drawn and colored. The third
# is the sera that cannot be assigned because they lack a titer against one of the two
# selected strains; they are counted but not drawn.
GROUP_ABOVE = "above threshold"
GROUP_AT_OR_BELOW = "at or below threshold"
GROUP_UNMEASURED = "not measured against both"

# width of the age density panel drawn beside the titers
DENSITY_WIDTH = 130

# `log2` units for the threshold slider, since a fold change is read in doublings, over a
# range wide enough to put every serum on one side of the split at either end
THRESHOLD_SLIDER_MIN = -8
THRESHOLD_SLIDER_MAX = 8
THRESHOLD_SLIDER_STEP = 0.1

group_colors = stratified_config["group_colors"]
if set(group_colors) != {"above", "at_or_below"}:
    raise ValueError(
        f"`group_colors` must name exactly 'above' and 'at_or_below', "
        f"got {sorted(group_colors)}"
    )

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

# "All" first, so it is what the cohort dropdown starts on
all_cohorts = ["All"] + sorted(set(sera_multicohort["cohort"]) - {"All"})
print(f"Cohorts: {all_cohorts}")

# the cohorts each serum belongs to, looked up per serum and filtered on by the dropdown
cohorts_per_serum = sera_multicohort.groupby("serum")["cohort"].apply(list)
metadata = metadata.merge(
    cohorts_per_serum.rename("cohorts").reset_index(), on="serum", how="left"
)
sera_missing_cohorts = metadata.loc[metadata["cohorts"].isna(), "serum"].tolist()
if sera_missing_cohorts:
    raise ValueError(
        f"Some sera missing cohort assignments: {sera_missing_cohorts[:10]}"
    )

# ---- controls shared by every chart -------------------------------------------------

# a dropdown rather than the summary charts' selectable legend: a legend-bound selection
# cannot drive a selection across the panels of a concatenated chart, and with the sera
# pooled rather than faceted a multi-select legend would draw a serum twice when both
# "All" and the serum's own cohort were selected
cohort_param = alt.param(
    name="serum_cohort",
    value="All",
    bind=alt.binding_select(options=all_cohorts, name="serum cohort"),
)

min_age_slider, max_age_slider = titer_charts.age_sliders(metadata)

# the age axis of the density panel, held fixed so the two groups' curves are comparable
MAX_AGE = float(metadata["age_numeric"].max())

# sera fields looked up per serum; the lookup frame is cut to these so no chart embeds
# columns it never draws
METADATA_LOOKUP_FIELDS = [
    "cohorts",
    "serum_collection_date",
    "age",
    "age_numeric",
    "sex",
]

VALUE_TITER = titer_charts.VALUE_TITER

group_scale = alt.Scale(
    domain=[GROUP_ABOVE, GROUP_AT_OR_BELOW],
    range=[group_colors["above"], group_colors["at_or_below"]],
)

group_color = alt.Color(
    "group:N", title="comparator titer relative to reference", scale=group_scale
)

# the density panel shares the titers' legend rather than drawing a second one
density_color = alt.Color("group:N", scale=group_scale, legend=None)

GROUP_TOOLTIP = alt.Tooltip("group:N", title="serum group")
FOLD_CHANGE_TOOLTIP = alt.Tooltip(
    "fold_change:Q", title="comparator / reference titer", format=".3g"
)

# the serum's own annotations, tooltipped on its line
SERUM_TOOLTIPS = [
    GROUP_TOOLTIP,
    FOLD_CHANGE_TOOLTIP,
    alt.Tooltip(
        "median_titer_serum:Q",
        title="serum median titer",
        format=VALUE_TITER["format"],
    ),
    alt.Tooltip("serum_collection_date:N", title="serum date"),
    alt.Tooltip("age:N", title="age"),
    alt.Tooltip("sex:N"),
]

# the aggregates split by group as well, since the color is what separates the groups
GROUP_GROUPBY = [*titer_charts.VIRUS_GROUPBY, "group"]


def strain_dropdowns(axis_labels, reference, comparator):
    """Return the dropdowns naming the strains the titer ratio is taken between."""
    return tuple(
        alt.param(
            name=name,
            value=value,
            bind=alt.binding_select(options=axis_labels, name=label),
        )
        for name, value, label in [
            ("reference_strain", reference, "reference strain"),
            ("comparator_strain", comparator, "comparator strain"),
        ]
    )


def threshold_slider(fold_change_threshold):
    """Return the slider setting the fold change the sera are split at."""
    return alt.param(
        name="log2_threshold",
        value=math.log2(fold_change_threshold),
        bind=alt.binding_range(
            min=THRESHOLD_SLIDER_MIN,
            max=THRESHOLD_SLIDER_MAX,
            step=THRESHOLD_SLIDER_STEP,
            name="log2 fold change threshold",
        ),
    )


def add_lookups_and_filters(chart, chart_viruses, serum_medians, median_sliders):
    """Look up the serum and virus annotations, then filter to the sera drawn.

    Every filter here is a property of the serum rather than of one of its titers, so a
    serum is either kept whole or dropped whole and the ratio computed below still sees
    all of its strains.

    """
    return (
        chart.transform_lookup(
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
        .transform_filter(alt.expr.indexof(alt.datum["cohorts"], cohort_param) >= 0)
        .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
        .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
        .transform_filter(
            alt.datum["median_titer_serum"] >= alt.expr.pow(10, median_sliders[0])
        )
        .transform_filter(
            alt.datum["median_titer_serum"] <= alt.expr.pow(10, median_sliders[1])
        )
    )


def add_group(chart):
    """Assign each serum to a group by its comparator-to-reference titer ratio.

    The two titers the ratio is taken between are picked in the browser, so they are
    gathered onto every row of a serum with a `joinaggregate` over a field that holds the
    titer on the selected strain's row and nothing elsewhere; `vega`'s `max` skips the
    nulls. Written as `vega` expressions rather than through `alt.expr` so the group names
    are unambiguously string literals.

    """
    return (
        chart.transform_calculate(
            _reference_titer="datum.axis_label === reference_strain ? datum.titer : null",
            _comparator_titer=(
                "datum.axis_label === comparator_strain ? datum.titer : null"
            ),
        )
        .transform_joinaggregate(
            reference_titer="max(_reference_titer)",
            comparator_titer="max(_comparator_titer)",
            groupby=["serum"],
        )
        .transform_calculate(
            fold_change="datum.comparator_titer / datum.reference_titer",
            group=(
                "isValid(datum.reference_titer) && isValid(datum.comparator_titer)"
                f" ? (datum.fold_change > pow(2, log2_threshold) ? '{GROUP_ABOVE}'"
                f" : '{GROUP_AT_OR_BELOW}') : '{GROUP_UNMEASURED}'"
            ),
        )
    )


def with_group(chart, chart_viruses, serum_medians, median_sliders):
    """The lookups, filters, and group assignment every panel of a chart shares."""
    return add_group(
        add_lookups_and_filters(chart, chart_viruses, serum_medians, median_sliders)
    )


def drawn_only(chart):
    """Drop the sera that cannot be assigned to a group; the counts still report them."""
    return chart.transform_filter(f"datum.group !== '{GROUP_UNMEASURED}'")


def one_row_per_serum(chart):
    """Collapse a serum's strains to one row carrying its age and its group.

    The age distributions and the counts are per serum, so they must not weight a serum by
    how many strains it was measured against.

    """
    return chart.transform_aggregate(
        age_numeric="max(age_numeric)", groupby=["serum", "group"]
    )


def median_points(base):
    """Median titer per strain and group, as points."""
    return titer_charts.median_points(
        base, VALUE_TITER, groupby=GROUP_GROUPBY, color=group_color
    )


def serum_lines(base):
    """One line per serum across the strains, colored by the serum's group."""
    return titer_charts.serum_lines(
        base, VALUE_TITER, tooltip_extras=SERUM_TOOLTIPS, color=group_color
    )


def interquartile_range(base):
    """Shaded interquartile range of each group's titers for each strain."""
    return titer_charts.interquartile_range(
        base, VALUE_TITER, groupby=["axis_label", "group"], color=group_color
    )


CHART_TYPES = {
    "individual_sera": {
        "title": "median (points) and per-serum (lines) titers",
        "build": lambda base: serum_lines(base) + median_points(base),
    },
    "interquartile_range": {
        "title": "median (points) and interquartile range (bands) titers",
        "build": lambda base: interquartile_range(base) + median_points(base),
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


def age_density(chart):
    """Distribution of the ages of each group's subjects, beside the titers.

    Each group's density integrates to one, so the two shapes stay comparable when the
    groups are very different sizes; the counts text is what gives their sizes.

    """
    return (
        drawn_only(one_row_per_serum(chart))
        .transform_density(
            "age_numeric",
            groupby=["group"],
            as_=["age", "density"],
            counts=False,
            extent=[0, MAX_AGE],
        )
        # `line` outlines each area, which an overlapping pair of them needs to stay
        # readable; `stack=False` overlays them rather than summing one onto the other
        .mark_area(opacity=0.5, line=True)
        .encode(
            x=alt.X("age:Q", title="subject age (years)"),
            y=alt.Y(
                "density:Q",
                title="density",
                stack=False,
                axis=alt.Axis(labels=False, ticks=False),
            ),
            color=density_color,
            tooltip=[alt.Tooltip("group:N", title="serum group")],
        )
        .properties(width=DENSITY_WIDTH)
    )


def counts_readout(chart):
    """Text giving how many sera fall in each group, and how many could not be assigned.

    Computed by transforms downstream of the same filters as the plot, so it follows the
    dropdowns and the sliders; the counts are of sera, not of titers.

    """
    return (
        one_row_per_serum(chart)
        .transform_calculate(
            _above=f"datum.group === '{GROUP_ABOVE}' ? 1 : 0",
            _at_or_below=f"datum.group === '{GROUP_AT_OR_BELOW}' ? 1 : 0",
            _unmeasured=f"datum.group === '{GROUP_UNMEASURED}' ? 1 : 0",
        )
        .transform_aggregate(
            n_above="sum(_above)",
            n_at_or_below="sum(_at_or_below)",
            n_unmeasured="sum(_unmeasured)",
        )
        .transform_calculate(
            readout=(
                "'n=' + datum.n_above + ' above threshold, n=' + datum.n_at_or_below"
                " + ' at or below, n=' + datum.n_unmeasured"
                " + ' not measured against both strains'"
            )
        )
        .mark_text(align="center", fontSize=12)
        .encode(
            text="readout:N",
            x=alt.value(titer_charts.READOUT_WIDTH / 2),
            y=alt.value(7),
        )
        # no stroke: this is a line of text, not a panel, and `configure_view` boxes
        # every view in the chart
        .properties(
            width=titer_charts.READOUT_WIDTH,
            height=14,
            view=alt.ViewBackground(stroke=None),
        )
    )


def threshold_readout():
    """Text naming the strains compared and the fold change the sera are split at.

    Its own one-row frame, with none of the plot's filters, so the line still reads when
    the controls exclude every serum.

    """
    return (
        alt.Chart(pd.DataFrame({"row": [0]}))
        .transform_calculate(
            # 4 significant digits, and `~` to trim the trailing zeros that leaves
            readout="'splitting sera at a ' + format(pow(2, log2_threshold), '.4~r')"
            " + '-fold titer ratio of ' + comparator_strain + ' to ' + reference_strain"
        )
        .mark_text(align="center", fontSize=12)
        .encode(
            text="readout:N",
            x=alt.value(titer_charts.READOUT_WIDTH / 2),
            y=alt.value(7),
        )
        .properties(
            width=titer_charts.READOUT_WIDTH,
            height=14,
            view=alt.ViewBackground(stroke=None),
        )
    )


# ---- build every chart, slicing the data once per subtype and strain set -----------

for (subtype, strain_set), records in itertools.groupby(
    sorted(charts_to_make, key=lambda r: (r["subtype"], r["strain_set"])),
    key=lambda r: (r["subtype"], r["strain_set"]),
):
    records = list(records)
    chart_titers, chart_viruses, strain_order, _, serum_medians = (
        titer_charts.slice_chart_data(
            titers,
            viruses,
            subtype,
            strain_set,
            STRAIN_SETS[strain_set],
            subtype_params,
        )
    )

    # config names the strains by `derived_haplotype`; the vaccine set labels a strain
    # `"D.3.1 (cell vaccine)"`, so each is resolved to the axis label it is drawn under
    entry = stratified_config["charts"][subtype][strain_set]
    haplotype_to_label = dict(
        chart_viruses.merge(
            viruses[["virus", "derived_haplotype"]], on="virus", validate="one_to_one"
        )[["derived_haplotype", "axis_label"]].itertuples(index=False, name=None)
    )
    for key in ["reference", "comparator"]:
        if entry[key] not in haplotype_to_label:
            raise ValueError(
                f"`{key}` {entry[key]!r} is not a derived_haplotype of any {subtype} "
                f"{strain_set} strain; those are {sorted(haplotype_to_label)}"
            )
    if entry["reference"] == entry["comparator"]:
        raise ValueError(
            f"{subtype} {strain_set} `reference` and `comparator` are both "
            f"{entry['reference']!r}, but must be different strains"
        )
    if entry["fold_change_threshold"] <= 0:
        raise ValueError(
            f"{subtype} {strain_set} `fold_change_threshold` must be > 0, got "
            f"{entry['fold_change_threshold']}"
        )

    reference_param, comparator_param = strain_dropdowns(
        strain_order,
        haplotype_to_label[entry["reference"]],
        haplotype_to_label[entry["comparator"]],
    )
    log2_threshold_param = threshold_slider(entry["fold_change_threshold"])
    median_sliders = titer_charts.median_titer_sliders(serum_medians)
    print(
        f"{subtype} {strain_set}: {len(chart_viruses)} strains, {len(chart_titers)} "
        f"titers, splitting at {entry['fold_change_threshold']}-fold of "
        f"{entry['comparator']} to {entry['reference']}"
    )

    facet_size = plot_titer_summaries_params["facet_size"]
    base = titer_charts.base_chart(
        chart_titers,
        strain_order,
        facet_size,
        [
            titer_charts.virus_selection,
            titer_charts.serum_selection,
            cohort_param,
            min_age_slider,
            max_age_slider,
            reference_param,
            comparator_param,
            log2_threshold_param,
        ],
    )

    titer_base = drawn_only(
        with_group(base, chart_viruses, serum_medians, median_sliders)
    )
    # the same frame object as the titer panel, so `altair` embeds the rows only once
    side_base = with_group(
        alt.Chart(chart_titers), chart_viruses, serum_medians, median_sliders
    )
    density = age_density(side_base).properties(height=facet_size)
    counts = counts_readout(side_base)

    for record in records:
        chart_type = CHART_TYPES[record["chart_type"]]
        titer_panel = chart_type["build"](titer_base)
        title = (
            f"{chart_type['title']} for {subtype} {strain_set} strains, "
            "sera split by relative titer to two strains"
        )
        if strain_set == "recent":
            titer_panel = titer_charts.add_tree(
                titer_panel,
                tree_jsons[subtype],
                subtype_params[subtype],
                record["color_label"],
            )
            title += f", tree colored by {record['color_label']}"
        chart = titer_charts.finalize(
            alt.hconcat(titer_panel, density, spacing=14).resolve_scale(
                color="independent"
            ),
            title,
            "",
            above=[
                threshold_readout(),
                counts,
                titer_charts.median_titer_readout(*median_sliders),
            ],
        )

        print(f"Saving to {record['path']!r}")
        chart.save(record["path"])
