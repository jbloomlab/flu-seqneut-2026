"""Interactive Altair plots comparing paired pre- and post-vaccination titers.

Writes a standalone HTML per subtype, strain set, tree coloring, and chart type. Two
comparisons are drawn: the pre- and post-vaccination titers overlaid in two colors, and
each subject's fold change from its pre- to its post-vaccination titer. Each facets by
the compared vaccination arm, and the recent-strain charts are drawn alongside their
subtype's tree as the summary charts are.

Only subjects whose pre- and post-vaccination sera both passed QC are drawn, so the
overlaid titers and the fold changes describe the same people.

"""

import itertools
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
pre_post_config = snakemake.params.pre_post_config

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

# the two sera compared, in the order they are drawn and colored; the fold change is of
# the second relative to the first
CONDITIONS = ["pre", "post"]

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

# ---- pair the sera being compared -------------------------------------------------

comparisons = pre_post_config["comparisons"]
pair_by = pre_post_config["pair_by"]
colors = pre_post_config["colors"]

if set(colors) != set(CONDITIONS):
    raise ValueError(f"`colors` must name exactly {CONDITIONS}, got {sorted(colors)}")

if pair_by not in sera_multicohort.columns:
    raise ValueError(
        f"`pair_by` {pair_by!r} is not a column of the sera metadata; columns are "
        f"{sorted(sera_multicohort.columns)}"
    )

available_cohorts = set(sera_multicohort["cohort"])

# one row per serum drawn, naming the comparison it is in and whether it is the pre- or
# the post-vaccination serum of its subject
paired_sera = []
for comparison, cohorts in comparisons.items():
    if set(cohorts) != set(CONDITIONS):
        raise ValueError(
            f"comparison {comparison!r} must name exactly {CONDITIONS}, "
            f"got {sorted(cohorts)}"
        )
    if cohorts["pre"] == cohorts["post"]:
        raise ValueError(
            f"comparison {comparison!r} compares cohort {cohorts['pre']!r} with itself"
        )

    subjects = {}
    for condition in CONDITIONS:
        cohort = cohorts[condition]
        if cohort not in available_cohorts:
            raise ValueError(
                f"comparison {comparison!r} names cohort {cohort!r}, which is not in "
                f"the sera metadata; cohorts are {sorted(available_cohorts)}"
            )
        cohort_sera = sera_multicohort[sera_multicohort["cohort"] == cohort]

        # a serum being compared with no value to pair it by is a gap in the sera
        # metadata, not a serum to quietly leave out
        unpaired = cohort_sera.loc[cohort_sera[pair_by].isnull(), "serum"].tolist()
        if unpaired:
            raise ValueError(
                f"sera in cohort {cohort!r} with no {pair_by}: {sorted(unpaired)}"
            )

        per_subject = cohort_sera.groupby(pair_by)["serum"].nunique()
        repeated = per_subject[per_subject > 1]
        if len(repeated):
            raise ValueError(
                f"subjects with more than one serum in cohort {cohort!r}:\n"
                f"{repeated.to_string()}"
            )

        subjects[condition] = cohort_sera.set_index(pair_by)["serum"]

    drawn = subjects["pre"].index.intersection(subjects["post"].index)
    if not len(drawn):
        raise ValueError(
            f"comparison {comparison!r} has no subject with both a "
            f"{cohorts['pre']!r} and a {cohorts['post']!r} serum"
        )
    for condition in CONDITIONS:
        dropped = sorted(set(subjects[condition].index) - set(drawn))
        if dropped:
            print(
                f"{comparison}: dropping {len(dropped)} subject(s) with only a "
                f"{condition}-vaccination serum: {dropped}"
            )
        paired_sera.append(
            pd.DataFrame(
                {
                    "serum": subjects[condition][drawn].values,
                    "subject": drawn,
                    "comparison": comparison,
                    "condition": condition,
                }
            )
        )
    print(f"{comparison}: drawing {len(drawn)} paired subjects")

paired_sera = pd.concat(paired_sera, ignore_index=True)

# a serum in two comparisons would be drawn twice over, and its subject counted twice
if not paired_sera["serum"].is_unique:
    raise ValueError(
        "these sera are in more than one comparison: "
        f"{paired_sera.loc[paired_sera['serum'].duplicated(), 'serum'].tolist()}"
    )

# ---- what the charts draw ---------------------------------------------------------

VALUE_FOLD_CHANGE = {
    "field": "fold_change",
    "title": ["titer fold change", "post- vs pre-vaccination"],
    "format": ".3g",
}

condition_color = alt.Color(
    "condition:N",
    title="vaccination",
    scale=alt.Scale(domain=CONDITIONS, range=[colors[c] for c in CONDITIONS]),
    legend=alt.Legend(orient="left", columns=1),
)

# sera fields looked up per serum; the lookup frame is cut to these so no chart embeds
# columns it never draws
METADATA_LOOKUP_FIELDS = ["serum_collection_date", "age", "age_numeric", "sex"]

SUBJECT_TOOLTIP = alt.Tooltip("subject:N", title=pair_by)

# the serum's own annotations, tooltipped on its line
SERUM_TOOLTIPS = [
    alt.Tooltip("serum_collection_date:N", title="serum date"),
    alt.Tooltip("age:N", title="age"),
    alt.Tooltip("sex:N"),
]

# the titers behind a fold change, and whether either was a bound rather than a
# measurement, which a ratio of the two cannot show on its own
FOLD_CHANGE_AGGREGATE_EXTRAS = {
    "median_pre_titer": "median(pre_titer)",
    "median_post_titer": "median(post_titer)",
}
FOLD_CHANGE_TOOLTIP_EXTRAS = [
    alt.Tooltip("median_pre_titer:Q", title="median pre titer", format=".1f"),
    alt.Tooltip("median_post_titer:Q", title="median post titer", format=".1f"),
]
FOLD_CHANGE_SERUM_TOOLTIPS = [
    SUBJECT_TOOLTIP,
    alt.Tooltip("pre_titer:Q", title="pre titer", format=".1f"),
    alt.Tooltip("pre_titer_bound:N", title="pre titer bound"),
    alt.Tooltip("post_titer:Q", title="post titer", format=".1f"),
    alt.Tooltip("post_titer_bound:N", title="post titer bound"),
    *SERUM_TOOLTIPS,
]

OVERLAY_SERUM_TOOLTIPS = [
    SUBJECT_TOOLTIP,
    alt.Tooltip("condition:N", title="vaccination"),
    alt.Tooltip("titer_bound:N", title="titer bound"),
    *SERUM_TOOLTIPS,
]

# the aggregates split by the compared arm as well, since the facet is what separates it
COMPARISON_GROUPBY = [*titer_charts.VIRUS_GROUPBY, "comparison"]


def overlay_median_points(base):
    """Median pre- and post-vaccination titer per strain, as points."""
    return titer_charts.median_points(
        base,
        titer_charts.VALUE_TITER,
        groupby=[*COMPARISON_GROUPBY, "condition"],
        color=condition_color,
    )


def overlay_serum_lines(base):
    """One line per serum, colored by whether it is the pre- or post-vaccination one."""
    return titer_charts.serum_lines(
        base,
        titer_charts.VALUE_TITER,
        tooltip_extras=OVERLAY_SERUM_TOOLTIPS,
        color=condition_color,
    )


def overlay_interquartile_range(base):
    """Interquartile range of the pre- and of the post-vaccination titers."""
    return titer_charts.interquartile_range(
        base,
        titer_charts.VALUE_TITER,
        groupby=["axis_label", "condition"],
        color=condition_color,
    )


def fold_change_median_points(base):
    """Median fold change per strain, as points."""
    return titer_charts.median_points(
        base,
        VALUE_FOLD_CHANGE,
        groupby=COMPARISON_GROUPBY,
        aggregate_extras=FOLD_CHANGE_AGGREGATE_EXTRAS,
        tooltip_extras=FOLD_CHANGE_TOOLTIP_EXTRAS,
    )


def fold_change_serum_lines(base):
    """One line per subject, keyed on its post-vaccination serum."""
    return titer_charts.serum_lines(
        base, VALUE_FOLD_CHANGE, tooltip_extras=FOLD_CHANGE_SERUM_TOOLTIPS
    )


def fold_change_interquartile_range(base):
    """Interquartile range of the fold changes for each strain."""
    return titer_charts.interquartile_range(
        base,
        VALUE_FOLD_CHANGE,
        groupby=["axis_label"],
        aggregate_extras=FOLD_CHANGE_AGGREGATE_EXTRAS,
        tooltip_extras=FOLD_CHANGE_TOOLTIP_EXTRAS,
    )


def no_change_line(chart_data):
    """Rule at a fold change of one, where a subject's titer did not change.

    Built from the same frame object as the rest of the layer so `altair` hoists the data
    to the layer, which puts this mark inside each facet alongside the marks it marks the
    baseline of.

    """
    return (
        alt.Chart(chart_data)
        .encode(y=alt.datum(1))
        .mark_rule(color="#888888", strokeWidth=1, strokeDash=[4, 3])
    )


# each chart type names its title, whether it draws the fold changes or the titers
# themselves, and how it is built from the mark builders
CHART_TYPES = {
    "individual_sera": {
        "title": (
            "median (points) and per-serum (lines) titers before and after vaccination"
        ),
        "fold_change": False,
        "build": lambda base: overlay_serum_lines(base) + overlay_median_points(base),
    },
    "interquartile_range": {
        "title": (
            "median (points) and interquartile range (bands) titers before and after "
            "vaccination"
        ),
        "fold_change": False,
        "build": lambda base: overlay_interquartile_range(base)
        + overlay_median_points(base),
    },
    "individual_sera_fold_change": {
        "title": (
            "median (points) and per-subject (lines) fold change from before to after "
            "vaccination"
        ),
        "fold_change": True,
        "build": lambda base: fold_change_serum_lines(base)
        + fold_change_median_points(base),
    },
    "interquartile_range_fold_change": {
        "title": (
            "median (points) and interquartile range (bands) fold change from before to "
            "after vaccination"
        ),
        "fold_change": True,
        "build": lambda base: fold_change_interquartile_range(base)
        + fold_change_median_points(base),
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

min_age_slider, max_age_slider = titer_charts.age_sliders(metadata)


def paired_chart_data(chart_titers, chart_viruses):
    """Return the drawn titers and the fold changes of the sera being compared.

    The titers arrive sliced to the strains the chart draws but covering every serum, so
    that the strain checks in `slice_chart_data` see all of them; here they are cut to
    the paired sera. `titer_bound` comes along so a tooltip can say that a titer, or a
    fold change resting on one, is a bound rather than a measurement.

    """
    bounds = titers[["serum", "virus", "titer_bound"]].merge(
        chart_viruses[["virus", "axis_label"]], on="virus", validate="many_to_one"
    )
    overlay = (
        chart_titers.merge(
            bounds[["serum", "axis_label", "titer_bound"]],
            on=["serum", "axis_label"],
            validate="one_to_one",
        )
        .merge(paired_sera, on="serum", validate="many_to_one")
        .reset_index(drop=True)
    )

    drawn_strains = set(overlay["axis_label"])
    missing = set(chart_titers["axis_label"]) - drawn_strains
    if missing:
        raise ValueError(
            f"strains with no titer among the sera drawn: {sorted(missing)}; they would "
            "be an empty row against a tree tip"
        )

    keys = ["comparison", "subject", "axis_label"]
    sides = {
        condition: overlay[overlay["condition"] == condition]
        .drop(columns=["condition"])
        .rename(
            columns={
                "titer": f"{condition}_titer",
                "titer_bound": f"{condition}_titer_bound",
            }
        )
        for condition in CONDITIONS
    }
    # `serum` is kept from the post-vaccination side alone, so each subject is one line
    fold_change = sides["post"].merge(
        sides["pre"].drop(columns=["serum"]), on=keys, validate="one_to_one"
    )
    fold_change["fold_change"] = fold_change["post_titer"] / fold_change["pre_titer"]

    # a serum kept by the QC may still lack a titer against a given strain, leaving its
    # subject with only one side of that strain's comparison and so no fold change.
    # Counted per side, as a shortfall on one says nothing about the other.
    unmatched = {
        condition: len(sides[condition]) - len(fold_change) for condition in CONDITIONS
    }
    if any(unmatched.values()):
        print(
            f"  no fold change for {unmatched['pre']} subject-strain pair(s) measured "
            f"only before and {unmatched['post']} measured only after vaccination"
        )

    return overlay, fold_change


# ---- build every chart, slicing the data once per subtype and strain set -----------

for (subtype, strain_set), records in itertools.groupby(
    sorted(charts_to_make, key=lambda r: (r["subtype"], r["strain_set"])),
    key=lambda r: (r["subtype"], r["strain_set"]),
):
    chart_titers, chart_viruses, strain_order, _, _ = titer_charts.slice_chart_data(
        titers, viruses, subtype, strain_set, STRAIN_SETS[strain_set], subtype_params
    )
    overlay_data, fold_change_data = paired_chart_data(chart_titers, chart_viruses)
    print(
        f"{subtype} {strain_set}: {len(chart_viruses)} strains, "
        f"{len(overlay_data)} titers of {overlay_data['serum'].nunique()} sera, "
        f"{len(fold_change_data)} fold changes"
    )

    bases = {
        fold_change: titer_charts.base_chart(
            fold_change_data if fold_change else overlay_data,
            strain_order,
            plot_titer_summaries_params["facet_size"],
            [
                titer_charts.virus_selection,
                titer_charts.serum_selection,
                min_age_slider,
                max_age_slider,
            ],
        )
        for fold_change in [False, True]
    }

    for record in records:
        chart_type = CHART_TYPES[record["chart_type"]]
        fold_change = chart_type["fold_change"]
        layer = chart_type["build"](bases[fold_change])
        subtitle = pre_post_config["title"]
        if fold_change:
            # layered last so the thin rule draws over the interquartile band
            layer += no_change_line(fold_change_data)
            subtitle += "; dashed gray line marks no change in titer"
        chart = (
            layer.facet(row=alt.Row("comparison_n:N", title=None))
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
            # age is a property of the subject, so it takes both of its sera together and
            # cannot leave half a pair drawn
            .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
            .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
            # facet labels carry the number of subjects, which falls with the age sliders
            .transform_joinaggregate(
                n_subjects="distinct(subject)", groupby=["comparison"]
            )
            .transform_calculate(
                comparison_n="datum.comparison + ' (n=' + datum.n_subjects + ')'"
            )
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
        chart = titer_charts.finalize(chart, title, subtitle)

        print(f"Saving to {record['path']!r}")
        chart.save(record["path"])
