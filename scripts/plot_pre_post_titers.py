"""Interactive Altair plots comparing paired pre- and post-vaccination titers.

Writes a standalone HTML per subtype, strain set, tree coloring, and chart type. Two
comparisons are drawn: the pre- and post-vaccination titers overlaid in two colors, and
each subject's fold change from its pre- to its post-vaccination titer. Each facets by
the compared vaccination arm, and the recent-strain charts are drawn alongside their
subtype's tree as the summary charts are. A further chart draws the distribution of the
fold changes as box plots for the reference strain and a strain picked from a dropdown.

Only subjects whose pre- and post-vaccination sera both passed QC are drawn, so the
overlaid titers and the fold changes describe the same people. A comparison may restrict
itself to subjects recording a given metadata value, which draws one vaccination arm as
several facets.

"""

import itertools
import json
import sys

import altair as alt
import numpy as np
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
lower_titer_limit = plot_titer_summaries_params["draw_lower_titer_limit_line"]

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

# the facet label is rotated and so bounded by the facet height rather than the chart
# width, so a comparison's name may carry this character to say where it breaks onto a
# second line; the subject count follows the name, on whichever line ends it
FACET_LABEL_SPLIT = "|"

titers, metadata, sera_multicohort, viruses = titer_charts.load_and_validate(
    snakemake.input.titers_csv,
    snakemake.input.sera_csv,
    snakemake.input.sera_multicohort_csv,
    snakemake.input.viruses_csv,
    recent_vaccine_strains=recent_vaccine_strains,
    circulating_strain_type=circulating_strain_type,
    subtypes=subtypes,
    subtype_params=subtype_params,
    lower_titer_limit=lower_titer_limit,
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
    if comparison.count(FACET_LABEL_SPLIT) > 1:
        raise ValueError(
            f"comparison {comparison!r} carries more than one {FACET_LABEL_SPLIT!r}, "
            "which marks the one place its facet label breaks onto a second line"
        )
    subset_by = cohorts.get("subset_by")  # optional; absent draws the whole arm
    if set(cohorts) - {"subset_by"} != set(CONDITIONS):
        raise ValueError(
            f"comparison {comparison!r} must name exactly {CONDITIONS}, and may name "
            f"`subset_by`; got {sorted(cohorts)}"
        )
    if cohorts["pre"] == cohorts["post"]:
        raise ValueError(
            f"comparison {comparison!r} compares cohort {cohorts['pre']!r} with itself"
        )
    if subset_by is not None:
        if set(subset_by) != {"column", "value"}:
            raise ValueError(
                f"`subset_by` of comparison {comparison!r} must name exactly "
                f"['column', 'value'], got {sorted(subset_by)}"
            )
        if subset_by["column"] not in sera_multicohort.columns:
            raise ValueError(
                f"`subset_by` of comparison {comparison!r} names column "
                f"{subset_by['column']!r}, which is not in the sera metadata; columns "
                f"are {sorted(sera_multicohort.columns)}"
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

        if subset_by is not None:
            # a serum not recording the value subset on cannot be placed on either side
            # of the split, so it is a gap in the metadata rather than a serum to drop
            unrecorded = cohort_sera.loc[
                cohort_sera[subset_by["column"]].isnull(), "serum"
            ].tolist()
            if unrecorded:
                raise ValueError(
                    f"comparison {comparison!r} subsets cohort {cohort!r} on "
                    f"{subset_by['column']!r}, but these sera do not record it: "
                    f"{sorted(unrecorded)}"
                )
            cohort_sera = cohort_sera[
                cohort_sera[subset_by["column"]] == subset_by["value"]
            ]
            if cohort_sera.empty:
                raise ValueError(
                    f"comparison {comparison!r} draws no {cohort!r} serum with "
                    f"{subset_by['column']} == {subset_by['value']!r}"
                )

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

    if subset_by is not None:
        # a subject whose two sera disagree on the column subset on would have one serum
        # selected into this comparison and the other into its counterpart, leaving each
        # to drop it as unpaired; that is a gap in the sera metadata, and dropping it
        # would be indistinguishable from a serum lost to the QC
        arm_sera = sera_multicohort[
            sera_multicohort["cohort"].isin([cohorts["pre"], cohorts["post"]])
        ]
        per_subject = arm_sera.groupby(pair_by)[subset_by["column"]].nunique()
        inconsistent = per_subject[per_subject > 1]
        if len(inconsistent):
            disagreeing = arm_sera[arm_sera[pair_by].isin(inconsistent.index)]
            raise ValueError(
                f"comparison {comparison!r} subsets on {subset_by['column']!r}, but "
                f"these subjects do not record the same value for both of their sera:\n"
                + disagreeing.set_index(pair_by)[
                    ["serum", "cohort", subset_by["column"]]
                ].to_string()
            )

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

boxplot_comparison_haplotypes = pre_post_config["boxplot_comparison_haplotype"]
if set(boxplot_comparison_haplotypes) != set(subtypes):
    raise ValueError(
        f"`boxplot_comparison_haplotype` must name exactly the subtypes {sorted(subtypes)}"
        f", got {sorted(boxplot_comparison_haplotypes)}"
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
    baseline of. The aggregate is what makes that one rule per facet rather than one per
    titer: a `datum` encoding fixes where a mark sits, not how many are drawn.

    """
    return (
        alt.Chart(chart_data)
        .transform_aggregate(_n_rows="count()")
        .encode(y=alt.datum(1))
        .mark_rule(color="#888888", strokeWidth=1, strokeDash=[4, 3])
    )


def add_subject_metadata(chart):
    """Look up each serum's metadata and drop the subjects the age sliders exclude.

    Age is a property of the subject, so it takes both of that subject's sera together
    and cannot leave half a pair drawn.

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
        .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
        .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
    )


def label_facets_with_subject_count(chart):
    """Name each comparison with the number of subjects its facet draws.

    Counted over whatever rows have reached this point, so a chart that drops rows before
    calling this reports the count it actually draws rather than the whole arm's.

    """
    return chart.transform_joinaggregate(
        n_subjects="distinct(subject)", groupby=["comparison"]
    ).transform_calculate(
        comparison_n="datum.comparison + ' (n=' + datum.n_subjects + ')'"
    )


# the two boxes of a facet, in the order drawn and colored
VACCINE_ROLE = "recent vaccine strain"
COMPARISON_ROLE = "comparison strain"
BOXPLOT_COLORS = ["#1b5e8a", "#d97b29"]  # dark blue and orange
BOXPLOT_FACET_WIDTH = 170  # px
BOXPLOT_FACET_HEIGHT = 260  # px
BOXPLOT_BAND = BOXPLOT_FACET_WIDTH / 2  # px of x axis per box, as the bands abut
BOXPLOT_POINT_SPREAD = 22  # px each side of a box's center the subject points span

# spreads the subject points across their box; drawn once per point here rather than by a
# Vega expression, so that redrawing the chart does not move the points around
JITTER_SEED = 0
JITTER_SD = 4  # the spread covers this many standard deviations

BOXPLOT_TOOLTIPS = [
    alt.Tooltip("axis_label:N", title="strain"),
    alt.Tooltip("fold_change:Q", title="fold change", format=".3g"),
    *FOLD_CHANGE_SERUM_TOOLTIPS,
]


def fold_change_boxplot_chart(fold_change_data, chart_viruses, ref_axis_label, subtype):
    """Box plots of the fold changes for the reference strain and a selectable strain.

    Each facet draws one box for `ref_axis_label` and one for the strain picked in the
    dropdown, with a point per subject over each box. The dropdown starts on the
    subtype's `boxplot_comparison_haplotype`.

    A subject is drawn only where it has a fold change for both of the strains drawn, so
    the two boxes of a facet describe the same people and the facet's subject count is
    the count behind each of its boxes.

    """
    if ref_axis_label is None:
        raise ValueError(
            f"{subtype} has no `draw_titer_line` strain to compare against"
        )
    default = boxplot_comparison_haplotypes[subtype]
    options = sorted(set(chart_viruses["axis_label"]) - {ref_axis_label})
    if default not in options:
        raise ValueError(
            f"`boxplot_comparison_haplotype` {default!r} for {subtype} is not a recent "
            f"strain other than the reference {ref_axis_label!r}; options are {options}"
        )

    comparison_strain = alt.param(
        name="comparison_strain",
        value=default,
        bind=alt.binding_select(options=options, name="comparison strain: "),
    )
    # only the columns drawn, with fold changes rounded as the ratio is otherwise 18 digits
    frame = fold_change_data[
        [
            "comparison",
            "subject",
            "serum",
            "axis_label",
            "pre_titer",
            "pre_titer_bound",
            "post_titer",
            "post_titer_bound",
            "fold_change",
        ]
    ].copy()
    frame["fold_change"] = frame["fold_change"].round(3)
    frame["jitter"] = (
        np.random.default_rng(JITTER_SEED).standard_normal(len(frame)).round(3)
    )

    ref = json.dumps(ref_axis_label)
    base = alt.Chart(frame).add_params(
        comparison_strain, min_age_slider, max_age_slider
    )
    encodings = {
        "x": alt.X(
            "axis_label:N",
            sort=alt.EncodingSortField("role_order", op="min"),
            title=None,
            scale=alt.Scale(paddingInner=0, paddingOuter=0),
            axis=alt.Axis(labelAngle=0, labelLimit=200),
        ),
        "y": alt.Y(
            "fold_change:Q",
            title=VALUE_FOLD_CHANGE["title"],
            scale=titer_charts.titer_scale,
        ),
        "color": alt.Color(
            "role:N",
            title=None,
            scale=alt.Scale(
                domain=[VACCINE_ROLE, COMPARISON_ROLE], range=BOXPLOT_COLORS
            ),
            legend=alt.Legend(orient="left", columns=1),
        ),
    }
    box = base.mark_boxplot(extent="min-max", size=45, opacity=0.5).encode(**encodings)
    points = base.mark_circle(size=25, opacity=0.8).encode(
        **encodings,
        xOffset=alt.XOffset(
            "jitter:Q",
            title=None,
            scale=alt.Scale(
                domain=[-JITTER_SD, JITTER_SD],
                # an offset is measured from the start of the box's band
                range=[
                    BOXPLOT_BAND / 2 - BOXPLOT_POINT_SPREAD,
                    BOXPLOT_BAND / 2 + BOXPLOT_POINT_SPREAD,
                ],
            ),
        ),
        tooltip=BOXPLOT_TOOLTIPS,
    )
    layer = alt.layer(box, points, no_change_line(frame)).properties(
        width=BOXPLOT_FACET_WIDTH, height=BOXPLOT_FACET_HEIGHT
    )
    drawn = add_subject_metadata(
        layer.facet(column=alt.Column("comparison_n:N", title=None))
    ).transform_filter(
        f"datum.axis_label == {ref} || datum.axis_label == comparison_strain"
    )
    # the two boxes compare the strains within a subject, so they are only comparable
    # where they rest on the same subjects: one with a fold change for just one of the two
    # strains drawn is dropped from both boxes rather than from only the other
    paired = drawn.transform_joinaggregate(
        n_strains_drawn="distinct(axis_label)", groupby=["comparison", "subject"]
    ).transform_filter("datum.n_strains_drawn == 2")
    return label_facets_with_subject_count(paired).transform_calculate(
        role=f"datum.axis_label == {ref} ? {json.dumps(VACCINE_ROLE)} : "
        f"{json.dumps(COMPARISON_ROLE)}",
        role_order=f"datum.axis_label == {ref} ? 0 : 1",
    )


# each chart type names its title and how it is drawn. Most build a layer of marks onto a
# base faceted by vaccination arm and aligned to the tree, and so say whether they draw the
# fold changes or the titers themselves; a `build_chart` entry returns the whole chart
# instead, drawing two selected strains rather than a row per strain, and carries the
# drawing conventions its subtitle has to explain.
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
    "fold_change_boxplot": {
        "title": "box plots of the fold change from before to after vaccination",
        "build_chart": fold_change_boxplot_chart,
        "subtitle": (
            "the recent vaccine strain against a strain picked from the dropdown; only "
            "subjects with a fold change for both strains drawn are shown; boxes span "
            "the interquartile range and whiskers the full range; dashed gray line marks "
            "no change in titer"
        ),
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
    chart_titers, chart_viruses, strain_order, ref_axis_label, _ = (
        titer_charts.slice_chart_data(
            titers,
            viruses,
            subtype,
            strain_set,
            STRAIN_SETS[strain_set],
            subtype_params,
        )
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
        title = f"{chart_type['title']} for {subtype} {strain_set} strains"

        if "build_chart" in chart_type:
            chart = titer_charts.finalize(
                chart_type["build_chart"](
                    fold_change_data, chart_viruses, ref_axis_label, subtype
                ),
                title,
                f"{pre_post_config['title']}; {chart_type['subtitle']}",
                facet_label_split=FACET_LABEL_SPLIT,
            ).configure_headerColumn(labelOrient="top", labelAngle=0, labelPadding=4)
            print(f"Saving to {record['path']!r}")
            chart.save(record["path"])
            continue

        fold_change = chart_type["fold_change"]
        layer = chart_type["build"](bases[fold_change])
        subtitle = pre_post_config["title"]
        # layered last so the thin rule draws over the interquartile band
        if fold_change:
            layer += no_change_line(fold_change_data)
            subtitle += "; dashed gray line marks no change in titer"
        elif lower_titer_limit is not None:
            layer += titer_charts.lower_limit_line(overlay_data, lower_titer_limit)
            subtitle += "; " + titer_charts.LOWER_LIMIT_SUBTITLE.format(
                lower_titer_limit
            )
        chart = label_facets_with_subject_count(
            add_subject_metadata(
                layer.facet(row=alt.Row("comparison_n:N", title=None))
            ).transform_lookup(
                lookup="axis_label",
                from_=alt.LookupData(
                    data=chart_viruses,
                    key="axis_label",
                    fields=["virus", "strain_type", "subclade"],
                ),
            )
        )
        if strain_set == "recent":
            chart = titer_charts.add_tree(
                chart,
                tree_jsons[subtype],
                subtype_params[subtype],
                record["color_label"],
            )
            title += f", tree colored by {record['color_label']}"
        chart = titer_charts.finalize(
            chart, title, subtitle, facet_label_split=FACET_LABEL_SPLIT
        )

        print(f"Saving to {record['path']!r}")
        chart.save(record["path"])
