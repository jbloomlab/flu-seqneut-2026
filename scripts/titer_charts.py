"""Shared machinery for the interactive Altair charts of the final titers.

Holds what those charts have in common: reading and validating the final titer data,
slicing the data to the strains one chart draws, the marks the chart types are built
from, aligning a tree to the strain axis, and the styling of a finished chart.

What each script keeps for itself is which sera a chart draws, what it facets by, and
the chart types it assembles from these marks.

"""

import altair as alt
import pandas as pd
import tree_annotated_plot

# pixels per strain along the strain axis
STRAIN_STEP = 11

# no explicit domain: `altair` fits it to whatever each chart draws, which is the spread of
# the individual values for the per-serum charts but only the interquartile band for the
# interquartile-range ones. Leaving it to `altair` also means the axis follows the chart's
# filters, so a filtered view can never fall outside the axis. Fold changes span orders of
# magnitude for the same reason titers do, so both use this scale.
titer_scale = alt.Scale(type="log", nice=False, padding=4)

# what a chart plots on the titer axis: the field, its axis title, and its tooltip
# format. A list of strings is an axis title split over that many lines.
VALUE_TITER = {"field": "titer", "title": "titer", "format": ".1f"}

# fields carried per strain rather than repeated on every titer row
virus_tooltips = [
    alt.Tooltip("axis_label:N", title="strain"),
    alt.Tooltip("virus:N"),
    alt.Tooltip("strain_type:N"),
    alt.Tooltip("subclade:N"),
]

# groupby for the aggregates, which must carry through the looked-up virus annotations. A
# chart appends whatever else it splits the aggregate by, such as its facet field.
VIRUS_GROUPBY = ["axis_label", "virus", "strain_type", "subclade"]

virus_selection = alt.selection_point(
    fields=["axis_label"], on="mouseover", empty=False, clear="mouseout", nearest=False
)

serum_selection = alt.selection_point(
    fields=["serum"], on="mouseover", empty=False, clear="mouseout", nearest=False
)


def strain_sets(circulating_strain_type):
    """Return the `strain_type` values plotted in each strain set.

    The recent vaccine strains are in both, labeled by haplotype alone among the recent
    strains and by vaccine type among the others.

    """
    return {
        "recent": [circulating_strain_type, "recent_vaccine"],
        "vaccine": ["vaccine", "recent_vaccine"],
    }


def load_and_validate(
    titers_csv,
    sera_csv,
    sera_multicohort_csv,
    viruses_csv,
    *,
    recent_vaccine_strains,
    circulating_strain_type,
    subtypes,
    subtype_params,
):
    """Read the final titer data, validating it and the configuration against it.

    Returns the titers, the sera metadata (one row per serum), the multicohort sera
    assignments (one row per serum and cohort), and the viruses, with the recent vaccine
    strains marked `recent_vaccine` in `strain_type`.

    """
    titers = pd.read_csv(titers_csv, dtype={"serum": str, "virus": str})
    print(f"Read {len(titers)=} titers from {titers_csv=}")

    required_titer_cols = {"serum", "virus", "titer"}
    missing_titer_cols = required_titer_cols - set(titers.columns)
    if missing_titer_cols:
        raise ValueError(f"titers_csv missing required columns: {missing_titer_cols}")

    duplicate_pairs = titers.groupby(["serum", "virus"]).size()
    duplicate_pairs = duplicate_pairs[duplicate_pairs > 1]
    if len(duplicate_pairs) > 0:
        raise ValueError(
            f"Found {len(duplicate_pairs)} duplicate serum-virus pairs in titers:\n"
            f"{duplicate_pairs.head(10).to_string()}"
        )

    metadata = pd.read_csv(sera_csv, dtype={"serum": str})
    print(f"Read {len(metadata)=} sera from {sera_csv=}")

    required_sera_cols = {
        "serum",
        "cohort",
        "age_numeric",
        "serum_collection_date",
        "age",
        "sex",
    }
    missing_sera_cols = required_sera_cols - set(metadata.columns)
    if missing_sera_cols:
        raise ValueError(f"sera_csv missing required columns: {missing_sera_cols}")

    titers_sera = set(titers["serum"])
    metadata_sera = set(metadata["serum"])
    if titers_sera != metadata_sera:
        raise ValueError(
            f"Serum mismatch between titers and sera metadata.\n"
            f"  In titers but not metadata: {titers_sera - metadata_sera}\n"
            f"  In metadata but not titers: {metadata_sera - titers_sera}"
        )

    sera_multicohort = pd.read_csv(sera_multicohort_csv, dtype={"serum": str})
    print(f"Read {len(sera_multicohort)=} rows from {sera_multicohort_csv=}")

    multicohort_sera = set(sera_multicohort["serum"])
    if titers_sera != multicohort_sera:
        raise ValueError(
            f"Serum mismatch between titers and sera_multicohort.\n"
            f"  In titers but not multicohort: {titers_sera - multicohort_sera}\n"
            f"  In multicohort but not titers: {multicohort_sera - titers_sera}"
        )

    viruses = pd.read_csv(viruses_csv, dtype={"virus": str})
    print(f"Read {len(viruses)=} viruses from {viruses_csv=}")

    required_virus_cols = {
        "virus",
        "subtype",
        "strain_type",
        "subclade",
        "derived_haplotype",
        "vaccine_type",
        "virus_collection_date",
    }
    missing_virus_cols = required_virus_cols - set(viruses.columns)
    if missing_virus_cols:
        raise ValueError(f"viruses_csv missing required columns: {missing_virus_cols}")

    if not viruses["virus"].is_unique:
        raise ValueError(
            "viruses_csv has more than one row per virus: "
            f"{viruses.loc[viruses['virus'].duplicated(), 'virus'].tolist()}"
        )

    # `derived_haplotype` is the key joining each chart to its tree and is the axis label, so
    # a strain sharing one with another would silently collapse two strains onto one row
    if viruses["derived_haplotype"].isnull().any():
        raise ValueError(
            "viruses_csv has a null derived_haplotype for: "
            f"{viruses.loc[viruses['derived_haplotype'].isnull(), 'virus'].tolist()}"
        )
    shared_haplotypes = viruses[viruses.duplicated("derived_haplotype", keep=False)]
    if len(shared_haplotypes):
        raise ValueError(
            "derived_haplotype must be unique per strain, but these strains share one:\n"
            f"{shared_haplotypes[['virus', 'derived_haplotype']].sort_values('derived_haplotype')}"
        )

    missing_vaccine_strains = set(recent_vaccine_strains) - set(viruses["virus"])
    if missing_vaccine_strains:
        raise ValueError(
            f"recent_vaccine_strains not in viruses: {missing_vaccine_strains}"
        )

    valid_strain_types = {circulating_strain_type, "vaccine"}
    invalid_strain_types = set(viruses["strain_type"]) - valid_strain_types
    if invalid_strain_types:
        raise ValueError(
            f"Invalid strain_type values: {invalid_strain_types}. Expected: {valid_strain_types}"
        )

    # Mark recent vaccine strains as "recent_vaccine" strain_type
    viruses["strain_type"] = viruses["strain_type"].where(
        ~viruses["virus"].isin(recent_vaccine_strains), "recent_vaccine"
    )

    data_subtypes = set(viruses["subtype"].unique())
    if not set(subtypes).issubset(data_subtypes):
        raise ValueError(
            f"subtypes param {set(subtypes)} not all in viruses data. "
            f"Available subtypes: {data_subtypes}"
        )

    # a `draw_titer_line` haplotype naming no strain would silently draw no line
    for subtype in subtypes:
        haplotype = subtype_params[subtype]["draw_titer_line"]
        if haplotype is not None and haplotype not in set(
            viruses.loc[viruses["subtype"] == subtype, "derived_haplotype"]
        ):
            raise ValueError(
                f"`draw_titer_line` {haplotype!r} is not a derived_haplotype of any "
                f"{subtype} strain; use null to draw no line"
            )

    return titers, metadata, sera_multicohort, viruses


def age_sliders(metadata):
    """Return the sliders bounding which sera are drawn, by the subject's age."""
    max_age = 5 * int(metadata["age_numeric"].max() // 5) + 5
    assert all(metadata["age_numeric"] <= max_age)
    return tuple(
        alt.param(
            value=value,
            bind=alt.binding_range(min=0, max=max_age, step=5, name=name),
        )
        for value, name in [
            (0, "minimum subject age (years)"),
            (max_age, "maximum subject age (years)"),
        ]
    )


def slice_chart_data(
    titers, viruses, subtype, strain_set, strain_types, subtype_params
):
    """Return the titers, lookup tables, axis order, and reference strain for a set.

    Both frames are cut down to exactly the strains the chart draws. `tree_annotated_plot`
    compares the chart's strains against the tree's tips, and where no explicit axis sort
    is present it gathers them by walking every inline frame in the spec -- the lookup
    table included -- so an unsliced frame would contribute strains the chart never draws.
    Slicing also keeps the embedded data, and so the saved HTML, as small as the chart.

    """
    chart_viruses = viruses[
        (viruses["subtype"] == subtype) & viruses["strain_type"].isin(strain_types)
    ].copy()
    if not len(chart_viruses):
        raise ValueError(f"no {strain_set} strains for {subtype=}")

    if strain_set == "recent":
        # matches the tree tips, which are labeled by haplotype
        chart_viruses["axis_label"] = chart_viruses["derived_haplotype"]
        # `tree_annotated_plot` replaces this sort with the tree's tip order, so it serves
        # only to declare the strains the chart draws, as described above
        strain_order = chart_viruses.sort_values(["subclade", "virus"])["axis_label"]
    else:
        if chart_viruses["vaccine_type"].isnull().any():
            raise ValueError(
                "vaccine strains with no vaccine_type: "
                f"{chart_viruses.loc[chart_viruses['vaccine_type'].isnull(), 'virus'].tolist()}"
            )
        chart_viruses["axis_label"] = (
            chart_viruses["derived_haplotype"]
            + " ("
            + chart_viruses["vaccine_type"]
            + " vaccine)"
        )
        strain_order = chart_viruses.sort_values("virus_collection_date")["axis_label"]

    # `axis_label` is the axis field and the `transform_lookup` key, so it must be 1:1
    if not chart_viruses["axis_label"].is_unique:
        raise ValueError(
            f"axis labels are not unique for {subtype} {strain_set} strains:\n"
            f"{chart_viruses.loc[chart_viruses['axis_label'].duplicated(keep=False), ['virus', 'axis_label']]}"
        )

    # the vaccine set labels a strain `"D.3.1 (cell vaccine)"`, not by haplotype alone,
    # so the reference line's strain is resolved per strain set
    ref_haplotype = subtype_params[subtype]["draw_titer_line"]
    ref_axis_label = None
    if ref_haplotype is not None:
        ref_rows = chart_viruses[chart_viruses["derived_haplotype"] == ref_haplotype]
        if len(ref_rows) != 1:
            raise ValueError(
                f"`draw_titer_line` {ref_haplotype!r} matches {len(ref_rows)} of the "
                f"{subtype} {strain_set} strains, but must match exactly one"
            )
        ref_axis_label = ref_rows["axis_label"].item()

    chart_titers = titers.merge(
        chart_viruses[["virus", "axis_label"]],
        on="virus",
        how="inner",
        validate="many_to_one",
    )[["serum", "axis_label", "titer"]]

    # a strain with no titers would otherwise become an empty row on the axis
    strains_without_titers = set(chart_viruses["axis_label"]) - set(
        chart_titers["axis_label"]
    )
    if strains_without_titers:
        raise ValueError(
            f"{subtype} {strain_set} strains with no titers: {sorted(strains_without_titers)}"
        )

    # the median titer of each serum over the strains this chart draws, which the fold
    # changes are relative to and the median-titer sliders filter on. Kept out of
    # `chart_titers` and looked up instead, as it is one value per serum rather than per
    # titer, and rounded because a median over an even number of titers is otherwise 18
    # digits of precision far below a pixel.
    serum_medians = (
        chart_titers.groupby("serum")["titer"]
        .median()
        .round(1)
        .rename("median_titer_serum")
        .reset_index()
    )
    assert serum_medians["serum"].is_unique  # the key of a `transform_lookup`

    return (
        chart_titers,
        chart_viruses[["axis_label", "virus", "strain_type", "subclade"]],
        strain_order.tolist(),
        ref_axis_label,
        serum_medians,
    )


def base_chart(chart_titers, strain_order, facet_size, params):
    """Return the chart every chart type is built from, sized one step per strain.

    `params` are declared here because a param may be declared only once per spec, so
    every mark layered onto this chart shares the ones it needs.

    """
    return (
        alt.Chart(chart_titers)
        .add_params(*params)
        .encode(
            x=alt.X(
                "axis_label",
                sort=strain_order,
                title=None,
                axis=alt.Axis(labelLimit=500),
            )
        )
        .properties(width=alt.Step(STRAIN_STEP), height=facet_size)
    )


def median_points(
    base, value, *, groupby, aggregate_extras=None, tooltip_extras=(), color=None
):
    """Median of the plotted value per strain, as points.

    `aggregate_extras` names further aggregates to compute, and `tooltip_extras` what to
    show alongside the median, so a chart plotting something other than the titer can
    still tooltip the titers behind it. Where `color` is given it replaces the hover
    highlight on that channel, which then shows only in the point size.

    """
    field = value["field"]
    aggregates = {f"median_{field}": f"median({field})", **(aggregate_extras or {})}
    return (
        base.transform_aggregate(**aggregates, groupby=groupby)
        .encode(
            y=alt.Y(f"median_{field}:Q", title=value["title"], scale=titer_scale),
            tooltip=[
                *virus_tooltips,
                alt.Tooltip(f"median_{field}:Q", format=value["format"]),
                *tooltip_extras,
            ],
            color=(
                color
                if color is not None
                else alt.condition(
                    virus_selection, alt.value("red"), alt.value("black")
                )
            ),
            size=alt.condition(virus_selection, alt.value(92), alt.value(46)),
        )
        .mark_circle(opacity=1)
    )


def serum_lines(base, value, *, tooltip_extras=(), color=None):
    """One line per serum across the strains, thickened where hovered."""
    field = value["field"]
    return base.encode(
        y=alt.Y(f"{field}:Q", title=value["title"], scale=titer_scale),
        detail=alt.Detail("serum"),
        tooltip=[
            *virus_tooltips,
            alt.Tooltip("serum:N"),
            alt.Tooltip(f"{field}:Q", format=value["format"]),
            *tooltip_extras,
        ],
        size=alt.condition(serum_selection, alt.value(3), alt.value(1.5)),
        opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.2)),
        **({} if color is None else {"color": color}),
    ).mark_line()


def interquartile_range(
    base, value, *, groupby, aggregate_extras=None, tooltip_extras=(), color=None
):
    """Shaded interquartile range of the plotted value for each strain."""
    field = value["field"]
    aggregates = {
        f"median_{field}": f"median({field})",
        f"{field}_q1": f"q1({field})",
        f"{field}_q3": f"q3({field})",
    }
    tooltips = [alt.Tooltip(f"{name}:Q", format=value["format"]) for name in aggregates]
    aggregates.update(aggregate_extras or {})
    return (
        base.transform_joinaggregate(**aggregates, groupby=groupby)
        .encode(
            y=alt.Y(f"{field}:Q", title=value["title"], scale=titer_scale),
            tooltip=[*virus_tooltips, *tooltips, *tooltip_extras],
            **({} if color is None else {"color": color}),
        )
        .mark_errorband(extent="iqr", opacity=0.5, interpolate="linear")
    )


def add_tree(chart, tree_json, params, color_label):
    """Return `chart` with its subtype's tree drawn alongside the strain axis.

    Tree tips absent from the chart are pruned, since the tree is built from the whole
    viral library while the titer QC may have dropped a strain. A chart strain absent
    from the tree is left to raise.

    """
    return tree_annotated_plot.plot(
        tree_json,
        chart,
        chart_strain_field="axis_label",
        tree_strain_field="derived_haplotype",
        branch_length="div",
        tree_size=params["tree_size"],
        color_tree_by=params["color_trees_by"][color_label],
        prune_tree_to_chart=True,
        prune_chart_to_tree=False,
        connect_leader_to_label=True,
        strain_label_font_size=params["strain_label_font_size"],
        shift_tree_loc=params["shift_tree_loc"],
        tree_color_legend_format={"orient": "left", "columns": 1},
    )


def finalize(chart, title, subtitle, *, above=(), below=()):
    """Stack `above` and `below` around `chart`, and style it.

    Must run after `add_tree`, which hoists the top-level attributes off the chart it is
    given so it can nest it, and so would drop anything set here.

    """
    return (
        alt.vconcat(*above, chart, *below, spacing=1, center=True)
        .resolve_scale(fill="independent", color="independent")
        .configure_axis(
            grid=False,
            titleFontWeight="normal",
            titleFontSize=13,
            labelOverlap=True,
        )
        .configure_header(
            title=None,
            labelOrient="right",
            labelFontSize=13,
            labelPadding=2,
        )
        .configure_view(stroke="black")
        .configure_facet(spacing=8)
        .configure_legend(
            labelFontSize=12,
            titleFontSize=13,
            titleLimit=0,  # no truncation; the cohort legend's title is wider than the 180px default
            symbolStrokeWidth=1,
            symbolOpacity=1,
            symbolStrokeColor="black",
            columns=12,
            orient="bottom",
        )
        .properties(
            title=alt.TitleParams(
                title, subtitle=subtitle, anchor="middle", fontSize=13
            ),
            # recompute the layout on every view update, not just at load: toggling a
            # cohort adds a facet, and without this the view keeps its original size and
            # clips the tree. `pad` is required -- the `fit*` types do not support
            # concatenated or faceted views.
            autosize=alt.AutoSizeParams(type="pad", resize=True),
        )
    )
