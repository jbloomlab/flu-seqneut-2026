"""Interactive Altair plots summarizing the titers of one serum group.

Writes a standalone HTML per subtype, strain set, tree coloring, and chart type. The
recent-strain charts are drawn alongside the
`nextstrain-prot-titers-tree` tree for their subtype, which supplies the strain order,
the strain labels, and the coloring. The vaccine-strain charts have no tree: they are
ordered by collection date and their labels name the vaccine type.

Chart and tree are joined on `derived_haplotype`.

"""

import itertools
import sys

import altair as alt
import pandas as pd
import tree_annotated_plot

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

_ = alt.data_transformers.disable_max_rows()

titers_csv = snakemake.input.titers_csv
sera_csv = snakemake.input.sera_csv
sera_multicohort_csv = snakemake.input.sera_multicohort_csv
viruses_csv = snakemake.input.viruses_csv
charts_to_make = snakemake.params.charts
tree_jsons = snakemake.params.tree_jsons
tree_params = snakemake.params.trees
recent_vaccine_strains = snakemake.params.recent_vaccine_strains
circulating_strain_type = snakemake.params.circulating_strain_type
plot_titer_summaries_params = snakemake.params.plot_titer_summaries_params
subtypes = snakemake.params.subtypes
facet_orientation = snakemake.params.facet_orientation
group = snakemake.wildcards.group

# pixels per strain along the strain axis
STRAIN_STEP = 11

# `strain_type` values plotted in each strain set. The recent vaccine strains are in both,
# labeled by haplotype alone among the recent strains and by vaccine type among the others.
STRAIN_SETS = {
    "recent": [circulating_strain_type, "recent_vaccine"],
    "vaccine": ["vaccine", "recent_vaccine"],
}

if facet_orientation not in {"vertical", "horizontal"}:
    raise ValueError(f"invalid {facet_orientation=}")

# Read titers
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

# Read sera metadata (one row per serum)
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

# Read sera multicohort (multiple rows per serum, one per cohort assignment)
sera_multicohort = pd.read_csv(sera_multicohort_csv, dtype={"serum": str})
print(f"Read {len(sera_multicohort)=} rows from {sera_multicohort_csv=}")

multicohort_sera = set(sera_multicohort["serum"])
if titers_sera != multicohort_sera:
    raise ValueError(
        f"Serum mismatch between titers and sera_multicohort.\n"
        f"  In titers but not multicohort: {titers_sera - multicohort_sera}\n"
        f"  In multicohort but not titers: {multicohort_sera - titers_sera}"
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

# Read viruses
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

missing_trees = set(subtypes) - set(tree_jsons)
if missing_trees:
    raise ValueError(f"no tree configured for subtype(s): {missing_trees}")

# ---- selections and sliders, shared by every chart --------------------------------

virus_selection = alt.selection_point(
    fields=["axis_label"], on="mouseover", empty=False, clear="mouseout", nearest=False
)

serum_selection = alt.selection_point(
    fields=["serum"], on="mouseover", empty=False, clear="mouseout", nearest=False
)

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

max_age = 5 * int(metadata["age_numeric"].max() // 5) + 5
assert all(metadata["age_numeric"] <= max_age)
min_age_slider = alt.param(
    value=0,
    bind=alt.binding_range(
        min=0, max=max_age, step=5, name="minimum subject age (years)"
    ),
)
max_age_slider = alt.param(
    value=max_age,
    bind=alt.binding_range(
        min=0, max=max_age, step=5, name="maximum subject age (years)"
    ),
)

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

# no explicit domain: `altair` fits it to whatever each chart draws, which is the spread of
# the individual titers for the per-serum charts but only the interquartile band for the
# interquartile-range ones. Leaving it to `altair` also means the axis follows the cohort
# and age filters, so a filtered view can never fall outside the axis.
titer_scale = alt.Scale(type="log", nice=False, padding=4)

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

# the strain axis is on y when facets run in columns, and on x when they run in rows
if facet_orientation == "vertical":
    StrainChannel, TiterChannel, FacetChannel = alt.Y, alt.X, alt.Column
    facet_channel = "column"
    panel_size = {
        "height": alt.Step(STRAIN_STEP),
        "width": plot_titer_summaries_params["facet_size"],
    }
else:
    StrainChannel, TiterChannel, FacetChannel = alt.X, alt.Y, alt.Row
    facet_channel = "row"
    panel_size = {
        "width": alt.Step(STRAIN_STEP),
        "height": plot_titer_summaries_params["facet_size"],
    }
strain_channel = "y" if facet_orientation == "vertical" else "x"
titer_channel = "x" if facet_orientation == "vertical" else "y"


def slice_chart_data(subtype, strain_set):
    """Return the titers, the virus lookup table, and the axis order for one strain set.

    Both frames are cut down to exactly the strains the chart draws. `tree_annotated_plot`
    compares the chart's strains against the tree's tips, and where no explicit axis sort
    is present it gathers them by walking every inline frame in the spec -- the lookup
    table included -- so an unsliced frame would contribute strains the chart never draws.
    Slicing also keeps the embedded data, and so the saved HTML, as small as the chart.

    """
    chart_viruses = viruses[
        (viruses["subtype"] == subtype)
        & viruses["strain_type"].isin(STRAIN_SETS[strain_set])
    ].copy()
    if not len(chart_viruses):
        raise ValueError(f"no {strain_set} strains for {subtype=}")

    if strain_set == "recent":
        # matches the tree tips, which are labeled by haplotype
        chart_viruses["axis_label"] = chart_viruses["derived_haplotype"]
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

    return (
        chart_titers,
        chart_viruses[["axis_label", "virus", "strain_type", "subclade"]],
        strain_order.tolist(),
    )


def base_chart(chart_titers, strain_order):
    """Return the chart every chart type is built from, sized one step per strain."""
    return (
        alt.Chart(chart_titers)
        .add_params(
            virus_selection,
            serum_selection,
            cohort_selection,
            min_age_slider,
            max_age_slider,
        )
        .encode(
            **{
                strain_channel: StrainChannel(
                    "axis_label",
                    sort=strain_order,
                    title=None,
                    axis=alt.Axis(labelLimit=500),
                )
            }
        )
        .properties(**panel_size)
    )


# fields carried per strain rather than repeated on every titer row
virus_tooltips = [
    alt.Tooltip("axis_label:N", title="strain"),
    alt.Tooltip("virus:N"),
    alt.Tooltip("strain_type:N"),
    alt.Tooltip("subclade:N"),
]

# groupby for the aggregates, which must carry through the looked-up virus annotations
virus_groupby = ["axis_label", "virus", "strain_type", "subclade", "cohort"]


def median_points(base):
    """Median titer per strain and cohort, as points."""
    return (
        base.transform_aggregate(median_titer="median(titer)", groupby=virus_groupby)
        .encode(
            **{
                titer_channel: TiterChannel(
                    "median_titer:Q", title="titer", scale=titer_scale
                )
            },
            tooltip=[*virus_tooltips, alt.Tooltip("median_titer:Q", format=".1f")],
            color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
            size=alt.condition(virus_selection, alt.value(80), alt.value(40)),
        )
        .mark_circle(opacity=1)
    )


def serum_lines(base):
    """One line per serum across the strains."""
    return base.encode(
        **{titer_channel: TiterChannel("titer", scale=titer_scale)},
        detail=alt.Detail("serum"),
        tooltip=[
            *virus_tooltips,
            alt.Tooltip("serum:N"),
            alt.Tooltip("titer", format=".1f"),
            alt.Tooltip("serum_collection_date:N", title="serum date"),
            alt.Tooltip("age:N", title="age"),
            alt.Tooltip("sex:N"),
        ],
        size=alt.condition(serum_selection, alt.value(3), alt.value(1.5)),
        opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.2)),
    ).mark_line()


def interquartile_range(base):
    """Shaded interquartile range of the titers for each strain."""
    return (
        base.transform_joinaggregate(
            median_titer="median(titer)",
            titer_q1="q1(titer)",
            titer_q3="q3(titer)",
            groupby=["axis_label"],
        )
        .encode(
            **{titer_channel: TiterChannel("titer", scale=titer_scale)},
            tooltip=[
                *virus_tooltips,
                alt.Tooltip("median_titer:Q", format=".1f"),
                alt.Tooltip("titer_q1:Q", format=".1f"),
                alt.Tooltip("titer_q3:Q", format=".1f"),
            ],
        )
        .mark_errorband(extent="iqr", opacity=0.5, interpolate="linear")
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
            **{
                titer_channel: TiterChannel(
                    "frac_below_cutoff:Q", title="fraction below cutoff"
                )
            },
            tooltip=[
                *virus_tooltips,
                alt.Tooltip("frac_below_cutoff:Q", format=".2f"),
            ],
            color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
        )
        .mark_bar(opacity=0.8)
    )


CHART_TYPES = {
    "individual_sera": (
        "median (points) and per-serum (lines) titers",
        lambda base: serum_lines(base) + median_points(base),
    ),
    "interquartile_range": (
        "median (points) and interquartile range titers",
        lambda base: interquartile_range(base) + median_points(base),
    ),
    "frac_below_cutoff": (
        "fraction sera below titer cutoff",
        frac_below_cutoff,
    ),
}


def facet_and_add_lookups(chart, chart_viruses):
    """Facet `chart` by cohort and look up the serum and virus annotations.

    Scoping when layering and faceting charts with `transform_lookup` requires the
    faceting to be done before the lookups, so both happen here.

    """
    return (
        chart.facet({facet_channel: FacetChannel("cohort_n:N", title=None)})
        .transform_lookup(
            lookup="serum",
            from_=alt.LookupData(
                data=metadata,
                key="serum",
                fields=[
                    "cohorts",
                    "serum_collection_date",
                    "age",
                    "age_numeric",
                    "sex",
                ],
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
        # flatten cohorts list (from sera_multicohort) to one row per cohort
        .transform_flatten(["cohorts"], as_=["cohort"])
        # filter by cohort and age
        .transform_filter(cohort_selection)
        .transform_filter(alt.datum["age_numeric"] >= min_age_slider)
        .transform_filter(alt.datum["age_numeric"] <= max_age_slider)
        # make facet labels w n per cohort
        .transform_joinaggregate(n_per_cohort="distinct(serum)", groupby=["cohort"])
        .transform_calculate(
            cohort_n="datum.cohort + ' (n=' + datum.n_per_cohort + ')'"
        )
    )


def add_tree(chart, subtype, color_label):
    """Return `chart` with its subtype's tree drawn alongside the strain axis.

    Tree tips absent from the chart are pruned, since the tree is built from the whole
    viral library while the titer QC may have dropped a strain. A chart strain absent
    from the tree is left to raise.

    """
    params = tree_params[subtype]
    return tree_annotated_plot.plot(
        tree_jsons[subtype],
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


def finalize(chart, title):
    """Stack the cohort legend under `chart` and apply the shared styling."""
    return (
        alt.vconcat(chart, cohort_legend, spacing=1)
        .resolve_scale(fill="independent", color="independent")
        .configure_axis(
            grid=False,
            titleFontWeight="normal",
            titleFontSize=13,
            labelOverlap=True,
        )
        .configure_header(
            title=None,
            labelOrient="top" if facet_orientation == "vertical" else "right",
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
            title=alt.TitleParams(title, anchor="middle", fontSize=13),
            # recompute the layout on every view update, not just at load: toggling a
            # cohort adds a facet, and without this the view keeps its original size and
            # clips the tree. `pad` is required -- the `fit*` types do not support
            # concatenated or faceted views.
            autosize=alt.AutoSizeParams(type="pad", resize=True),
        )
    )


# ---- build every chart, slicing the data once per subtype and strain set -----------

made_charts = []
for (subtype, strain_set), records in itertools.groupby(
    sorted(charts_to_make, key=lambda r: (r["subtype"], r["strain_set"])),
    key=lambda r: (r["subtype"], r["strain_set"]),
):
    chart_titers, chart_viruses, strain_order = slice_chart_data(subtype, strain_set)
    print(
        f"{subtype} {strain_set}: {len(chart_viruses)} strains, {len(chart_titers)} titers"
    )
    base = base_chart(chart_titers, strain_order)

    for record in records:
        chart_title, build = CHART_TYPES[record["chart_type"]]
        chart = facet_and_add_lookups(build(base), chart_viruses)
        title = f"{chart_title} for {subtype} {strain_set} strains"
        if strain_set == "recent":
            chart = add_tree(chart, subtype, record["color_label"])
            title += f", tree colored by {record['color_label']}"
        chart = finalize(chart, title)

        print(f"Saving to {record['path']!r}")
        chart.save(record["path"])
        made_charts.append(record["path"])

if set(made_charts) != set(snakemake.output.chart_htmls):
    raise ValueError(
        "charts written do not match the rule's declared outputs:\n"
        f"  not written: {sorted(set(snakemake.output.chart_htmls) - set(made_charts))}\n"
        f"  not declared: {sorted(set(made_charts) - set(snakemake.output.chart_htmls))}"
    )
