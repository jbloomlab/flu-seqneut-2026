"""Interactive Altair plots summarizing the titers of one serum group.

Writes a standalone HTML per subtype, strain type, and chart type, plus a report
embedding the first of them.

"""

import itertools
import json
import os
import sys

import altair as alt
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# `snakemake` puts only this script's own directory on the path, so the report module in
# the pipeline submodule has to be added before it can be imported
sys.path.append(os.path.dirname(snakemake.input.report_module))

from seqneut_report import Report

_ = alt.data_transformers.disable_max_rows()

titers_csv = snakemake.input.titers_csv
sera_csv = snakemake.input.sera_csv
sera_multicohort_csv = snakemake.input.sera_multicohort_csv
viruses_csv = snakemake.input.viruses_csv
chart_htmls = snakemake.output.chart_htmls
recent_vaccine_strains = snakemake.params.recent_vaccine_strains
circulating_strain_type = snakemake.params.circulating_strain_type
plot_titer_summaries_params = snakemake.params.plot_titer_summaries_params
subtypes = snakemake.params.subtypes
facet_orientation = snakemake.params.facet_orientation
group = snakemake.wildcards.group

report = Report(title=f"Titer summary plots for {group} ({facet_orientation} facets)")

report.md("""
    Interactive Altair plots showing median titers, individual serum titers,
    interquartile ranges, and fraction of sera below titer cutoffs.

    ## Read data
    """)

# Read titers
titers = pd.read_csv(titers_csv, dtype={"serum": str, "virus": str})
report.md(f"Read {len(titers)=} titers from {titers_csv=}", log=True)

# Validate required columns in titers
required_titer_cols = {"serum", "virus", "titer"}
missing_titer_cols = required_titer_cols - set(titers.columns)
if missing_titer_cols:
    raise ValueError(f"titers_csv missing required columns: {missing_titer_cols}")

# Read sera metadata (one row per serum)
metadata = pd.read_csv(sera_csv, dtype={"serum": str})
report.md(f"Read {len(metadata)=} sera from {sera_csv=}", log=True)

# Validate required columns in sera metadata
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

# Validate sera match between titers and metadata
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
report.md(f"Read {len(sera_multicohort)=} rows from {sera_multicohort_csv=}", log=True)

# Validate sera match between titers and multicohort
multicohort_sera = set(sera_multicohort["serum"])
if titers_sera != multicohort_sera:
    raise ValueError(
        f"Serum mismatch between titers and sera_multicohort.\n"
        f"  In titers but not multicohort: {titers_sera - multicohort_sera}\n"
        f"  In multicohort but not titers: {multicohort_sera - titers_sera}"
    )

# Validate "All" cohort exists in multicohort
if "All" not in sera_multicohort["cohort"].values:
    raise ValueError(
        "Expected 'All' cohort in sera_multicohort but not found. "
        f"Available cohorts: {sera_multicohort['cohort'].unique().tolist()}"
    )

# Get list of all cohorts (for legend), with "All" first
all_cohorts = sera_multicohort["cohort"].unique().tolist()
all_cohorts = ["All"] + sorted([c for c in all_cohorts if c != "All"])
report.md(f"Cohorts: {all_cohorts}", log=True)

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
        + (
            f" (and {len(sera_missing_cohorts) - 10} more)"
            if len(sera_missing_cohorts) > 10
            else ""
        )
    )

# Read viruses
viruses = pd.read_csv(viruses_csv, dtype={"virus": str})
report.md(f"Read {len(viruses)=} viruses from {viruses_csv=}", log=True)

# Validate required columns in viruses
required_virus_cols = {
    "virus",
    "subtype",
    "strain_type",
    "subclade",
    "derived_haplotype",
    "vaccine_type",
}
missing_virus_cols = required_virus_cols - set(viruses.columns)
if missing_virus_cols:
    raise ValueError(f"viruses_csv missing required columns: {missing_virus_cols}")

# Validate recent_vaccine_strains are in viruses
missing_vaccine_strains = set(recent_vaccine_strains) - set(viruses["virus"])
if missing_vaccine_strains:
    raise ValueError(
        f"recent_vaccine_strains not in viruses: {missing_vaccine_strains}"
    )

# Validate strain_type values
valid_strain_types = {circulating_strain_type, "vaccine"}
invalid_strain_types = set(viruses["strain_type"]) - valid_strain_types
if invalid_strain_types:
    raise ValueError(
        f"Invalid strain_type values: {invalid_strain_types}. "
        f"Expected: {valid_strain_types}"
    )

# Mark recent vaccine strains as "recent_vaccine" strain_type
viruses["strain_type"] = viruses["strain_type"].where(
    ~viruses["virus"].isin(recent_vaccine_strains), "recent_vaccine"
)

# Validate subtypes param matches virus data
data_subtypes = set(viruses["subtype"].unique())
param_subtypes = set(subtypes)
if not param_subtypes.issubset(data_subtypes):
    raise ValueError(
        f"subtypes param {param_subtypes} not all in viruses data. "
        f"Available subtypes: {data_subtypes}"
    )

# Compute strain plot order: sort by subclade then alphabetically
# Strains without subclade (vaccine strains) sort to end
viruses_sorted = viruses.sort_values(
    by=["subclade", "virus"],
    key=lambda col: col.fillna("zzz") if col.name == "subclade" else col,
)
viral_strain_plot_order = viruses_sorted["virus"].tolist()
assert set(viruses["virus"]) == set(viral_strain_plot_order)
report.md(
    f"Computed strain plot order: {len(viral_strain_plot_order)} viruses", log=True
)

report.md("""
    ## Plot all the titers

    ### Assign label colors by subclade / vaccine type
    Define color mapping from subclade (circulating strains) or vaccine type (vaccine
    strains) to colors for label coloring, then create an expression that can be passed
    to `altair` *labelColor*:
    """)

strain_color_prop = viruses.assign(
    strain=lambda x: pd.Categorical(x["virus"], viral_strain_plot_order, ordered=True),
    color_prop=lambda x: x["subclade"].where(
        x["strain_type"] == circulating_strain_type, x["vaccine_type"] + " vaccine"
    ),
).sort_values("strain")

assert strain_color_prop["color_prop"].notnull().all()
assert set(viruses["virus"]) == set(strain_color_prop["strain"])

viruses["color_prop"] = viruses["virus"].map(
    strain_color_prop.set_index("strain")["color_prop"].to_dict()
)

prop_colors = dict(plot_titer_summaries_params["prop_colors"])

other_prop_colors = plot_titer_summaries_params["other_prop_colors"]

for _subtype in strain_color_prop["subtype"].unique():
    subtype_color_props = (
        strain_color_prop[strain_color_prop["subtype"] == _subtype]["color_prop"]
        .unique()
        .tolist()
    )
    props_not_yet_colored = [p for p in subtype_color_props if p not in prop_colors]
    if len(props_not_yet_colored) > len(other_prop_colors):
        raise ValueError(
            f"props_not_yet_colored={props_not_yet_colored!r} longer than "
            f"other_prop_colors={other_prop_colors!r}"
        )
    prop_colors.update(dict(zip(props_not_yet_colored, other_prop_colors)))

report.table(pd.Series(prop_colors).rename("color").rename_axis("property").to_frame())
assert set(strain_color_prop["color_prop"]).issubset(prop_colors)

strain_color_prop = strain_color_prop.assign(
    color=lambda x: x["color_prop"].map(prop_colors)
)

color_mapping = strain_color_prop.set_index("virus")["color"].to_dict()

# make a different color map for each subtype as they are plotted separately
labelColor_expr = f"({json.dumps(color_mapping)})[datum.label] || 'black'"

# Create label mapping: derived_haplotype if exists, else strain (without subtype suffix)
label_mapping = {
    row["virus"]: (
        row["derived_haplotype"]
        if pd.notna(row["derived_haplotype"])
        else row["virus"].rsplit("_", 1)[0]  # remove _H1N1 or _H3N2 suffix
    )
    for _, row in viruses.iterrows()
}
labelText_expr = f"({json.dumps(label_mapping)})[datum.label] || datum.label"

report.md("""
    ### Now make nicely formatted charts
    """)

# First set up the base chart and selections

facet_size = plot_titer_summaries_params["facet_size"]
if facet_orientation not in {"vertical", "horizontal"}:
    raise ValueError(
        f"facet_orientation must be 'vertical' or 'horizontal', got {facet_orientation!r}"
    )
titer_encoding = "x" if facet_orientation == "vertical" else "y"

# Validate no duplicate serum-virus pairs in titers
duplicate_pairs = titers.groupby(["serum", "virus"]).size()
duplicate_pairs = duplicate_pairs[duplicate_pairs > 1]
if len(duplicate_pairs) > 0:
    raise ValueError(
        f"Found {len(duplicate_pairs)} duplicate serum-virus pairs in titers:\n"
        f"{duplicate_pairs.head(10).to_string()}"
    )

# Validate each virus has unique subtype and strain_type
virus_groups = viruses.groupby("virus")[["subtype", "strain_type"]].nunique()
inconsistent = virus_groups[(virus_groups > 1).any(axis=1)]
if len(inconsistent) > 0:
    raise ValueError(
        f"Found viruses with inconsistent subtype/strain_type:\n{inconsistent}"
    )

report.md(f"Plotting cohorts={all_cohorts} with {facet_orientation=}", log=True)

virus_selection = alt.selection_point(
    fields=["virus"], on="mouseover", empty=False, clear="mouseout", nearest=False
)

serum_selection = alt.selection_point(
    fields=["serum"], on="mouseover", empty=False, clear="mouseout", nearest=False
)

cohort_selection = alt.selection_point(
    fields=["cohort"], bind="legend", empty="all", toggle="true", clear=False
)

# select by color used to color strain labels
color_prop_selection = alt.selection_point(
    fields=["color_prop"], bind="legend", empty="all", toggle="true", clear=False
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

# make the chart base, using transform_lookup to make it as small as possible
# by looking up serum-specific and virus-specific annotations
titers_base_nolookup = (
    alt.Chart(titers[["serum", "virus", "titer"]])
    .add_params(
        virus_selection,
        serum_selection,
        cohort_selection,
        color_prop_selection,
        min_age_slider,
        max_age_slider,
    )
    .encode(
        **{
            ("y" if facet_orientation == "vertical" else "x"): alt.Y(
                "virus",
                sort=list(reversed(viral_strain_plot_order)),
                axis=alt.Axis(
                    labelLimit=500,
                    labelColor={"expr": labelColor_expr},
                    labelFontWeight=600,  # make a bit bolder so colors show
                    labelExpr=labelText_expr,
                ),
            ),
        },
    )
    .properties(
        **(
            {"height": alt.Step(11), "width": facet_size}
            if facet_orientation == "vertical"
            else {"width": alt.Step(11), "height": facet_size}
        )
    )
)

# dummy chart to bind the selectable legend for serum cohort
dummy_cohort_chart = (
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


def facet_and_add_lookups(chart):
    """Facet `chart` by cohort and look up the serum and virus annotations.

    Scoping when layering and faceting charts with `transform_lookup` requires the
    faceting to be done before the lookups, so both happen here.

    """
    return (
        chart
        # facet
        .facet(
            {
                ("column" if facet_orientation == "vertical" else "row"): alt.Column(
                    "cohort_n:N", title=None
                )
            }
        )
        # lookup additional data
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
            lookup="virus",
            from_=alt.LookupData(
                data=viruses,
                key="virus",
                fields=[
                    "subtype",
                    "strain_type",
                    "subclade",
                    "color_prop",
                    "derived_haplotype",
                ],
            ),
        )
        # flatten cohorts list (from sera_multicohort) to one row per cohort
        .transform_flatten(["cohorts"], as_=["cohort"])
        # filter by property used to color strain labels
        .transform_filter(color_prop_selection)
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


# set titer scale
titer_lower_limit = plot_titer_summaries_params["titer_lower_limit"]
report.md(f"Using {titer_lower_limit=}", log=True)
titer_scale = alt.Scale(type="log", nice=False, domainMin=titer_lower_limit, padding=4)

# make median titer point chart
median_points = (
    titers_base_nolookup.transform_aggregate(
        median_titer="median(titer)",
        groupby=[
            "virus",
            "subtype",
            "derived_haplotype",
            "strain_type",
            "subclade",
            "cohort",
        ],
    )
    .encode(
        **{titer_encoding: alt.X("median_titer:Q", title="titer", scale=titer_scale)},
        tooltip=[
            "virus",
            "derived_haplotype:N",
            alt.Tooltip("median_titer:Q", format=".1f"),
            "strain_type:N",
            "subclade:N",
        ],
        color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
        size=alt.condition(virus_selection, alt.value(80), alt.value(40)),
    )
    .mark_circle(opacity=1)
)

# make per-serum lines
serum_lines = titers_base_nolookup.encode(
    **{titer_encoding: alt.X("titer", scale=titer_scale)},
    detail=alt.Detail("serum"),
    tooltip=[
        "virus",
        "derived_haplotype:N",
        "serum",
        alt.Tooltip("titer", format=".1f"),
        alt.Tooltip("serum_collection_date:N", title="serum date"),
        alt.Tooltip("age:N", title="age"),
        "sex:N",
    ],
    size=alt.condition(serum_selection, alt.value(3), alt.value(1.5)),
    opacity=alt.condition(serum_selection, alt.value(1), alt.value(0.2)),
).mark_line()

# make interquartile range chart
interquartile_range = (
    titers_base_nolookup.transform_joinaggregate(
        median_titer="median(titer)",
        titer_q1="q1(titer)",
        titer_q3="q3(titer)",
        groupby=["virus"],
    )
    .encode(
        **{titer_encoding: alt.X("titer", scale=titer_scale)},
        tooltip=[
            "virus",
            "derived_haplotype:N",
            alt.Tooltip("median_titer:Q", format=".1f"),
            alt.Tooltip("titer_q1:Q", format=".1f"),
            alt.Tooltip("titer_q3:Q", format=".1f"),
            "strain_type:N",
            "subclade:N",
        ],
    )
    .mark_errorband(extent="iqr", opacity=0.5, interpolate="linear")
)

# make fraction below titer cutoff chart

titer_cutoff = plot_titer_summaries_params["titer_cutoff"]
report.md(f"Setting initial {titer_cutoff=}", log=True)

titer_cutoff_slider = alt.param(
    value=titer_cutoff,
    bind=alt.binding_range(
        min=titer_lower_limit,
        max=1000,
        step=5,
        name="fraction sera below this cutoff",
    ),
)

# make titer cutoff chart
frac_below_cutoff = (
    titers_base_nolookup.add_params(titer_cutoff_slider)
    .transform_calculate(below_cutoff=alt.datum["titer"] < titer_cutoff_slider)
    .transform_aggregate(
        n_below_cutoff="sum(below_cutoff)",
        n_total="distinct(serum)",
        groupby=[
            "virus",
            "subtype",
            "derived_haplotype",
            "strain_type",
            "subclade",
            "cohort",
        ],
    )
    .transform_calculate(
        frac_below_cutoff=alt.datum["n_below_cutoff"] / alt.datum["n_total"]
    )
    .encode(
        **{titer_encoding: alt.X("frac_below_cutoff:Q", title="fraction below cutoff")},
        tooltip=[
            "virus",
            "derived_haplotype:N",
            alt.Tooltip("frac_below_cutoff:Q", format=".2f"),
            "strain_type:N",
            "subclade:N",
        ],
        color=alt.condition(virus_selection, alt.value("red"), alt.value("black")),
    )
    .mark_bar(opacity=0.8)
)

made_chart = {c: False for c in chart_htmls}

for _subtype, strain_type, (chart_obj, chart_desc, title) in itertools.product(
    subtypes,
    ["recent", "vaccine"],
    [
        (
            serum_lines + median_points,
            "individual_sera",
            "median (points) and per-serum (lines) titers",
        ),
        (
            interquartile_range + median_points,
            "interquartile_range",
            "median (points) and interquartile range titers",
        ),
        (
            frac_below_cutoff,
            "frac_below_cutoff",
            "fraction sera below titer cutoff",
        ),
    ],
):
    filepattern = f"{_subtype}_{strain_type}_{chart_desc}"
    filename = [c for c in chart_htmls if filepattern in c]
    assert (
        len(filename) == 1
    ), f"did not find one filepattern={filepattern!r} in chart_htmls={chart_htmls!r}"
    filename = filename[0]

    # strain types to plot
    strain_types = {
        "recent": [circulating_strain_type, "recent_vaccine"],
        "vaccine": ["vaccine", "recent_vaccine"],
    }[strain_type]

    # ---- Make the legend for the colored strain labels ------------------------
    # get the virus colors plotted for the labels
    plotted_colors = strain_color_prop[
        (strain_color_prop["subtype"] == _subtype)
        & (strain_color_prop["strain_type"].isin(strain_types))
    ][["color_prop", "color"]].drop_duplicates()

    label_color_legend = (
        alt.Chart(plotted_colors)
        .add_params(color_prop_selection)
        .mark_point(opacity=0)  # invisible mark; we just want the legend
        .encode(
            fill=alt.Fill(
                "color_prop",
                title="virus type (click to select)",
                scale=alt.Scale(
                    domain=list(reversed(plotted_colors["color_prop"].tolist())),
                    range=list(reversed(plotted_colors["color"].tolist())),
                ),
                legend=alt.Legend(symbolType="square"),
            )
        )
        .properties(width=1, height=1)  # tiny plot; legend renders outside
    )
    # ---- Finished making the legend for the colored strain labels -------------

    chart = (
        alt.vconcat(
            (
                facet_and_add_lookups(chart_obj)
                .transform_filter(alt.datum["subtype"] == _subtype)
                .transform_filter(alt.FieldOneOfPredicate("strain_type", strain_types))
            ),
            label_color_legend,
            dummy_cohort_chart,
            spacing=1,
        )
        .resolve_scale(fill="independent")
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
            symbolStrokeWidth=1,
            symbolOpacity=1,
            symbolStrokeColor="black",
            columns=12,
            orient="bottom",
        )
        .properties(
            title=alt.TitleParams(
                f"{title} for {_subtype} {strain_type} strains",
                anchor="middle",
                fontSize=13,
            )
        )
    )
    if not any(made_chart.values()):
        report.md("Displaying just the first chart here (since they are large).")
        report.chart(chart)

    report.md(f"Saving to filename={filename!r}", log=True)
    chart.save(filename)

    made_chart[filename] = True

assert all(made_chart.values()), f"made_chart={made_chart!r}"

report.write(snakemake.output.html)
print(f"Wrote the report to {snakemake.output.html}")
