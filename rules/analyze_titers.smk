"""Aggregation of sera metadata, QC of the final titers, and plots of those titers."""

import os
import re

groups_to_analyze = config["groups_to_analyze"]
_unknown_groups = set(groups_to_analyze) - set(groups)
if _unknown_groups:
    raise ValueError(
        f"`groups_to_analyze` names groups with no plates: {sorted(_unknown_groups)}"
    )


rule aggregate_sera_metadata:
    """Aggregate and validate sera metadata from multiple cohorts."""
    input:
        csvs=config["sera_metadata"],
    output:
        csv="results/sera_metadata/all_sera_metadata.csv",
    log:
        "results/logs/aggregate_sera_metadata.txt",
    conda:
        "../seqneut-pipeline/environment.yml"
    script:
        "../scripts/aggregate_sera_metadata.py"


rule process_final_titer_data:
    """Process and QC final titer data for each group."""
    input:
        sera_metadata=rules.aggregate_sera_metadata.output.csv,
        viral_library=lambda wc: config["viral_libraries"][
            config["process_final_titer_data"]["viral_library"]
        ],
        titers="results/aggregated_titers/titers_{group}.csv",
    output:
        titers="results/final_titer_data/{group}_titers.csv",
        sera="results/final_titer_data/{group}_sera.csv",
        sera_multicohort="results/final_titer_data/{group}_sera_multicohort.csv",
        viruses="results/final_titer_data/{group}_viruses.csv",
        titers_summarized="results/final_titer_data/{group}_titers_summarized_by_virus.csv",
        summary="results/final_titer_data/{group}_summary.txt",
    log:
        "results/logs/process_final_titer_data_{group}.txt",
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        config=config["process_final_titer_data"],
    script:
        "../scripts/process_final_titer_data.py"


# Strain sets each chart type is made for. The fold-change charts are recent-strain
# only: they plot each titer relative to the serum's median over the strains the chart
# draws, which is not a baseline worth plotting against for the few vaccine strains.
titer_chart_types = {
    "individual_sera": ["recent", "vaccine"],
    "interquartile_range": ["recent", "vaccine"],
    "frac_below_cutoff": ["recent", "vaccine"],
    "individual_sera_fold_change": ["recent"],
    "interquartile_range_fold_change": ["recent"],
}

# Define set of titer charts to make, as (output path template, chart record) pairs
titer_charts = []
for _subtype in config["subtypes"]:
    _colorings = config["plot_titer_summaries_params"]["subtype_params"][_subtype][
        "color_trees_by"
    ]
    for _chart_type, _strain_sets in titer_chart_types.items():
        for _strain_set, _color_label in [
            *(("recent", _label) for _label in _colorings),
            ("vaccine", None),
        ]:
            if _strain_set not in _strain_sets:
                continue
            if _color_label and not re.fullmatch(r"[A-Za-z0-9_-]+", _color_label):
                raise ValueError(
                    f"`color_trees_by` key {_color_label!r} has invalid characters."
                )
            _stem = "_".join(
                part
                for part in (_subtype, _strain_set, _chart_type, _color_label)
                if part
            )
            titer_charts.append(
                (
                    f"results/titer_plots/{{group}}_{_stem}_{{orientation}}.html",
                    {
                        "subtype": _subtype,
                        "strain_set": _strain_set,
                        "color_label": _color_label,
                        "chart_type": _chart_type,
                    },
                )
            )


rule plot_titer_summaries:
    """Create interactive Altair titer summary plots for each serum group."""
    input:
        titers_csv=rules.process_final_titer_data.output.titers,
        sera_csv=rules.process_final_titer_data.output.sera,
        sera_multicohort_csv=rules.process_final_titer_data.output.sera_multicohort,
        viruses_csv=rules.process_final_titer_data.output.viruses,
        # one tree per subtype, in the same order as `params.subtypes`
        trees=[
            config["nextstrain-prot-titers-tree_config"][subtype]["auspice_json"]
            for subtype in config["subtypes"]
        ],
    output:
        chart_htmls=[template for template, _ in titer_charts],
    log:
        "results/logs/plot_titer_summaries_{group}_{orientation}.txt",
    wildcard_constraints:
        group="|".join(groups),
        orientation="vertical|horizontal",
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        # in the same order as `output.chart_htmls`, which names each chart's path
        charts=[record for _, record in titer_charts],
        recent_vaccine_strains=config["recent_vaccine_strains"],
        circulating_strain_type=config["circulating_strain_type"],
        plot_titer_summaries_params=config["plot_titer_summaries_params"],
        subtypes=config["subtypes"],
        facet_orientation=lambda wc: wc.orientation,
    script:
        "../scripts/plot_titer_summaries.py"


# Add titer summary plots to docs HTMLs generated by pipeline, one section per tree
# coloring plus one for the tree-less vaccine-strain charts
for _group in groups_to_analyze:
    for _orientation in ["vertical", "horizontal"]:
        for _template, _chart in titer_charts:
            if _chart["strain_set"] == "recent":
                _section = (
                    f"Interactive charts of {_group} titers, "
                    f"tree colored by {_chart['color_label']}"
                )
                _subsection = f"{_chart['subtype']} recent strains ({_orientation})"
            else:
                _section = f"Interactive charts of {_group} titers, vaccine strains"
                _subsection = f"{_chart['subtype']} ({_orientation})"
            add_htmls_to_docs.setdefault(_section, {}).setdefault(_subsection, {})[
                _chart["chart_type"]
            ] = _template.format(group=_group, orientation=_orientation)


analyze_titers_outputs = [
    # aggregated sera metadata
    rules.aggregate_sera_metadata.output.csv,
    # final processed titer data
    *expand(
        "results/final_titer_data/{group}_{output_type}.csv",
        group=groups_to_analyze,
        output_type=[
            "titers",
            "sera",
            "sera_multicohort",
            "viruses",
            "titers_summarized_by_virus",
        ],
    ),
    *expand(rules.process_final_titer_data.output.summary, group=groups_to_analyze),
    # titer summary plots
    *expand(
        rules.plot_titer_summaries.output.chart_htmls,
        group=groups_to_analyze,
        orientation=["vertical", "horizontal"],
    ),
]
