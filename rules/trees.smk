"""Building of the `nextstrain-prot-titers-tree` protein trees for each subtype."""

import os

# Check if any subtype has titers configured (non-null titers key in config)
_any_tree_has_titers = any(
    config["nextstrain-prot-titers-tree_config"][subtype].get("titers")
    for subtype in config["subtypes"]
)


rule nextstrain_prot_titers_tree_alignment_and_metadata:
    """Build alignment, metadata, and titers TSV used by `nextstrain-prot-titers-tree`."""
    input:
        viral_libraries_csv=config["viral_libraries"][
            config["nextstrain-prot-titers-tree_viral_library"]
        ],
        # Only include titer inputs if titers are configured for any tree
        summarized_titers_csv=(
            f"results/final_titer_data/{config['nextstrain-prot-titers-tree_titers_from']}_titers_summarized_by_virus.csv"
            if _any_tree_has_titers
            else []
        ),
        titers_csv=(
            f"results/final_titer_data/{config['nextstrain-prot-titers-tree_titers_from']}_titers.csv"
            if _any_tree_has_titers
            else []
        ),
        sera_metadata_csv=(
            f"results/final_titer_data/{config['nextstrain-prot-titers-tree_titers_from']}_sera_multicohort.csv"
            if _any_tree_has_titers
            else []
        ),
    output:
        **{
            f"alignment_{subtype}": f"results/nextstrain-prot-titers-tree/{subtype}/alignment.fa"
            for subtype in config["subtypes"]
        },
        **{
            f"metadata_{subtype}": f"results/nextstrain-prot-titers-tree/{subtype}/metadata.tsv"
            for subtype in config["subtypes"]
        },
        # Only output titers TSV for subtypes that have titers configured
        **{
            f"titers_{subtype}": f"results/nextstrain-prot-titers-tree/{subtype}/titers.tsv"
            for subtype in config["subtypes"]
            if config["nextstrain-prot-titers-tree_config"][subtype].get("titers")
        },
    log:
        "results/logs/nextstrain_prot_titers_tree_alignment_and_metadata.txt",
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        subtypes=config["subtypes"],
        circulating_strain_type=config["circulating_strain_type"],
        recent_vaccine_strains=config["recent_vaccine_strains"],
        prefix_alignment=config["nextstrain-prot-titers-tree_prefix_alignment"],
        frac_below_cols=[
            f"frac_w_titer_below_{cutoff}" for cutoff in config["titer_cutoffs"]
        ],
        serum_cohorts_for_tree=(
            config["serum_cohorts_for_tree"] if _any_tree_has_titers else []
        ),
        has_titers=_any_tree_has_titers,
    script:
        "../scripts/nextstrain_prot_titers_tree_alignment_and_metadata.py"


# run the nextstrain-prot-titers-tree submodule on each lineage
for subtype in config["subtypes"]:
    module_name = f"nextstrain-prot-titers-tree_{subtype}"

    module:
        name: module_name
        snakefile:
            "../nextstrain-prot-titers-tree/Snakefile"
        config:
            config["nextstrain-prot-titers-tree_config"][subtype]

    use rule * from module_name as module_name*


# auspice JSONs (and their measurements JSONs) built by the module above
trees_outputs = []
for _tree_config in config.get("nextstrain-prot-titers-tree_config", {}).values():
    trees_outputs.append(_tree_config["auspice_json"])
    if _tree_config.get("titers"):  # only add measurements.json if titers configured
        trees_outputs.append(
            os.path.splitext(_tree_config["auspice_json"])[0] + "_measurements.json"
        )
