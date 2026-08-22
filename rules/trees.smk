"""Building of the `nextstrain-prot-titers-tree` protein trees for each subtype."""

import os

# subtypes whose tree displays titers (the others just get an alignment and metadata)
_subtypes_w_titers = [
    subtype
    for subtype in config["subtypes"]
    if config["nextstrain-prot-titers-tree_config"][subtype]["titers"]
]

_titers_from = config["nextstrain-prot-titers-tree_titers_from"]
if _subtypes_w_titers and _titers_from not in groups_to_analyze:
    raise ValueError(
        f"`nextstrain-prot-titers-tree_titers_from` is {_titers_from!r}, which is not in "
        "`groups_to_analyze`, so its final titer data is never built"
    )


rule nextstrain_prot_titers_tree_alignment_and_metadata:
    """Build alignment, metadata, and titers TSV used by `nextstrain-prot-titers-tree`."""
    input:
        viral_libraries_csv=config["viral_libraries"][
            config["nextstrain-prot-titers-tree_viral_library"]
        ],
        # Only include titer inputs if titers are configured for any tree
        summarized_titers_csv=(
            f"results/final_titer_data/{_titers_from}_titers_summarized_by_virus.csv"
            if _subtypes_w_titers
            else []
        ),
        titers_csv=(
            f"results/final_titer_data/{_titers_from}_titers.csv"
            if _subtypes_w_titers
            else []
        ),
        sera_metadata_csv=(
            f"results/final_titer_data/{_titers_from}_sera.csv"
            if _subtypes_w_titers
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
        # titers and the colorings derived from them only for trees that display titers
        **{
            f"titers_{subtype}": f"results/nextstrain-prot-titers-tree/{subtype}/titers.tsv"
            for subtype in _subtypes_w_titers
        },
        **{
            f"color_by_metadata_{subtype}": config[
                "nextstrain-prot-titers-tree_config"
            ][subtype]["color_by_metadata_file"]
            for subtype in _subtypes_w_titers
        },
    log:
        "results/logs/nextstrain_prot_titers_tree_alignment_and_metadata.txt",
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        subtypes=config["subtypes"],
        subtypes_w_titers=_subtypes_w_titers,
        circulating_strain_type=config["circulating_strain_type"],
        recent_vaccine_strains=config["recent_vaccine_strains"],
        prefix_alignment=config["nextstrain-prot-titers-tree_prefix_alignment"],
        titer_cutoffs=config["titer_cutoffs"],
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
for _tree_config in config["nextstrain-prot-titers-tree_config"].values():
    trees_outputs.append(_tree_config["auspice_json"])
    if _tree_config["titers"]:  # only add measurements.json if titers configured
        trees_outputs.append(
            os.path.splitext(_tree_config["auspice_json"])[0] + "_measurements.json"
        )
