"""QC of the viral library itself, as opposed to the neutralization titers.

Currently the composition of the library pools and the balance of the re-pools made from
them; the single-well infections will also go here.

"""

# miscellaneous plate whose barcode counts each pool is analyzed from
pool_plates = {
    pool: pool_d["miscellaneous_plate"]
    for (pool, pool_d) in config["analyze_pools"].items()
}

# the same for each re-pool, which is analyzed by `analyze_repool` below
repool_plates = {
    repool: repool_d["miscellaneous_plate"]
    for (repool, repool_d) in config["analyze_repools"].items()
}

# `analyze_pool` and `analyze_repool` both write `results/library_qc/{name}_*.csv`, so a
# name used in both sections gives two rules the same output file. Snakemake reports that
# as an ambiguous-rule error naming only the file, so catch it here where the cause is
# obvious.
_shared_names = set(config["analyze_pools"]) & set(config["analyze_repools"])
if _shared_names:
    raise ValueError(
        f"{sorted(_shared_names)} names both an `analyze_pools` and an `analyze_repools` "
        "entry in `config.yml`; the two write the same files, so the names must differ"
    )


rule analyze_pool:
    """Analyze the composition of a library pool and get the re-pooling volumes."""
    input:
        counts=lambda wc: expand(
            rules.miscellaneous_plate_count_barcodes.output.counts,
            misc_plate=pool_plates[wc.pool],
            well=miscellaneous_plates[pool_plates[wc.pool]]["wells"],
        ),
        fates=lambda wc: expand(
            rules.miscellaneous_plate_count_barcodes.output.fates,
            misc_plate=pool_plates[wc.pool],
            well=miscellaneous_plates[pool_plates[wc.pool]]["wells"],
        ),
        invalid=lambda wc: expand(
            rules.miscellaneous_plate_count_barcodes.output.invalid,
            misc_plate=pool_plates[wc.pool],
            well=miscellaneous_plates[pool_plates[wc.pool]]["wells"],
        ),
        samples_csv=lambda wc: config["miscellaneous_plates"][pool_plates[wc.pool]][
            "samples_csv"
        ],
        viral_library=lambda wc: config["viral_libraries"][
            miscellaneous_plates[pool_plates[wc.pool]]["viral_library"]
        ],
        neut_standard_set=lambda wc: config["neut_standard_sets"][
            miscellaneous_plates[pool_plates[wc.pool]]["neut_standard_set"]
        ],
    output:
        html="results/library_qc/{pool}_analyze_pool.html",
        repooling_math="results/library_qc/{pool}_repooling_math.csv",
        dropped_strains="results/library_qc/{pool}_dropped_strains.csv",
    log:
        "results/logs/analyze_pool_{pool}.txt",
    wildcard_constraints:
        pool="|".join(config["analyze_pools"]),
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        date=lambda wc: config["miscellaneous_plates"][pool_plates[wc.pool]]["date"],
        pool_config=lambda wc: config["analyze_pools"][wc.pool],
    script:
        "../scripts/analyze_pool.py"


rule analyze_repool:
    """Analyze how well a balanced re-pool came out and get corrective volumes."""
    input:
        counts=lambda wc: expand(
            rules.miscellaneous_plate_count_barcodes.output.counts,
            misc_plate=repool_plates[wc.repool],
            well=miscellaneous_plates[repool_plates[wc.repool]]["wells"],
        ),
        fates=lambda wc: expand(
            rules.miscellaneous_plate_count_barcodes.output.fates,
            misc_plate=repool_plates[wc.repool],
            well=miscellaneous_plates[repool_plates[wc.repool]]["wells"],
        ),
        samples_csv=lambda wc: config["miscellaneous_plates"][
            repool_plates[wc.repool]
        ]["samples_csv"],
        viral_library=lambda wc: config["viral_libraries"][
            miscellaneous_plates[repool_plates[wc.repool]]["viral_library"]
        ],
        neut_standard_set=lambda wc: config["neut_standard_sets"][
            miscellaneous_plates[repool_plates[wc.repool]]["neut_standard_set"]
        ],
        # the volumes that made this pool, needed to recover each strain's stock titer
        previous_repooling_math=lambda wc: rules.analyze_pool.output.repooling_math.format(
            pool=config["analyze_repools"][wc.repool]["previous_pool"]
        ),
    output:
        html="results/library_qc/{repool}_analyze_repool.html",
        repooling_math="results/library_qc/{repool}_repooling_math.csv",
        dropped_strains="results/library_qc/{repool}_dropped_strains.csv",
        # One pipetting CSV per subpool being remade from the strain stocks, written into
        # this directory. Which subpools those are comes from the config and so differs
        # between re-pools, while `output` has to be static and cannot consult the
        # `{repool}` wildcard -- hence a directory rather than a file each.
        subpool_repooling_math=directory(
            "results/library_qc/{repool}_subpool_repooling_math"
        ),
    log:
        "results/logs/analyze_repool_{repool}.txt",
    wildcard_constraints:
        repool="|".join(config["analyze_repools"]),
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        date=lambda wc: config["miscellaneous_plates"][repool_plates[wc.repool]]["date"],
        repool_config=lambda wc: config["analyze_repools"][wc.repool],
    script:
        "../scripts/analyze_repool.py"


# Add the library QC reports to docs HTMLs generated by pipeline. Accumulated rather than
# assigned, so that the subsections do not depend on each other's presence or order.
_library_qc_docs = {}

if config["analyze_pools"]:
    _library_qc_docs["Composition of the library pools"] = {
        f"{_pool} ({config['miscellaneous_plates'][pool_plates[_pool]]['date']})": (
            rules.analyze_pool.output.html.format(pool=_pool)
        )
        for _pool in config["analyze_pools"]
    }

if config["analyze_repools"]:
    _library_qc_docs["Balance of the library re-pools"] = {
        f"{_repool} ({config['miscellaneous_plates'][repool_plates[_repool]]['date']})": (
            rules.analyze_repool.output.html.format(repool=_repool)
        )
        for _repool in config["analyze_repools"]
    }

if _library_qc_docs:
    add_htmls_to_docs["Library quality-control"] = _library_qc_docs


# the HTMLs are not listed as they are already inputs to `build_docs`
library_qc_outputs = [
    *expand(rules.analyze_pool.output.repooling_math, pool=config["analyze_pools"]),
    *expand(rules.analyze_pool.output.dropped_strains, pool=config["analyze_pools"]),
    *expand(
        rules.analyze_repool.output.repooling_math, repool=config["analyze_repools"]
    ),
    *expand(
        rules.analyze_repool.output.dropped_strains, repool=config["analyze_repools"]
    ),
    *expand(
        rules.analyze_repool.output.subpool_repooling_math,
        repool=config["analyze_repools"],
    ),
]
