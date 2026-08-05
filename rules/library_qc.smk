"""QC of the viral library itself, as opposed to the neutralization titers.

The composition of the library pools, the balance of the re-pools made from them, and the
plates that infect each well with a single virus.

"""

# The three analyses this file can run, keyed by the `config.yml` section that configures
# each. Read with a default of empty rather than required, as a project need not do every
# kind of library QC: a section that is absent, or present with nothing under it, simply
# contributes no jobs. This is the only place these three keys are read, so it is the only
# place that default has to be applied. Everything *inside* one of these sections is still
# required, and the rules and scripts below raise on anything missing there.
analyze_pools = config.get("analyze_pools") or {}
analyze_repools = config.get("analyze_repools") or {}
analyze_single_well_infections = config.get("analyze_single_well_infections") or {}

# miscellaneous plate whose barcode counts each pool is analyzed from
pool_plates = {
    pool: pool_d["miscellaneous_plate"] for (pool, pool_d) in analyze_pools.items()
}

# the same for each re-pool, which is analyzed by `analyze_repool` below
repool_plates = {
    repool: repool_d["miscellaneous_plate"]
    for (repool, repool_d) in analyze_repools.items()
}

# The same again for each single-well analysis, except that one of those reads *several*
# plates and reports them together, so this maps to a list rather than to one plate.
single_well_plates = {
    single_well: single_well_d["miscellaneous_plates"]
    for (single_well, single_well_d) in analyze_single_well_infections.items()
}

# A re-pool is analyzed against the volumes that made the pool it came from, so that pool has
# to be one of the pools analyzed here. Without this the missing pool surfaces as a
# missing-input-file error naming a path, rather than as the configuration problem it is.
for _repool, _repool_d in analyze_repools.items():
    if _repool_d["previous_pool"] not in analyze_pools:
        raise ValueError(
            f"`analyze_repools` {_repool} names the previous pool "
            f"'{_repool_d['previous_pool']}', which is not an `analyze_pools` entry in "
            f"`config.yml`; the pools configured there are {sorted(analyze_pools)}"
        )

# The barcode counts of every plate in a single-well analysis are read against one
# whitelist and compared with each other, so the plates of a group have to agree on which
# viral library and neutralization standard set they used. Checked here rather than in the
# script, since it is what makes the rule able to name one of each as an input.
for _single_well, _plates in single_well_plates.items():
    if not _plates:
        raise ValueError(
            f"`analyze_single_well_infections` {_single_well} has no plates"
        )
    for _key in ["viral_library", "neut_standard_set"]:
        _values = {miscellaneous_plates[_plate][_key] for _plate in _plates}
        if len(_values) != 1:
            raise ValueError(
                f"the plates of `analyze_single_well_infections` {_single_well} have "
                f"{len(_values)} different values of '{_key}' ({sorted(_values)}), but "
                "are analyzed together and so must share one"
            )

# The three analyses here all write `results/library_qc/{name}_*`, so a name used in more
# than one section gives two rules the same output file. Snakemake reports that as an
# ambiguous-rule error naming only the file, so catch it here where the cause is obvious.
_names_seen = {}
for _section, _entries in [
    ("analyze_pools", analyze_pools),
    ("analyze_repools", analyze_repools),
    ("analyze_single_well_infections", analyze_single_well_infections),
]:
    for _name in _entries:
        if _name in _names_seen:
            raise ValueError(
                f"'{_name}' names both an `{_names_seen[_name]}` and an `{_section}` "
                "entry in `config.yml`; the two write the same files, so the names must "
                "differ"
            )
        _names_seen[_name] = _section


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
        pool="|".join(analyze_pools),
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        date=lambda wc: config["miscellaneous_plates"][pool_plates[wc.pool]]["date"],
        pool_config=lambda wc: analyze_pools[wc.pool],
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
            pool=analyze_repools[wc.repool]["previous_pool"]
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
        repool="|".join(analyze_repools),
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        date=lambda wc: config["miscellaneous_plates"][repool_plates[wc.repool]]["date"],
        repool_config=lambda wc: analyze_repools[wc.repool],
    script:
        "../scripts/analyze_repool.py"


rule analyze_single_well_infections:
    """Analyze plates that infect each well with a single virus of the library.

    Reports where each well's reads went, how much of a well is material that does not
    belong in it and which well on the plate that material came from, the reads matching
    neither whitelist, and the titer of each strain grown on its own relative to the
    neutralization standard.

    """
    input:
        counts=lambda wc: [
            rules.miscellaneous_plate_count_barcodes.output.counts.format(
                misc_plate=plate, well=well
            )
            for plate in single_well_plates[wc.single_well]
            for well in miscellaneous_plates[plate]["wells"]
        ],
        fates=lambda wc: [
            rules.miscellaneous_plate_count_barcodes.output.fates.format(
                misc_plate=plate, well=well
            )
            for plate in single_well_plates[wc.single_well]
            for well in miscellaneous_plates[plate]["wells"]
        ],
        invalid=lambda wc: [
            rules.miscellaneous_plate_count_barcodes.output.invalid.format(
                misc_plate=plate, well=well
            )
            for plate in single_well_plates[wc.single_well]
            for well in miscellaneous_plates[plate]["wells"]
        ],
        samples_csvs=lambda wc: [
            config["miscellaneous_plates"][plate]["samples_csv"]
            for plate in single_well_plates[wc.single_well]
        ],
        # one of each, the plates of a group having been checked above to share them
        viral_library=lambda wc: config["viral_libraries"][
            miscellaneous_plates[single_well_plates[wc.single_well][0]][
                "viral_library"
            ]
        ],
        neut_standard_set=lambda wc: config["neut_standard_sets"][
            miscellaneous_plates[single_well_plates[wc.single_well][0]][
                "neut_standard_set"
            ]
        ],
        # the library the strains assayed here are marked against, which is not the
        # whitelist they are counted against
        final_viral_library=lambda wc: config["viral_libraries"][
            analyze_single_well_infections[wc.single_well]["final_viral_library"]
        ],
    output:
        html="results/library_qc/{single_well}_analyze_single_well_infections.html",
        well_composition="results/library_qc/{single_well}_well_composition.csv",
    log:
        "results/logs/analyze_single_well_infections_{single_well}.txt",
    wildcard_constraints:
        single_well="|".join(analyze_single_well_infections),
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        plate_dates=lambda wc: {
            plate: config["miscellaneous_plates"][plate]["date"]
            for plate in single_well_plates[wc.single_well]
        },
        single_well_config=lambda wc: analyze_single_well_infections[wc.single_well],
    script:
        "../scripts/analyze_single_well_infections.py"


# Add the library QC reports to docs HTMLs generated by pipeline. Accumulated rather than
# assigned, so that the subsections do not depend on each other's presence or order.
_library_qc_docs = {}

if analyze_pools:
    _library_qc_docs["Composition of the library pools"] = {
        f"{_pool} ({config['miscellaneous_plates'][pool_plates[_pool]]['date']})": (
            rules.analyze_pool.output.html.format(pool=_pool)
        )
        for _pool in analyze_pools
    }

if analyze_repools:
    _library_qc_docs["Balance of the library re-pools"] = {
        f"{_repool} ({config['miscellaneous_plates'][repool_plates[_repool]]['date']})": (
            rules.analyze_repool.output.html.format(repool=_repool)
        )
        for _repool in analyze_repools
    }

if analyze_single_well_infections:
    _library_qc_docs["Single virus per well infections"] = {
        _single_well: rules.analyze_single_well_infections.output.html.format(
            single_well=_single_well
        )
        for _single_well in analyze_single_well_infections
    }

if _library_qc_docs:
    add_htmls_to_docs["Library quality-control"] = _library_qc_docs


# the HTMLs are not listed as they are already inputs to `build_docs`
library_qc_outputs = [
    *expand(rules.analyze_pool.output.repooling_math, pool=analyze_pools),
    *expand(rules.analyze_pool.output.dropped_strains, pool=analyze_pools),
    *expand(rules.analyze_repool.output.repooling_math, repool=analyze_repools),
    *expand(rules.analyze_repool.output.dropped_strains, repool=analyze_repools),
    *expand(
        rules.analyze_repool.output.subpool_repooling_math,
        repool=analyze_repools,
    ),
    *expand(
        rules.analyze_single_well_infections.output.well_composition,
        single_well=analyze_single_well_infections,
    ),
]
