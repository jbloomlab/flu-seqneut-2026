"""Build alignment, metadata, and titers for nextstrain-prot-titers-tree.

This script prepares input files for the nextstrain-prot-titers-tree pipeline:
- Alignment FASTA with protein sequences (with optional prefix for H1N1)
- Metadata TSV with strain information and titer summary columns
- Titers TSV with per-serum titer data for tree overlay (when titers configured)
- YAML of the per-cohort titer colorings, read via `color_by_metadata_file`, so the
  cohorts do not have to be enumerated in the configuration file

The titers TSV has one row per serum and strain, with the serum's cohort as a column
that the measurements panel groups by. The colorings instead cover every cohort of the
summarized titers, including the "All" cohort spanning all sera.
"""

import datetime
import sys

import numpy as np
import pandas as pd
import yaml

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

subtypes = snakemake.params.subtypes
subtypes_w_titers = snakemake.params.subtypes_w_titers
circulating_strain_type = snakemake.params.circulating_strain_type
recent_vaccine_strains = snakemake.params.recent_vaccine_strains
prefix_alignment = snakemake.params.prefix_alignment
titer_cutoffs = snakemake.params.titer_cutoffs

frac_below_cols = [f"frac_w_titer_below_{cutoff}" for cutoff in titer_cutoffs]

viruses = pd.read_csv(snakemake.input.viral_libraries_csv)[
    [
        "strain",
        "subtype",
        "derived_haplotype",
        "strain_type",
        "protein_sequence_HA_ectodomain",
        "subclade",
        "collection_date",
    ]
].drop_duplicates()

assert len(viruses) == viruses["strain"].nunique(), "Duplicate strain entries found"

# Validate that recent_vaccine_strains are in the viral library
if recent_vaccine_strains:
    assert set(recent_vaccine_strains).issubset(
        viruses["strain"]
    ), f"recent_vaccine_strains not found in viral library: {set(recent_vaccine_strains) - set(viruses['strain'])}"

# Filter to circulating strains and recent vaccine strains
df = viruses[
    (viruses["strain_type"] == circulating_strain_type)
    | viruses["strain"].isin(recent_vaccine_strains)
].copy()

# Relabel strain_type for vaccine strains using the label from recent_vaccine_strains dict
if recent_vaccine_strains:
    df["strain_type"] = df.apply(
        lambda x: recent_vaccine_strains.get(x["strain"], x["strain_type"]),
        axis=1,
    )

print(
    f"{len(df)=} of {len(viruses)} are {circulating_strain_type=} or in {recent_vaccine_strains=}"
)

# Ensure collection_date is in valid format (numerical year)
year = datetime.datetime.now(datetime.timezone.utc).year
if all((df["collection_date"] > year - 100) & (df["collection_date"] < year + 1)):
    df = df.rename(columns={"collection_date": "date"})
else:
    raise ValueError(f"Not valid numerical dates in {df['collection_date'].tolist()}")

# Process titer data if any tree displays titers
if subtypes_w_titers:
    print(f"\nProcessing titer data for {subtypes_w_titers=}")

    # Read summarized titers (for metadata columns)
    summarized_titers = pd.read_csv(snakemake.input.summarized_titers_csv).rename(
        columns={"virus": "strain"}
    )
    print(f"Read {len(summarized_titers)=} summarized titer rows")

    # Validate frac_below_cols exist
    assert set(frac_below_cols).issubset(
        summarized_titers.columns
    ), f"{frac_below_cols=} not all in {summarized_titers.columns.tolist()=}"

    # Filter summarized titers to only include strains in tree (may exclude historical vaccines)
    strains_not_in_tree = set(summarized_titers["strain"]) - set(df["strain"])
    if strains_not_in_tree:
        print(
            f"Filtering out {len(strains_not_in_tree)} strains not in tree: {strains_not_in_tree}"
        )
        summarized_titers = summarized_titers[
            summarized_titers["strain"].isin(df["strain"])
        ]
        print(f"After filtering: {len(summarized_titers)=} rows")

    # cohorts in the order they appear in the summarized titers ("All" first)
    cohorts = summarized_titers["cohort"].drop_duplicates().tolist()
    print(f"{cohorts=}")

    # The colorings the tree offers, one set per cohort. Generated here rather than
    # configured because the cohorts are set by the titer data, not by a choice.
    colorings = {}
    for cohort in cohorts:
        colorings[f"median_titer_{cohort}_sera"] = {
            "scale_type": "viridis_r_log",
            "title": f"median titer ({cohort} sera)",
        }
        for cutoff in titer_cutoffs:
            colorings[f"frac_w_titer_below_{cutoff}_{cohort}_sera"] = {
                "scale_type": "viridis_linear",
                "title": f"fraction of {cohort} sera with titer < {cutoff}",
            }

    # Pivot median_titer and frac_below columns by cohort and merge into metadata
    for col in ["median_titer"] + frac_below_cols:
        assert col not in df.columns, f"{col} already in df columns"
        pivoted = (
            summarized_titers.assign(
                cohort_col=lambda x, col=col: f"{col}_" + x["cohort"] + "_sera"
            )
            .pivot_table(index="strain", values=col, columns="cohort_col")
            .reset_index()
        )
        df = df.merge(pivoted, on="strain", how="left", validate="one_to_one")

    # Read individual titers (for per-serum titers TSV)
    titers = (
        pd.read_csv(snakemake.input.titers_csv)[["serum", "virus", "titer"]]
        .rename(columns={"virus": "strain"})
        .assign(log2_titer=lambda x: np.log2(x["titer"]))
        .drop(columns="titer")
    )
    print(f"Read {len(titers)=} individual titer rows")

    # Filter titers to only include strains in tree
    titers = titers[titers["strain"].isin(df["strain"])]
    print(f"After filtering to tree strains: {len(titers)=} rows")

    # Read sera metadata (one row per serum). The columns kept here are what the titers
    # TSV carries, so they are the ones the measurements panel can group by.
    sera = pd.read_csv(snakemake.input.sera_metadata_csv)[
        ["serum", "cohort", "serum_collection_date", "age_numeric", "sex"]
    ]
    assert len(sera) == sera["serum"].nunique(), "Duplicate sera in metadata"
    assert set(titers["serum"]).issubset(
        sera["serum"]
    ), f"Titers have sera not in metadata: {set(titers['serum']) - set(sera['serum'])}"
    # every cohort the measurements panel groups by needs one of the colorings above
    assert set(sera["cohort"]).issubset(
        cohorts
    ), f"Cohorts lacking colorings: {set(sera['cohort']) - set(cohorts)}"

    # nextstrain-prot-titers-tree allows only one titer per (serum, strain)
    assert len(titers) == len(
        titers[["serum", "strain"]].drop_duplicates()
    ), "Duplicate serum-strain pairs in titers"

    titers = titers.merge(sera, on="serum", validate="many_to_one")
    print(
        f"Titers span {titers['serum'].nunique()} sera in "
        f"{titers['cohort'].nunique()} cohorts"
    )

# Process each subtype
for subtype in subtypes:
    print(f"\nProcessing {subtype=}")
    subtype_df = df[df["subtype"] == subtype].drop(columns="subtype")
    print(f"{len(subtype_df)=} of {len(df)=} are {subtype=}")

    if len(subtype_df) == 0:
        raise ValueError(f"No strains found for {subtype=}")

    # Remove subtype suffix from strain names if present (e.g., "_H3N2")
    # Keep original names for filtering titers before renaming
    original_strain_names = set(subtype_df["strain"])
    strain_rename = {
        s: (s[: -len(subtype) - 1] if s.endswith(f"_{subtype}") else s)
        for s in subtype_df["strain"]
    }
    subtype_df["strain"] = subtype_df["strain"].map(strain_rename)
    assert len(subtype_df) == subtype_df["strain"].nunique()

    alignment_file = snakemake.output[f"alignment_{subtype}"]
    metadata_file = snakemake.output[f"metadata_{subtype}"]

    print(f"Writing alignment to {alignment_file=}")
    with open(alignment_file, "w") as f:
        for tup in subtype_df.itertuples():
            seq = prefix_alignment[subtype] + tup.protein_sequence_HA_ectodomain
            f.write(f">{tup.strain}\n{seq}\n")

    metadata = subtype_df.drop(columns=["protein_sequence_HA_ectodomain"])
    print(f"Writing metadata to {metadata_file=}")
    metadata.to_csv(metadata_file, index=False, sep="\t", float_format="%.6g")

    # Write titers TSV and the colorings of them if this tree displays titers
    if subtype in subtypes_w_titers:
        titers_file = snakemake.output[f"titers_{subtype}"]
        # Filter using original strain names (before suffix removal)
        subtype_titers = titers[titers["strain"].isin(original_strain_names)]
        print(
            f"{len(subtype_titers)} of {len(titers)} titer rows are for {subtype=} strains "
            f"({subtype_titers[['serum', 'strain']].drop_duplicates().shape[0]} unique serum-strain pairs)"
        )

        print(f"Writing titers to {titers_file=}")
        subtype_titers.assign(
            strain=lambda x, strain_rename=strain_rename: x["strain"].map(strain_rename)
        ).to_csv(titers_file, sep="\t", index=False, float_format="%.6g")

        colorings_file = snakemake.output[f"color_by_metadata_{subtype}"]
        missing_cols = [c for c in colorings if c not in metadata.columns]
        if missing_cols:
            raise ValueError(f"Metadata for {subtype=} lacks colorings {missing_cols}")
        print(f"Writing {len(colorings)} colorings to {colorings_file=}")
        with open(colorings_file, "w") as f:
            f.write(
                "# Titer colorings of the tree, one set per serum cohort, written by\n"
                "# `nextstrain_prot_titers_tree_alignment_and_metadata.py`.\n"
            )
            yaml.dump(colorings, f, sort_keys=False)

print("\nDone!")
