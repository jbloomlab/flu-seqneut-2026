"""Table of the circulating strains with the lowest titers, for each subtype.

Strains are ranked by their median titer over one cohort of the per-virus titer summary,
and each is given its median titer relative to the subtype's reference strain, which is
listed alongside them. The ratio is of the medians the summary reports, not a median of
per-serum ratios.

"""

import sys

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

n_strains = snakemake.params.config["n_strains"]
cohort = snakemake.params.config["cohort"]
circulating_strain_type = snakemake.params.circulating_strain_type
reference_haplotypes = snakemake.params.reference_haplotypes

VIRUS_COLUMNS = [
    "virus",
    "subtype",
    "strain_type",
    "subclade",
    "derived_haplotype",
    "virus_collection_date",
]
SUMMARY_COLUMNS = ["virus", "cohort", "n_sera", "median_titer", "titer_q1", "titer_q3"]
OUTPUT_COLUMNS = [
    "subtype",
    "row_type",
    "rank",
    "virus",
    "derived_haplotype",
    "subclade",
    "virus_collection_date",
    "n_sera",
    "median_titer",
    "titer_q1",
    "titer_q3",
    "median_titer_relative_to_reference",
]

summary = pd.read_csv(
    snakemake.input.titers_summarized, usecols=SUMMARY_COLUMNS, dtype={"virus": str}
)
viruses = pd.read_csv(
    snakemake.input.viruses, usecols=VIRUS_COLUMNS, dtype={"virus": str}
)

if cohort not in set(summary["cohort"]):
    raise ValueError(
        f"cohort {cohort!r} is not in the titer summary; cohorts are "
        f"{sorted(set(summary['cohort']))}"
    )
summary = summary[summary["cohort"] == cohort].drop(columns="cohort")

strains = summary.merge(viruses, on="virus", how="left", validate="one_to_one")
if strains["subtype"].isnull().any():
    raise ValueError(
        "viruses in the titer summary but not the viruses CSV: "
        f"{strains.loc[strains['subtype'].isnull(), 'virus'].tolist()}"
    )

tables = []
for subtype, reference_haplotype in reference_haplotypes.items():
    subtype_strains = strains[strains["subtype"] == subtype]
    if reference_haplotype is None:
        raise ValueError(f"{subtype} has no reference strain (`draw_titer_line`)")
    reference = subtype_strains[
        subtype_strains["derived_haplotype"] == reference_haplotype
    ]
    if len(reference) != 1:
        raise ValueError(
            f"reference {reference_haplotype!r} matches {len(reference)} {subtype} "
            "strains, but must match exactly one"
        )
    eligible = subtype_strains[
        subtype_strains["strain_type"] == circulating_strain_type
    ]
    if len(eligible) < n_strains:
        raise ValueError(
            f"{subtype} has {len(eligible)} {circulating_strain_type} strains, fewer "
            f"than the {n_strains} to tabulate"
        )

    # ties in the median titer are broken by name so the table is reproducible
    lowest = eligible.sort_values(["median_titer", "virus"]).head(n_strains)
    lowest = lowest.assign(row_type="lowest_titer", rank=range(1, n_strains + 1))
    table = pd.concat([reference.assign(row_type="reference"), lowest])
    table["median_titer_relative_to_reference"] = (
        table["median_titer"] / reference["median_titer"].item()
    )
    tables.append(table)
    print(
        f"{subtype}: reference {reference['virus'].item()} has median titer "
        f"{reference['median_titer'].item()}; lowest of {len(eligible)} "
        f"{circulating_strain_type} strains is {lowest['virus'].iloc[0]} at "
        f"{lowest['median_titer'].iloc[0]}"
    )

# `rank` is empty for the reference, which is not ranked among the circulating strains
output = pd.concat(tables).astype({"rank": "Int64"})[OUTPUT_COLUMNS]
output.to_csv(snakemake.output.csv, index=False, float_format="%.5g")
print(f"Wrote {len(output)} rows to {snakemake.output.csv}")
