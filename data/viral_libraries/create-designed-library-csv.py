"""Script to build the designed library CSV."""

import pandas as pd

input_csv = r"../../non-pipeline_analyses/library_design/construct_order/results/final_library/final_library.csv"
output_csv = "flu-seqneut-2026-barcode-to-strain-designed.csv"

# source columns read from the input library and carried through the build.
# `vaccine_annotation` drives `strain_type`/`vaccine_type` below, and
# `collection_date` is carried through from the library as-is.
columns = [
    "subtype",
    "derived_haplotype",
    "strain",
    "subclade",
    "genbank_accession",
    "vaccine_annotation",
    "passage_history_annotation",
    "shortname",
    "bloom_lab_plasmid_log_id",
    "barcode",
    "collection_date",
    "nt_sequence_HA_ectodomain",
    "protein_sequence_HA_ectodomain",
]

# final output column order (mirrors the previous round's designed library, using
# our source-column names where we chose not to rename, plus
# `passage_history_annotation` as an extra column).
output_columns = [
    "strain",
    "subtype",
    "strain_type",
    "vaccine_type",
    "subclade",
    "derived_haplotype",
    "genbank_accession",
    "shortname",
    "bloom_lab_plasmid_log_id",
    "passage_history_annotation",
    "barcode",
    "collection_date",
    "nt_sequence_HA_ectodomain",
    "protein_sequence_HA_ectodomain",
]

# `strain_type` value for non-vaccine (circulating) strains. Vaccine strains
# (vaccine_annotation == True in the library) are labeled "vaccine" instead.
circulating_strain_type = "circulating_2026"

def fill_derived_haplotype(df):
    """Validate that every row has a `derived_haplotype`.

    `derived_haplotype` is now supplied for all strains by the input library
    (including the testset strains, which previously required a manual fill), so
    this only fails fast if any value is still null.
    """
    df = df.copy()
    still_null = df.loc[df["derived_haplotype"].isnull(), "strain"].unique().tolist()
    if still_null:
        raise ValueError(
            f"derived_haplotype is null for strains: {still_null}"
        )
    return df


def fill_subclade(df):
    """Fill missing `subclade` values from `derived_haplotype` (text before ':').

    Fails fast if a row missing `subclade` also lacks a `derived_haplotype`, or if
    an existing `subclade` disagrees with the value implied by `derived_haplotype`.
    """
    df = df.copy()
    derived = df["derived_haplotype"].str.split(":").str[0]

    missing = df["subclade"].isnull()
    no_source = missing & df["derived_haplotype"].isnull()
    if no_source.any():
        raise ValueError(
            "Cannot derive subclade (missing derived_haplotype) for strains: "
            f"{df.loc[no_source, 'strain'].tolist()}"
        )

    check = (~missing) & df["derived_haplotype"].notnull()
    mismatch = check & (df["subclade"] != derived)
    if mismatch.any():
        raise ValueError(
            "Existing subclade disagrees with derived_haplotype-implied value:\n"
            f"{df.loc[mismatch, ['strain', 'subclade', 'derived_haplotype']]}"
        )

    df.loc[missing, "subclade"] = derived[missing]
    return df


def _is_vaccine(df):
    """Boolean mask of rows the library flags as vaccine strains."""
    return df["vaccine_annotation"].astype("string").str.strip() == "True"


def add_strain_type(df):
    """Set `strain_type` from the library's `vaccine_annotation`.

    Vaccine strains (vaccine_annotation == True) get "vaccine"; all others get
    the circulating label.
    """
    df = df.copy()
    df["strain_type"] = circulating_strain_type
    df.loc[_is_vaccine(df), "strain_type"] = "vaccine"
    return df


def add_vaccine_type(df):
    """Set `vaccine_type` to the passage (egg/cell) for vaccine strains, else null.

    Fails fast if a vaccine strain has no passage_history_annotation, since the
    egg/cell type is required to characterize a vaccine strain.
    """
    df = df.copy()
    df["vaccine_type"] = pd.NA
    vax = _is_vaccine(df)
    passage = df["passage_history_annotation"].astype("string").str.strip()
    missing = vax & (passage.isnull() | (passage == ""))
    if missing.any():
        raise ValueError(
            "Vaccine strain(s) missing passage_history_annotation "
            f"(needed for vaccine_type): {df.loc[missing, 'strain'].tolist()}"
        )
    df.loc[vax, "vaccine_type"] = passage[vax]
    return df


def add_collection_date(df):
    """Carry `collection_date` through from the input library unchanged.

    The library sources this from each strain's `latest_sequence` date; it is
    empty for strains without one (e.g. the older_* haplotype additions), which
    is preserved here rather than back-filled.
    """
    return df.copy()


def add_subtype_suffix(df):
    """Append ``_{subtype}`` to `strain` to match the previous library's naming.

    Must run AFTER any strain-based matching (e.g. `fill_derived_haplotype` and
    `add_collection_date`), which rely on the bare (unsuffixed) strain name.
    """
    df = df.copy()
    df["strain"] = df["strain"] + "_" + df["subtype"]
    return df


def tally_fields(df):
    """Print unique-value and null-value tallies for each column of `df`."""
    if df["strain"].isnull().any():
        raise ValueError(
            "Found rows with null strain, cannot tally other columns by strain:\n"
            f"{df[df['strain'].isnull()]}"
        )

    # sequence fields are too long to usefully print example values for
    sequence_fields = {
        "barcode",
        "nt_sequence_HA_ectodomain",
        "protein_sequence_HA_ectodomain",
    }

    n_rows = len(df)
    for col in df.columns:
        n_unique = df[col].nunique(dropna=True)
        null_strains = df.loc[df[col].isnull(), "strain"].tolist()

        print(f"=== {col} ===")
        print(f"unique values: {n_unique} / {n_rows}")
        if col not in sequence_fields:
            first_unique = df[col].dropna().unique()[:10].tolist()
            print(f"first 10 unique values: {first_unique}")
        print(f"null values: {len(null_strains)} / {n_rows}")
        print(f"strains with null values: {null_strains}")
        print()


if __name__ == "__main__":
    df = pd.read_csv(input_csv)[columns]
    df = fill_derived_haplotype(df)
    df = fill_subclade(df)
    df = add_strain_type(df)
    df = add_vaccine_type(df)
    df = add_collection_date(df)
    df = add_subtype_suffix(df)
    df = df[output_columns]
    tally_fields(df)
    df.to_csv(output_csv, index=False)
    print(f"Wrote {len(df)} rows to {output_csv}")
