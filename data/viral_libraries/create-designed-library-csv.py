"""Script to build the designed library CSV."""

import pandas as pd

input_csv = r"../../non-pipeline_analyses/library_design/construct_order/results/final_library/final_library.csv"
output_csv = "flu-seqneut-2026-barcode-to-strain-designed.csv"

# source columns read from the input library and carried through the build.
# `collection_date` and `vaccine_annotation` are intentionally NOT read here:
# `collection_date` is (re)built from the strain name below, and
# `vaccine_annotation` is dropped in favor of the `strain_type`/`vaccine_type`
# placeholders below.
columns = [
    "subtype",
    "derived_haplotype",
    "strain",
    "subclade",
    "genbank_accession",
    "passage_history_annotation",
    "shortname",
    "bloom_lab_plasmid_log_id",
    "barcode",
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

# Placeholder `strain_type` applied to every row.
# NOT YET VALID: this does not distinguish vaccine strains from circulating
# strains -- every row is marked circulating. It needs a proper vaccine-strain
# assignment before the output can be used.
strain_type_placeholder = "circulating_2026"

# manually fill missing derived haplotypes for these strains
manual_derived_haplotypes = {
    "A/Darwin/1454/2025": "K:F195Y",
    "A/Singapore/GP20238/2024": "J.2.4:F195Y",
    "A/Missouri/11/2025": "D.3.1:Q223R",
    "A/Nebraska/34/2026": "D.3.1.1:R205K",
    "A/Rhode_Island/11/2026": "D.3.1.1:D127N,D139N",
}


def fill_derived_haplotype(df):
    """Fill missing `derived_haplotype` values from `manual_derived_haplotypes`.

    Fails fast if a manual key matches no strain, if a manually specified strain
    already has a non-null `derived_haplotype` (would be an unexpected conflict),
    or if any `derived_haplotype` remains null after filling.
    """
    df = df.copy()
    for strain, haplotype in manual_derived_haplotypes.items():
        mask = df["strain"] == strain
        if not mask.any():
            raise ValueError(
                f"manual_derived_haplotypes strain {strain!r} matches no library row"
            )
        existing = df.loc[mask, "derived_haplotype"]
        if existing.notnull().any():
            raise ValueError(
                f"manual_derived_haplotypes strain {strain!r} already has a "
                f"derived_haplotype: {existing.dropna().unique().tolist()}"
            )
        df.loc[mask, "derived_haplotype"] = haplotype

    still_null = df.loc[df["derived_haplotype"].isnull(), "strain"].unique().tolist()
    if still_null:
        raise ValueError(
            f"derived_haplotype still null after manual fill for strains: {still_null}"
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


def add_strain_type(df):
    """Add the placeholder `strain_type` column.

    NOT YET VALID -- see the `strain_type_placeholder` note above: every row is
    marked circulating and vaccine strains are not yet distinguished.
    """
    df = df.copy()
    df["strain_type"] = strain_type_placeholder
    return df


def add_vaccine_type(df):
    """Add the placeholder `vaccine_type` column, null for every row.

    NOT YET VALID: cell/egg vaccine typing is not yet assigned.
    """
    df = df.copy()
    df["vaccine_type"] = pd.NA
    return df


def add_collection_date(df):
    """Add a rough `collection_date` = the 4-digit year parsed from the strain name.

    NOT YET REFINED: this is only the year embedded in the strain name (the
    previous round used a decimal year such as 2025.76); it needs to be replaced
    with an actual collection date. Fails fast if any strain has no parseable
    trailing year (expects names ending in ``/YYYY``). Must run BEFORE the subtype
    suffix is appended to `strain`.
    """
    df = df.copy()
    year = df["strain"].str.extract(r"/(\d{4})$", expand=False)
    missing = year.isnull()
    if missing.any():
        raise ValueError(
            "Could not parse a trailing /YYYY year from strain names: "
            f"{df.loc[missing, 'strain'].tolist()}"
        )
    df["collection_date"] = year.astype(int)
    return df


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
