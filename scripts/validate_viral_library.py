"""Validate a viral library CSV against the checks configured in `config.yml`.

Writes a per-column report of the library to both the log and the output file, and
raises `ValueError` on the first check that fails.

"""

import sys

import Bio.Seq
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# the columns this script names directly, and so cannot check a library without; every
# other entry in `columns` is purely declarative and may change without touching this file
SCRIPT_REQUIRED_COLS = {
    "strain",
    "subtype",
    "strain_type",
    "barcode",
    "collection_date",
    "nt_sequence_HA_ectodomain",
    "protein_sequence_HA_ectodomain",
}
# the values `unique` may take, and the keys a `columns` entry may have
UNIQUE_VALUES = {"barcode", "strain"}
SPEC_KEYS = {"nulls", "unique", "values", "non_null_only_when"}
# the keys each `by_subtype` entry must have
SUBTYPE_KEYS = {"protein_lengths", "protein_startswith", "protein_endswith"}

validations = snakemake.params.validations
col_specs = validations["columns"]
by_subtype = validations["by_subtype"]
barcode_length = validations["barcode_length"]


def require(condition, message):
    """Raise `ValueError` describing a malformed `viral_library_validations`."""
    if not condition:
        raise ValueError(f"`viral_library_validations` is invalid: {message}")


# the spec is checked before the CSV is read, so that a configuration error is never
# reported as though it were a problem with the data
require(
    SCRIPT_REQUIRED_COLS <= set(col_specs),
    "`columns` must list every column this script checks directly, but is missing "
    f"{sorted(SCRIPT_REQUIRED_COLS - set(col_specs))}",
)
require(
    "values" in col_specs["subtype"],
    "`columns: subtype` must declare `values`, which `by_subtype` is checked against",
)
for col, spec in col_specs.items():
    require(
        not set(spec) - SPEC_KEYS,
        f"`columns: {col}` has unknown keys "
        f"{sorted(set(spec) - SPEC_KEYS)}, expected some of {sorted(SPEC_KEYS)}",
    )
    require("nulls" in spec, f"`columns: {col}` must declare `nulls`")
    require(
        "unique" not in spec or spec["unique"] in UNIQUE_VALUES,
        f"`columns: {col}` has `unique: {spec.get('unique')}`, expected one of "
        f"{sorted(UNIQUE_VALUES)} or no `unique` at all",
    )
    for when_col in spec.get("non_null_only_when", {}):
        require(
            when_col in col_specs,
            f"`columns: {col}` has `non_null_only_when` naming {when_col!r}, which "
            "`columns` does not list",
        )

for subtype, spec in by_subtype.items():
    require(
        set(spec) == SUBTYPE_KEYS,
        f"`by_subtype: {subtype}` has keys {sorted(spec)}, expected exactly "
        f"{sorted(SUBTYPE_KEYS)}",
    )
    require(
        spec["protein_lengths"],
        f"`by_subtype: {subtype}: protein_lengths` must list at least one length",
    )

csv = snakemake.input.csv
print(f"Validating {csv}")
# every column is read as a string so that no identifier is coerced to a number;
# `collection_date` is converted to a number where it is checked
df = pd.read_csv(csv, dtype=str)
print(f"Read {len(df)} rows\n")

per_barcode_cols = [
    c for c, spec in col_specs.items() if spec.get("unique") == "barcode"
]


def describe(col):
    """Human-readable summary of the configured checks on column `col`."""
    spec = col_specs.get(col)
    if spec is None:
        return "NOT VALIDATED"
    parts = ["nulls allowed" if spec["nulls"] else "non-null"]
    if "unique" in spec:
        parts.append(f"unique per {spec['unique']}")
    if "values" in spec:
        parts.append(f"values {spec['values']}")
    if "non_null_only_when" in spec:
        ((when_col, when_val),) = spec["non_null_only_when"].items()
        parts.append(f"non-null only when {when_col} == {when_val}")
    return "; ".join(parts)


def summarize():
    """Emit lines describing the library and every column in it, listed or not."""
    lines = [
        f"=== summary of {csv} ===",
        f"{len(df)} barcodes for {df['strain'].nunique()} strains",
    ]
    strains = df.drop_duplicates("strain")
    for subtype, subtype_strains in strains.groupby("subtype"):
        by_type = subtype_strains["strain_type"].value_counts().sort_index()
        counts = ", ".join(f"{n} {t}" for t, n in by_type.items())
        lines.append(
            f"  {subtype}: {len(subtype_strains)} strains ({counts}), "
            f"{int((df['subtype'] == subtype).sum())} barcodes"
        )

    # one row per column of the CSV, in the order the CSV has them, so that a column the
    # configuration does not mention still shows up rather than passing unnoticed
    header = ("column", "null", "non-null", "distinct", "max per strain", "checks")
    rows = [
        (
            col,
            str(int(df[col].isnull().sum())),
            str(int(df[col].notnull().sum())),
            str(int(df[col].nunique())),  # counts distinct non-null values
            str(int(df.groupby("strain")[col].nunique(dropna=False).max())),
            describe(col),
        )
        for col in df.columns
    ]
    widths = [
        max(len(row[i]) for row in [header, *rows]) for i in range(len(header) - 1)
    ]
    lines.append("")
    for row in [header, *rows]:
        lines.append(
            "  "
            + row[0].ljust(widths[0])
            + "".join(f"  {v.rjust(w)}" for v, w in zip(row[1:-1], widths[1:]))
            + f"  {row[-1]}"
        )

    unlisted = [col for col in df.columns if col not in col_specs]
    lines.append(
        "\ncolumns not listed in `viral_library_validations: columns`: "
        + (str(unlisted) if unlisted else "none")
    )
    for line in lines:
        emit(line)


# width the check descriptions are padded to, so that the counts beside them line up
DESCRIPTION_WIDTH = 76

report_lines = []


def emit(line=""):
    """Print a line of the report, keeping it for the output file as well."""
    print(line)
    report_lines.append(line)


def check(description, examined, offenders):
    """Report a check of `examined`, raising `ValueError` if `offenders` is non-empty."""
    if len(offenders):
        raise ValueError(f"{description}: {len(offenders)} failures\n{offenders}")
    emit(f"  OK  {description.ljust(DESCRIPTION_WIDTH)}  ({examined})")


# the summary tabulates the listed columns, so their presence is settled first; every
# check of the library's contents then comes after the summary
emit("=== schema ===")
check(
    "all columns listed in `columns` are present",
    f"{len(col_specs)} columns",
    [col for col in col_specs if col not in df.columns],
)

emit()
summarize()

emit("\n=== columns ===")
non_null_cols = [col for col, spec in col_specs.items() if not spec["nulls"]]
check(
    "columns declared non-null have no null values",
    f"{len(non_null_cols)} columns",
    {
        col: int(df[col].isnull().sum())
        for col in non_null_cols
        if df[col].isnull().any()
    },
)
for col, spec in col_specs.items():
    if "non_null_only_when" in spec:
        ((when_col, when_val),) = spec["non_null_only_when"].items()
        check(
            f"`{col}` is non-null for exactly the rows with {when_col} == {when_val}",
            f"{len(df)} rows",
            df.loc[df[col].notnull() != (df[when_col] == when_val), "strain"]
            .unique()
            .tolist(),
        )
    if "values" in spec:
        check(
            f"`{col}` holds only {spec['values']}",
            f"{int(df[col].notnull().sum())} non-null rows",
            df.loc[~df[col].isin(spec["values"]) & df[col].notnull(), col]
            .unique()
            .tolist(),
        )

# checked before the uniqueness section below, so that a malformed identifier is
# reported as malformed rather than slipping past a comparison of exact strings
emit("\n=== identifiers ===")
# uppercase is required, not merely conventional: the uniqueness check below compares
# exact strings, so a lowercase barcode would not collide with its uppercase twin
check(
    f"`barcode` is exactly {barcode_length} uppercase ACGT nucleotides",
    f"{len(df)} barcodes",
    df.loc[
        ~df["barcode"].str.fullmatch(f"[ACGT]{{{barcode_length}}}", na=False),
        "barcode",
    ].tolist(),
)
check(
    "`strain` ends in an underscore followed by its `subtype`",
    f"{df['strain'].nunique()} strains",
    [
        row.strain
        for row in df.itertuples()
        if not row.strain.endswith(f"_{row.subtype}")
    ],
)

emit("\n=== uniqueness ===")
for col, spec in col_specs.items():
    if spec.get("unique") == "barcode":
        check(
            f"`{col}` is distinct in every row",
            f"{len(df)} rows",
            df.loc[df[col].duplicated(), col].tolist(),
        )
    elif spec.get("unique") == "strain":
        per_value = df.drop_duplicates(["strain", col]).groupby(col)["strain"]
        check(
            f"each `{col}` value names one strain",
            f"{df[col].nunique()} values",
            per_value.apply(list)[per_value.nunique() > 1].tolist(),
        )
per_strain = df[list(col_specs)].drop(columns=per_barcode_cols).drop_duplicates()
check(
    f"each strain has one set of values for all listed columns but {per_barcode_cols}",
    f"{df['strain'].nunique()} strains",
    per_strain.loc[per_strain["strain"].duplicated(keep=False), "strain"].tolist(),
)

emit("\n=== sequences ===")
# the `subtype` values themselves are checked against the column spec above; this
# catches a subtype whose ectodomain conventions the configuration forgot to give
check(
    "`by_subtype` covers every subtype the `subtype` column may hold",
    f"{len(col_specs['subtype']['values'])} subtypes",
    [s for s in col_specs["subtype"]["values"] if s not in by_subtype],
)
for subtype, spec in by_subtype.items():
    prot = df.loc[
        df["subtype"] == subtype, ["strain", "protein_sequence_HA_ectodomain"]
    ]
    for description, failed in [
        (
            f"{subtype} protein length is one of {spec['protein_lengths']}",
            ~prot["protein_sequence_HA_ectodomain"]
            .str.len()
            .isin(spec["protein_lengths"]),
        ),
        (
            f"{subtype} protein starts with {spec['protein_startswith']}",
            ~prot["protein_sequence_HA_ectodomain"].str.match(
                spec["protein_startswith"]
            ),
        ),
        (
            f"{subtype} protein ends with {spec['protein_endswith']}",
            ~prot["protein_sequence_HA_ectodomain"].str.contains(
                spec["protein_endswith"] + "$", regex=True
            ),
        ),
    ]:
        check(
            description,
            f"{prot['strain'].nunique()} strains",
            prot.loc[failed, "strain"].unique().tolist(),
        )

check(
    "`protein_sequence_HA_ectodomain` has no stop codons",
    f"{df['strain'].nunique()} strains",
    df.loc[df["protein_sequence_HA_ectodomain"].str.contains(r"\*"), "strain"].tolist(),
)
check(
    "`nt_sequence_HA_ectodomain` is uppercase ACGT and a whole number of codons",
    f"{df['strain'].nunique()} strains",
    df.loc[
        ~df["nt_sequence_HA_ectodomain"].str.fullmatch("([ACGT]{3})+", na=False),
        "strain",
    ].tolist(),
)
check(
    "`nt_sequence_HA_ectodomain` translates to `protein_sequence_HA_ectodomain`",
    f"{len(df)} rows",
    [
        row.strain
        for row in df.itertuples()
        if str(Bio.Seq.Seq(row.nt_sequence_HA_ectodomain).translate())
        != row.protein_sequence_HA_ectodomain
    ],
)

emit("\n=== dates ===")
check(
    "`collection_date` is a number",
    f"{df['strain'].nunique()} strains",
    df.loc[
        pd.to_numeric(df["collection_date"], errors="coerce").isna(), "strain"
    ].tolist(),
)
# one date per strain, so a strain carrying several barcodes is not weighted by them;
# two decimals is the precision the dates are recorded to
dates = pd.to_numeric(df.drop_duplicates("strain")["collection_date"])
emit(
    f"      `collection_date` over {len(dates)} strains:  min {dates.min():.2f}"
    f"  median {dates.median():.2f}  max {dates.max():.2f}"
)

with open(snakemake.output.validation, "w") as f:
    f.write("\n".join(report_lines) + "\n")
