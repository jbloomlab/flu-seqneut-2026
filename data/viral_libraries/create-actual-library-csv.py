"""Script to build the actual library CSV from the designed library CSV.

The *designed* library lists every barcode that was designed into the library.
The *actual* library lists only the barcodes used for titer measurements. The
dropped (strain, barcode) pairs are recorded in `drop_strains.csv`, which
also carries a free-text `note`.
"""

import pandas as pd

designed_csv = "flu-seqneut-2026-barcode-to-strain-designed.csv"
drop_csv = "drop_strains.csv"
output_csv = "flu-seqneut-2026-barcode-to-strain-actual.csv"

# columns required in `drop_csv`; `note` is provenance only and is not used for
# matching against the designed library.
drop_columns = ["strain", "barcode", "note"]

# columns matched on when removing dropped rows from the designed library
match_columns = ["strain", "barcode"]


def read_drops(csv_path):
    """Read `drop_csv`, validating its columns and checking for duplicates.

    Fails fast if a required column is missing, if any `strain`/`barcode` value
    is null, or if the same (strain, barcode) pair is listed more than once
    (which would mean the same drop was recorded twice, likely with conflicting
    notes).
    """
    df = pd.read_csv(csv_path)

    missing_columns = [col for col in drop_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(
            f"{csv_path} is missing required column(s): {missing_columns}; "
            f"it must have columns {drop_columns}"
        )
    df = df[drop_columns]

    null_rows = df[match_columns].isnull().any(axis=1)
    if null_rows.any():
        raise ValueError(
            f"{csv_path} has row(s) with null strain or barcode:\n{df[null_rows]}"
        )

    duplicated = df.duplicated(match_columns, keep=False)
    if duplicated.any():
        raise ValueError(
            f"{csv_path} lists the same (strain, barcode) more than once:\n"
            f"{df[duplicated].sort_values(match_columns)}"
        )

    return df


def check_drops_present(designed, drops):
    """Confirm every (strain, barcode) in `drops` is present in `designed`.

    Fails fast on any drop that does not match a designed row, since that means
    the drop list is stale or has a typo and the intended barcode would silently
    remain in the actual library.
    """
    designed_pairs = set(map(tuple, designed[match_columns].values))
    drop_pairs = list(map(tuple, drops[match_columns].values))
    not_found = [pair for pair in drop_pairs if pair not in designed_pairs]
    if not_found:
        raise ValueError(
            f"{len(not_found)} (strain, barcode) pair(s) in {drop_csv} are not in "
            f"{designed_csv}: {not_found}"
        )


def drop_rows(designed, drops):
    """Return `designed` with the `drops` rows removed, matching on strain+barcode.

    Fails fast if the number of removed rows does not equal the number of drops,
    which would mean the designed library has duplicate (strain, barcode) rows.
    """
    to_drop = designed.merge(
        drops[match_columns], on=match_columns, how="left", indicator=True
    )
    keep = (to_drop["_merge"] == "left_only").values

    n_dropped = len(designed) - keep.sum()
    if n_dropped != len(drops):
        raise ValueError(
            f"Dropped {n_dropped} rows from {designed_csv} but {drop_csv} lists "
            f"{len(drops)} drops; the designed library likely has duplicate "
            f"(strain, barcode) rows"
        )

    return designed[keep]


def report_drops(designed, actual, drops):
    """Print a per-strain summary of what was dropped."""
    print(
        f"Read {len(designed)} rows ({designed['strain'].nunique()} strains) "
        f"from {designed_csv}"
    )
    print(f"Read {len(drops)} drops from {drop_csv}")
    print()

    for strain, group in drops.groupby("strain", sort=True):
        n_designed = (designed["strain"] == strain).sum()
        n_actual = (actual["strain"] == strain).sum()
        print(f"=== {strain} ===")
        print(f"barcodes: {n_designed} designed -> {n_actual} actual")
        for _, row in group.iterrows():
            print(f"  dropped {row['barcode']}: {row['note']}")
        print()

    dropped_strains = sorted(set(designed["strain"]) - set(actual["strain"]))
    if dropped_strains:
        print(f"Strains with no remaining barcodes: {dropped_strains}")
        print()


if __name__ == "__main__":
    designed = pd.read_csv(designed_csv)
    drops = read_drops(drop_csv)
    check_drops_present(designed, drops)
    actual = drop_rows(designed, drops)
    report_drops(designed, actual, drops)
    actual.to_csv(output_csv, index=False)
    print(f"Wrote {len(actual)} rows to {output_csv}")
