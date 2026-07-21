"""
aggregate_final_library.py

Aggregation step run after `generate_constructs.py` has produced a plasmid-log
CSV and an excluded CSV for every order. Concatenates all plasmid-log rows
across orders into a single final-library CSV, with two modifications:

  1. The `contributor` column is replaced with `selection_file`, sourced from
     the original input TSV for the order that produced the row (joined on
     strain + protein_sequence_HA_ectodomain).
  2. A new boolean column `need_to_order` is appended:
       - True  for rows coming from the plasmid-log CSVs (newly ordered)
       - False for rows coming from the excluded CSVs (already exist)

Excluded rows are pulled in via `past_protein_sequences_to_avoid` (the same
config entry consumed by generate_constructs.py). For each excluded entry
with reason='existing_sequence', the first up-to-2 IDs from `matched_against`
are looked up in the past-sequences file(s) and emitted as final-library rows.
Intra-order duplicates are skipped.

Usage:
    python aggregate_final_library.py \\
        --config           config.yaml \\
        --plasmid-log-csv  results/plasmid_log/2026-05-12_H3-testset.csv \\
                           results/plasmid_log/2026-05-12_H1-testset.csv \\
        --excluded-csv     results/excluded_from_order/2026-05-12_H3-testset_excluded.csv \\
                           results/excluded_from_order/2026-05-12_H1-testset_excluded.csv \\
        --input-tsv        data/input/H3-testset.tsv \\
                           data/input/H1-testset.tsv \\
        --output           results/final_library/final_library.csv
"""

import argparse
import csv
import datetime
import logging
import re
from pathlib import Path

import yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# Final-library output column order. Mirrors the plasmid-log schema, with
# `contributor` replaced by `selection_file` and `need_to_order` appended.
# `derived_haplotype` is placed first for priority.
FINAL_LIBRARY_FIELDNAMES = [
    "derived_haplotype",
    "selection_file",
    "vector",
    "shortname",
    "strain",
    "subtype",
    "barcode",
    "nt_sequence_HA_ectodomain",
    "protein_sequence_HA_ectodomain",
    "genbank_accession",
    "vaccine_annotation",
    "passage_history_annotation",
    "bloom_lab_plasmid_log_id",
    "original_publication",
    "subclade",
    "library",
    "equivalent_strains",
    "collection_date",
    "need_to_order",
]

SN_ID_PATTERN = re.compile(r"^SN-(\d+)$")


# ---------------------------------------------------------------------------
# Config and CSV loading
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_vaccine_annotations(annotation_file: str) -> dict[str, dict]:
    """
    Load the vaccine annotation file keyed by protein_sequence_HA_ectodomain.

    This is the authoritative source of `vaccine_annotation` /
    `passage_history_annotation`: when a final-library row's HA ectodomain is
    present here, these values override whatever the row inherited from its
    source (plasmid log or the past-sequences construct log). Matching is by HA
    ectodomain sequence, not strain name, mirroring generate_constructs.py.
    """
    annotations: dict[str, dict] = {}
    with open(annotation_file, newline="") as f:
        for row in csv.DictReader(f):
            protein_seq = row["protein_sequence_HA_ectodomain"]
            if protein_seq:
                annotations[protein_seq] = {
                    "vaccine_annotation": row["vaccine_annotation"],
                    "passage_history_annotation": row["passage_history_annotation"],
                }
    log.info(f"Loaded {len(annotations)} vaccine annotations from {annotation_file}.")
    return annotations


def apply_vaccine_annotation_overrides(rows: list[dict], annotations: dict[str, dict]) -> None:
    """
    Override vaccine_annotation / passage_history_annotation on final-library
    rows from the vaccine annotation file (authoritative), matched by HA
    ectodomain sequence. Mutates rows in place; rows whose HA is not in the
    annotation file are left unchanged.
    """
    n = 0
    for row in rows:
        ann = annotations.get(row.get("protein_sequence_HA_ectodomain", ""))
        if ann is not None:
            row["vaccine_annotation"] = ann["vaccine_annotation"]
            row["passage_history_annotation"] = ann["passage_history_annotation"]
            n += 1
    log.info(f"Applied vaccine annotation overrides to {n} final-library rows.")


def _normalize_collection_date(value: str) -> str:
    """Normalize a collection date to YYYY-MM-DD.

    Accepts an already-ISO date (returned unchanged) or a decimal year such as
    "2025.76", which is converted to a calendar date via
    year + round(fraction * days_in_year), with day 1 = Jan 1. Empty stays empty.
    """
    value = (value or "").strip()
    if not value:
        return ""
    if "-" in value:  # already an ISO date
        return value
    dec = float(value)
    year = int(dec)
    start = datetime.date(year, 1, 1)
    days_in_year = (datetime.date(year + 1, 1, 1) - start).days
    return (start + datetime.timedelta(days=round((dec - year) * days_in_year))).isoformat()


# Column names a reference file may use for the collection date, in priority
# order. Different library rounds have named it differently (e.g. `collection_date`
# in 2025to2026, `num_date` in 2025).
_REFERENCE_DATE_COLUMNS = ("collection_date", "num_date")


def load_collection_date_references(reference_files: list[str]) -> dict[str, str]:
    """
    Build a {protein_sequence_HA_ectodomain -> collection_date (ISO)} lookup from
    reference CSVs used to back-fill missing collection dates. Files are consulted
    in order; the first file with a matching sequence wins (later files do not
    override an earlier hit). Rows without a sequence or date are skipped.

    Each file must have a `protein_sequence_HA_ectodomain` column and a date
    column named by one of `_REFERENCE_DATE_COLUMNS`; the date may be ISO or a
    decimal year and is normalized to YYYY-MM-DD.
    """
    lookup: dict[str, str] = {}
    for path in reference_files:
        n = 0
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames or []
            date_col = next((c for c in _REFERENCE_DATE_COLUMNS if c in fieldnames), None)
            if "protein_sequence_HA_ectodomain" not in fieldnames or date_col is None:
                raise ValueError(
                    f"Collection-date reference {path!r} must have a "
                    f"'protein_sequence_HA_ectodomain' column and one of "
                    f"{_REFERENCE_DATE_COLUMNS}. Found columns: {fieldnames}"
                )
            for row in reader:
                seq = (row.get("protein_sequence_HA_ectodomain") or "").strip()
                date_val = (row.get(date_col) or "").strip()
                if seq and date_val and seq not in lookup:
                    lookup[seq] = _normalize_collection_date(date_val)
                    n += 1
        log.info(f"Loaded {n} collection-date references from {path} (date column {date_col!r}).")
    return lookup


def fill_missing_collection_dates(rows: list[dict], references: dict[str, str]) -> None:
    """
    Back-fill an empty `collection_date` on final-library rows from the reference
    lookup, matched by HA ectodomain sequence. Mutates rows in place; only fills
    rows whose collection_date is empty (never overwrites an existing value).
    """
    n = 0
    for row in rows:
        if (row.get("collection_date") or "").strip():
            continue
        date_val = references.get(row.get("protein_sequence_HA_ectodomain", ""))
        if date_val:
            row["collection_date"] = date_val
            n += 1
    log.info(f"Back-filled collection_date for {n} final-library rows.")


def _read_csv_rows(path: Path, delimiter: str = ",") -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def read_plasmid_log_csvs(paths: list[str]) -> list[dict]:
    """Read and concatenate plasmid-log CSVs."""
    rows: list[dict] = []
    for p in paths:
        path_rows = _read_csv_rows(Path(p))
        log.info(f"Read {len(path_rows)} rows from {p}")
        rows.extend(path_rows)
    log.info(f"Total plasmid-log rows: {len(rows)}")
    return rows


def read_excluded_csvs(paths: list[str]) -> list[dict]:
    """Read and concatenate excluded CSVs."""
    rows: list[dict] = []
    for p in paths:
        path_rows = _read_csv_rows(Path(p))
        log.info(f"Read {len(path_rows)} rows from {p}")
        rows.extend(path_rows)
    log.info(f"Total excluded rows: {len(rows)}")
    return rows


# ---------------------------------------------------------------------------
# Input TSV → selection_file lookup
# ---------------------------------------------------------------------------

def build_selection_file_lookup(
    input_tsv_paths: list[str],
) -> tuple[dict[tuple[str, str], str], dict[str, str]]:
    """
    Build a dict keyed by (strain, protein_sequence_HA_ectodomain) → selection_file.

    Used to recover `selection_file` for rows in the plasmid-log CSV (which
    doesn't carry that column). Input TSVs are read raw — no ectodomain
    translation is done here — so the protein-sequence join key is whatever
    is in the TSV's `protein_sequence_HA_ectodomain` column if present.

    NOTE: The plasmid-log CSV stores the *translated* ectodomain. If the input
    TSV doesn't itself carry a protein_sequence_HA_ectodomain column, we fall
    back to keying on strain alone.
    """
    lookup: dict[tuple[str, str], str] = {}
    strain_only_lookup: dict[str, str] = {}
    for tsv in input_tsv_paths:
        path = Path(tsv)
        rows = _read_csv_rows(path, delimiter="\t")
        has_protein = bool(rows) and "protein_sequence_HA_ectodomain" in rows[0]
        for row in rows:
            strain = row.get("strain", "")
            selection_file = (row.get("selection_file") or "").strip()
            if has_protein:
                protein_seq = (row.get("protein_sequence_HA_ectodomain") or "").strip()
                lookup.setdefault((strain, protein_seq), selection_file)
            strain_only_lookup.setdefault(strain, selection_file)
    log.info(
        f"Built selection_file lookup from {len(input_tsv_paths)} TSV(s): "
        f"{len(lookup)} (strain, protein) keys, {len(strain_only_lookup)} strain-only keys."
    )
    return lookup, strain_only_lookup


def lookup_selection_file(
    strain: str,
    protein_seq: str,
    lookup: dict[tuple[str, str], str],
    strain_only_lookup: dict[str, str],
) -> str:
    """Return selection_file for a row, preferring the (strain, protein) key."""
    if (strain, protein_seq) in lookup:
        return lookup[(strain, protein_seq)]
    if strain in strain_only_lookup:
        return strain_only_lookup[strain]
    return ""


# ---------------------------------------------------------------------------
# past_protein_sequences_to_avoid lookup (by SN-XXX id)
# ---------------------------------------------------------------------------

def build_past_sequences_lookup(
    past_files: list[str],
) -> dict[str, dict]:
    """
    Build a dict keyed by construct_id → full row from past_protein_sequences_to_avoid.

    The id column is autodetected per file as either 'construct_id' or
    'bloom_lab_plasmid_log_id' (the same identifiers the excluded CSV records
    in `matched_against`).
    """
    lookup: dict[str, dict] = {}
    for path in past_files:
        rows = _read_csv_rows(Path(path))
        if not rows:
            continue
        fieldnames = list(rows[0].keys())
        if "construct_id" in fieldnames:
            id_col = "construct_id"
        elif "bloom_lab_plasmid_log_id" in fieldnames:
            id_col = "bloom_lab_plasmid_log_id"
        else:
            log.warning(
                f"{path}: no construct_id or bloom_lab_plasmid_log_id column "
                f"found; skipping."
            )
            continue
        for row in rows:
            ident = (row.get(id_col) or "").strip()
            if ident:
                lookup[ident] = row
    log.info(f"Loaded {len(lookup)} past-sequence rows for excluded lookup.")
    return lookup


def _check_sequential_sn_ids(ids: list[str], context: str) -> None:
    """Warn if `ids` aren't sequential SN-\\d+ identifiers."""
    parsed = []
    for ident in ids:
        m = SN_ID_PATTERN.match(ident)
        if not m:
            log.warning(
                f"{context}: matched ID {ident!r} doesn't follow SN-XXX format."
            )
            return
        parsed.append(int(m.group(1)))
    for a, b in zip(parsed, parsed[1:]):
        if b != a + 1:
            log.warning(
                f"{context}: matched IDs {ids} are not sequential "
                f"(found {a} then {b})."
            )
            return


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------

def build_kept_rows(
    plasmid_log_rows: list[dict],
    selection_lookup: dict[tuple[str, str], str],
    strain_only_lookup: dict[str, str],
) -> list[dict]:
    """Map plasmid-log rows to final-library rows (need_to_order=True)."""
    out: list[dict] = []
    missing_selection = 0
    for row in plasmid_log_rows:
        strain = row.get("strain", "")
        protein_seq = row.get("protein_sequence_HA_ectodomain", "")
        selection_file = lookup_selection_file(
            strain, protein_seq, selection_lookup, strain_only_lookup
        )
        if not selection_file:
            missing_selection += 1
        new_row = {col: row.get(col, "") for col in FINAL_LIBRARY_FIELDNAMES}
        new_row["selection_file"] = selection_file
        new_row["need_to_order"] = True
        out.append(new_row)
    if missing_selection:
        log.warning(
            f"{missing_selection} kept rows had no matching selection_file in "
            f"the input TSVs (filled with empty string)."
        )
    return out


def build_excluded_rows(
    excluded_rows: list[dict],
    past_lookup: dict[str, dict],
) -> list[dict]:
    """
    Map excluded CSV rows (reason='existing_sequence' only) to final-library
    rows by looking up the first up-to-2 IDs in `matched_against` against the
    past_protein_sequences_to_avoid lookup. Intra-order and cross-order
    duplicates (same strain + protein + construct_id) are skipped.
    """
    out: list[dict] = []
    seen_construct: set[tuple[str, str, str]] = set()
    n_skipped_intra = 0
    n_skipped_cross = 0
    n_missing_in_past = 0
    for excl in excluded_rows:
        reason = excl.get("reason", "")
        if reason == "intra_order_duplicate":
            n_skipped_intra += 1
            continue
        if reason != "existing_sequence":
            log.warning(
                f"Unknown exclusion reason {reason!r} for strain "
                f"{excl.get('strain', '')!r}; skipping."
            )
            continue

        matched_against = excl.get("matched_against", "") or ""
        ids = [i.strip() for i in matched_against.split(";") if i.strip()]
        ids_to_emit = ids[:2]
        context = (
            f"excluded strain {excl.get('strain', '')!r} "
            f"(matched_against={matched_against!r})"
        )
        _check_sequential_sn_ids(ids_to_emit, context)

        # The selection_file the strain *would have* been ordered for; carried
        # forward so we know which campaign asked for it. derived_haplotype is
        # also useful context — preserve it from the excluded row when present.
        excl_selection_file = (excl.get("selection_file") or "").strip()
        excl_derived_haplotype = (excl.get("derived_haplotype") or "").strip()

        for ident in ids_to_emit:
            past_row = past_lookup.get(ident)
            if past_row is None:
                log.warning(
                    f"{context}: ID {ident!r} not found in "
                    f"past_protein_sequences_to_avoid; skipping this ID."
                )
                n_missing_in_past += 1
                continue
            new_row = {col: past_row.get(col, "") for col in FINAL_LIBRARY_FIELDNAMES}
            # Override columns the past file doesn't supply consistently.
            new_row["selection_file"] = excl_selection_file
            if not new_row.get("derived_haplotype"):
                new_row["derived_haplotype"] = excl_derived_haplotype
            new_row["need_to_order"] = False
            
            # Skip cross-order duplicates: same (strain, protein, construct_id) already emitted.
            key = (
                new_row.get("strain", ""),
                new_row.get("protein_sequence_HA_ectodomain", ""),
                new_row.get("bloom_lab_plasmid_log_id", ""),
            )
            if key in seen_construct:
                n_skipped_cross += 1
                continue
            seen_construct.add(key)
            out.append(new_row)

    log.info(
        f"Excluded → final-library: emitted {len(out)} rows "
        f"(skipped {n_skipped_intra} intra-order duplicates, "
        f"{n_skipped_cross} cross-order duplicates; "
        f"{n_missing_in_past} IDs missing in past-sequences lookup)."
    )
    return out


# ---------------------------------------------------------------------------
# Summary logging
# ---------------------------------------------------------------------------

def log_library_summary(rows: list[dict]) -> None:
    """
    Log a breakdown of the aggregated final-library rows.

    Deduplicates on (strain, protein_sequence_HA_ectodomain) before counting.
    Reports per-subtype totals plus need_to_order split, then within each
    subtype a per-selection_file breakdown with the same split. Empty
    subtype/selection_file values are bucketed as "(none)".
    """
    NONE_BUCKET = "(none)"

    # Deduplicate on (strain, protein_sequence_HA_ectodomain), keeping first.
    # Because main() concatenates kept + excluded, kept (need_to_order=True)
    # wins any collision.
    seen: dict[tuple[str, str], dict] = {}
    n_collisions_need_order = 0
    for row in rows:
        key = (
            row.get("strain", ""),
            row.get("protein_sequence_HA_ectodomain", ""),
        )
        if key in seen:
            if bool(seen[key].get("need_to_order")) != bool(row.get("need_to_order")):
                n_collisions_need_order += 1
            continue
        seen[key] = row
    if n_collisions_need_order:
        log.warning(
            f"{n_collisions_need_order} (strain, protein) key(s) appeared with "
            f"both need_to_order=True and False; kept the first occurrence "
            f"(need_to_order=True wins since kept rows are listed first)."
        )

    unique_rows = list(seen.values())

    # Group by subtype → selection_file → list of rows.
    by_subtype: dict[str, dict[str, list[dict]]] = {}
    for row in unique_rows:
        subtype = (row.get("subtype") or "").strip() or NONE_BUCKET
        selection_file = (row.get("selection_file") or "").strip() or NONE_BUCKET
        by_subtype.setdefault(subtype, {}).setdefault(selection_file, []).append(row)

    def _split(rs: list[dict]) -> tuple[int, int]:
        to_order = sum(1 for r in rs if r.get("need_to_order"))
        return to_order, len(rs) - to_order

    # Subtype overview.
    log.info(
        f"Final library summary: {len(unique_rows)} unique strains "
        f"(by strain + protein_sequence_HA_ectodomain) across "
        f"{len(by_subtype)} subtype(s)."
    )
    log.info("Per-subtype totals (need_to_order True / False):")
    for subtype in sorted(by_subtype):
        all_rows = [r for sf_rows in by_subtype[subtype].values() for r in sf_rows]
        to_order, existing = _split(all_rows)
        log.info(
            f"  {subtype}: {len(all_rows)} total "
            f"({to_order} to order, {existing} already exist)"
        )

    # Per-subtype breakdown by selection_file.
    log.info("Per-subtype breakdown by selection_file:")
    for subtype in sorted(by_subtype):
        log.info(f"  {subtype}:")
        for selection_file in sorted(by_subtype[subtype]):
            sf_rows = by_subtype[subtype][selection_file]
            to_order, existing = _split(sf_rows)
            log.info(
                f"    {selection_file}: {len(sf_rows)} total "
                f"({to_order} to order, {existing} already exist)"
            )


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_final_library_csv(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=FINAL_LIBRARY_FIELDNAMES, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    log.info(f"Final library CSV written ({len(rows)} rows): {output_path}")


def write_haplotypes_by_subtype(rows: list[dict], output_dir: Path) -> None:
    """
    Write derived haplotypes to per-subtype text files, one haplotype per line.
    Files are named {subtype}_final_haplotypes.txt and written to output_dir.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Group unique derived_haplotypes by subtype.
    by_subtype: dict[str, set[str]] = {}
    for row in rows:
        subtype = (row.get("subtype") or "").strip() or "(none)"
        derived_haplotype = (row.get("derived_haplotype") or "").strip()
        if derived_haplotype:  # Only include non-empty haplotypes.
            by_subtype.setdefault(subtype, set()).add(derived_haplotype)
    
    for subtype in sorted(by_subtype):
        filename = f"{subtype}_final_haplotypes.txt"
        filepath = output_dir / filename
        haplotypes = sorted(by_subtype[subtype])
        with open(filepath, "w") as f:
            for haplotype in haplotypes:
                f.write(haplotype + "\n")
        log.info(
            f"Haplotypes file written ({len(haplotypes)} unique haplotypes): {filepath}"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-order plasmid-log and excluded CSVs into a "
                    "single final-library CSV."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file.")
    parser.add_argument(
        "--plasmid-log-csv",
        nargs="+",
        required=True,
        help="One or more per-order plasmid-log CSVs.",
    )
    parser.add_argument(
        "--excluded-csv",
        nargs="+",
        required=True,
        help="One or more per-order excluded CSVs.",
    )
    parser.add_argument(
        "--input-tsv",
        nargs="+",
        required=True,
        help="One or more per-order input TSVs (source of selection_file).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the aggregated final-library CSV.",
    )
    args = parser.parse_args()

    config = load_config(args.config)

    past_files = config.get("past_protein_sequences_to_avoid", []) or []
    if isinstance(past_files, str):
        past_files = [past_files]
    past_lookup = build_past_sequences_lookup(past_files)

    selection_lookup, strain_only_lookup = build_selection_file_lookup(args.input_tsv)

    plasmid_log_rows = read_plasmid_log_csvs(args.plasmid_log_csv)
    excluded_rows_in = read_excluded_csvs(args.excluded_csv)

    kept = build_kept_rows(plasmid_log_rows, selection_lookup, strain_only_lookup)
    excluded = build_excluded_rows(excluded_rows_in, past_lookup)

    final_rows = kept + excluded

    # vaccine_annotations.csv is authoritative for vaccine_annotation /
    # passage_history_annotation: apply it to all final rows (matched by HA
    # ectodomain) so it overrides values inherited from the past-sequences
    # construct log for excluded/pre-existing constructs.
    annotation_file = config.get("annotation_file")
    if annotation_file:
        apply_vaccine_annotation_overrides(final_rows, load_vaccine_annotations(annotation_file))

    # Back-fill missing collection dates from reference CSVs (matched by HA
    # ectodomain sequence), e.g. carrying dates over from a prior library round.
    # Only fills rows whose collection_date is empty.
    reference_files = config.get("collection_date_reference_files", []) or []
    if isinstance(reference_files, str):
        reference_files = [reference_files]
    if reference_files:
        fill_missing_collection_dates(
            final_rows, load_collection_date_references(reference_files)
        )

    output_path = Path(args.output)
    write_final_library_csv(final_rows, output_path)
    write_haplotypes_by_subtype(final_rows, output_path.parent)
    log_library_summary(final_rows)
    log.info("Done.")


if __name__ == "__main__":
    main()