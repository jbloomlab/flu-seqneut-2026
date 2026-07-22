"""
generate_constructs.py

Processes new orders specified in a config.yaml file and generates the respective order CSV (for 
upload to the GenScript ordering portal) and the plasmid log CSV for recording these new plasmids in
the flu-seqneut-plasmid-log.

Usage:
    python generate_constructs.py \
        --order               H3-testset \
        --config              config.yaml \
        --output-order-csv    results/orders/H3-testset.csv \
        --output-log-csv      results/plasmid_log/2026-05-12_H3-testset.csv \
        --output-excluded-csv results/excluded_from_order/2026-05-12_H3-testset_excluded.csv

Intended for use in a Snakemake pipeline (one invocation per order).
"""

import argparse
import csv
import logging
import random
from pathlib import Path

import yaml

from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

NUCLEOTIDES = list("ACGT")
MAX_BARCODE_ATTEMPTS = 1000
MIN_HAMMING_DISTANCE = 3


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """Load and return the YAML config as a dict."""
    with open(config_path) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Barcode loading
# ---------------------------------------------------------------------------

def _barcodes_from_csv(path: Path) -> set[str]:
    """Return all values in the 'barcode' column of a CSV file."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return {row["barcode"] for row in reader}


def load_past_barcodes(barcode_files: list[str]) -> set[str]:
    """
    Read one or more CSVs, each containing a 'barcode' column, and return
    all barcodes as a set. Using a set allows O(1) lookup when checking for
    duplicates during generation.
    """
    barcodes = set()
    for path in barcode_files:
        barcodes |= _barcodes_from_csv(Path(path))
    log.info(f"Loaded {len(barcodes)} unique past barcodes from {len(barcode_files)} file(s).")
    return barcodes


def load_existing_order_barcodes(log_dir: Path = Path("results/plasmid_log")) -> set[str]:
    """
    Read any existing log CSVs in results/plasmid_log/ and return their barcodes.
    This ensures we don't duplicate barcodes across orders, even across pipeline runs.
    """
    log_files = list(log_dir.glob("*.csv"))
    if not log_files:
        return set()

    barcodes = set()
    for path in log_files:
        barcodes |= _barcodes_from_csv(path)
    log.info(f"Loaded {len(barcodes)} unique barcodes from {len(log_files)} existing log file(s).")
    return barcodes


# ---------------------------------------------------------------------------
# Protein sequence (duplicate-detection) loading
# ---------------------------------------------------------------------------

def _protein_seqs_from_csv(path: Path, id_column: str) -> dict[str, list[str]]:
    """
    Read a CSV with arbitrary columns and return a dict mapping
    protein_sequence_HA_ectodomain -> list of identifier values from `id_column`.

    Rows with an empty protein sequence are skipped. Rows with a missing or empty
    identifier value get an empty-string placeholder so the match is still recorded.
    """
    seq_to_ids: dict[str, list[str]] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or "protein_sequence_HA_ectodomain" not in reader.fieldnames:
            raise ValueError(
                f"{path}: required column 'protein_sequence_HA_ectodomain' not found "
                f"(found: {reader.fieldnames})"
            )
        if id_column not in reader.fieldnames:
            raise ValueError(
                f"{path}: required identifier column {id_column!r} not found "
                f"(found: {reader.fieldnames})"
            )
        for row in reader:
            protein_seq = (row.get("protein_sequence_HA_ectodomain") or "").strip()
            if not protein_seq:
                continue
            identifier = (row.get(id_column) or "").strip()
            seq_to_ids.setdefault(protein_seq, []).append(identifier)
    return seq_to_ids


def _merge_seq_id_maps(
    target: dict[str, list[str]],
    other: dict[str, list[str]],
) -> None:
    """Merge `other` into `target` in place, concatenating identifier lists."""
    for seq, ids in other.items():
        target.setdefault(seq, []).extend(ids)


def load_past_protein_sequences(
    dedicated_files: list[str],
    log_dir: Path = Path("results/plasmid_log"),
) -> dict[str, list[str]]:
    """
    Build a dict mapping known protein_sequence_HA_ectodomain values to a list of
    identifiers (construct_id from dedicated file(s); bloom_lab_plasmid_log_id from
    existing log CSVs) that share that sequence.

    Identifier lists can contain multiple entries when several past constructs map
    to the same protein sequence.
    """
    seq_to_ids: dict[str, list[str]] = {}

    for path in dedicated_files:
        _merge_seq_id_maps(
            seq_to_ids,
            _protein_seqs_from_csv(Path(path), id_column="construct_id"),
        )
    log.info(
        f"Loaded {len(seq_to_ids)} unique protein sequences from "
        f"{len(dedicated_files)} dedicated file(s)."
    )

    log_files = list(log_dir.glob("*.csv"))
    if log_files:
        n_before = len(seq_to_ids)
        for path in log_files:
            _merge_seq_id_maps(
                seq_to_ids,
                _protein_seqs_from_csv(path, id_column="bloom_lab_plasmid_log_id"),
            )
        log.info(
            f"Loaded protein sequences from {len(log_files)} existing log file(s); "
            f"total unique sequences now {len(seq_to_ids)} (added {len(seq_to_ids) - n_before} new)."
        )

    return seq_to_ids


# ---------------------------------------------------------------------------
# Annotation loading
# ---------------------------------------------------------------------------

def load_annotations(annotation_file: str) -> dict[str, dict]:
    """
    Read the annotation CSV and return a dict keyed by protein_sequence_HA_ectodomain.
    Each value contains vaccine_annotation and passage_history_annotation.
    Only rows with a non-empty protein sequence are indexed.
    """
    annotations = {}
    with open(annotation_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            protein_seq = row["protein_sequence_HA_ectodomain"]
            if protein_seq:
                annotations[protein_seq] = {
                    "vaccine_annotation": row["vaccine_annotation"],
                    "passage_history_annotation": row["passage_history_annotation"],
                }
    log.info(f"Loaded annotations for {len(annotations)} sequences from {annotation_file}.")
    return annotations


# ---------------------------------------------------------------------------
# Barcode generation
# ---------------------------------------------------------------------------

def _hamming_distance(a: str, b: str) -> int:
    """Return the Hamming distance between two equal-length strings."""
    if len(a) != len(b):
        raise ValueError("Strings must have equal length")
    return sum(x != y for x, y in zip(a, b))


def _too_close_to_existing(barcode: str, past_barcodes: set[str]) -> bool:
    """
    Return True if barcode is within MIN_HAMMING_DISTANCE - 1 (i.e. <= 2) of
    any barcode already in past_barcodes.

    The exact-match check (Hamming == 0) is assumed to have been done already
    via `barcode in past_barcodes`, so this function only needs to catch
    near-matches (Hamming distance 1 or 2).
    """
    for existing in past_barcodes:
        if _hamming_distance(barcode, existing) < MIN_HAMMING_DISTANCE:
            return True
    return False


def generate_barcodes(n_barcodes: int, past_barcodes: set[str]) -> list[str]:
    """
    Generate n_barcodes unique 16-nt barcodes, excluding any in past_barcodes.

    Barcodes are rejected if they:
      - start with 'GG'
      - already exist in past_barcodes (Hamming distance == 0)
      - are within Hamming distance <= 2 of any barcode in past_barcodes

    The exact-match check is done first (O(1) set lookup) so the more expensive
    Hamming scan is only reached by candidates that would otherwise pass.

    Raises RuntimeError if a barcode cannot be generated after MAX_BARCODE_ATTEMPTS tries.
    """
    barcodes = []
    for _ in range(n_barcodes):
        for attempt in range(MAX_BARCODE_ATTEMPTS):
            barcode = "".join(random.choices(NUCLEOTIDES, k=16))
            if barcode[:2] == "GG":
                continue
            if barcode in past_barcodes:
                continue
            if _too_close_to_existing(barcode, past_barcodes):
                continue
            past_barcodes.add(barcode)  # claim immediately to avoid duplicates within this run
            barcodes.append(barcode)
            break
        else:
            raise RuntimeError(
                f"Failed to generate a unique barcode after {MAX_BARCODE_ATTEMPTS} attempts. "
                "Consider resetting the barcode index."
            )
    log.info(f"Generated {len(barcodes)} unique barcodes.")
    return barcodes


# ---------------------------------------------------------------------------
# Plate-based indexing
# ---------------------------------------------------------------------------

PLATE_ROWS = "ABCDEFGH"
PLATE_COLS = 12
WELLS_PER_PLATE = len(PLATE_ROWS) * PLATE_COLS  # 96


def _parse_well(well: str) -> int:
    """
    Convert a well ID like 'A1' or 'E9' to a 0-based row-major index (0..95).
    Raises ValueError if the well ID is malformed or out of range.
    """
    well = well.upper()
    if len(well) < 2 or well[0] not in PLATE_ROWS:
        raise ValueError(f"Invalid well ID: {well!r}")
    try:
        col = int(well[1:])
    except ValueError as e:
        raise ValueError(f"Invalid well ID: {well!r}") from e
    if not 1 <= col <= PLATE_COLS:
        raise ValueError(f"Well column out of range (1-{PLATE_COLS}): {well!r}")
    row = PLATE_ROWS.index(well[0])
    return row * PLATE_COLS + (col - 1)


def _format_well(well_idx: int) -> str:
    """Convert a 0-based row-major well index (0..95) back to a well ID like 'A1'."""
    row = well_idx // PLATE_COLS
    col = well_idx % PLATE_COLS + 1
    return f"{PLATE_ROWS[row]}{col}"


class PlateIndex:
    """
    Iterator-like helper for plate-format plasmid indices.

    Initialised from a config string like 'plate14' (starts at A1) or
    'plate14-E9' (starts at well E9 of plate 14). Each call to next_idx()
    returns the current 'plateXX-WellID' token and advances by one well,
    rolling over to plate XX+1 well A1 after H12.
    """

    def __init__(self, start: str):
        self.plate_num, self.well_idx = self._parse(start)
        self.start_plate_num = self.plate_num  # for computing per-order relative plate

    @staticmethod
    def _parse(start: str) -> tuple[int, int]:
        # Expected formats: 'plate14' or 'plate14-E9'
        if not start.startswith("plate"):
            raise ValueError(f"Plate index must start with 'plate': {start!r}")
        rest = start[len("plate"):]
        if "-" in rest:
            plate_str, well_str = rest.split("-", 1)
            well_idx = _parse_well(well_str)
        else:
            plate_str, well_idx = rest, 0
        if not plate_str.isdigit():
            raise ValueError(f"Plate number must be an integer: {start!r}")
        return int(plate_str), well_idx

    def next_idx(self) -> tuple[str, int, str]:
        """
        Return (token, relative_plate, well_id) for the current position
        and advance one well.

        - token: 'plateXX-WellID' (absolute plate number, used in filenames)
        - relative_plate: 1-based plate index relative to the order's starting plate
        - well_id: e.g. 'A1', 'E9'
        """
        well_id = _format_well(self.well_idx)
        token = f"plate{self.plate_num}-{well_id}"
        relative_plate = self.plate_num - self.start_plate_num + 1
        self.well_idx += 1
        if self.well_idx >= WELLS_PER_PLATE:
            self.well_idx = 0
            self.plate_num += 1
        return token, relative_plate, well_id


def _is_plate_format(value) -> bool:
    """Return True if value looks like a plate-format index (e.g. 'plate14')."""
    return isinstance(value, str) and value.startswith("plate")


# ---------------------------------------------------------------------------
# Filename generation
# ---------------------------------------------------------------------------

def build_plasmid_filename(idx, vector: str, strain: str, bc_idx: str) -> str:
    """
    Build the GenBank filename for a single construct.
    Format: {idx}_{vector}_{strain}_{bc_idx}.gb
    where strain has '/' replaced with '_'.

    `idx` may be an integer (e.g. 1234) or a plate-format token
    (e.g. 'plate14-E9').
    """
    sanitized_strain = strain.replace("/", "_")
    return f"{idx}_{vector}_{sanitized_strain}_{bc_idx}.gb"


# ---------------------------------------------------------------------------
# Input reading
# ---------------------------------------------------------------------------

def read_input_tsv(input_tsv: str) -> list[dict]:
    """Read the input TSV and return a list of rows as dicts."""
    with open(input_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    log.info(f"Read {len(rows)} strains from {input_tsv}.")
    return rows


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def extract_ectodomain(
    nt_sequence: str,
    strain: str,
    subtype: str,
    start_codon: str,
    special_start_codons: dict[str, str],
    subtype_params: dict,
) -> str:
    """
    Extract the ectodomain insert sequence from a nucleotide sequence.

    Finds the appropriate start codon (strain-specific if provided, otherwise
    the global default), then slices out the ectodomain based on subtype params.
    """
    params = subtype_params[subtype]
    ectodomain_start = params["ectodomain_start"]
    ectodomain_length = params["ectodomain_length"]

    codon_to_use = special_start_codons.get(strain, start_codon)
    start_position = nt_sequence.find(codon_to_use)
    if start_position == -1:
        raise ValueError(f"For {strain} - no start codon '{codon_to_use}' found in sequence")

    insert_start = start_position + ectodomain_start
    insert_end = start_position + ectodomain_length
    return nt_sequence[insert_start:insert_end]


# ---------------------------------------------------------------------------
# Per-row sequence preparation (shared by dedup filter and row building)
# ---------------------------------------------------------------------------

def prepare_input_row(
    row: dict,
    config: dict,
    order_config: dict,
    subtype_params: dict,
    subtype_offsets: dict[str, int],
    start_codon: str,
    special_start_codons: dict[str, str],
    full_sequence: bool,
) -> dict:
    """
    Compute the derived per-strain fields (stripped strain name, normalized
    subtype, ectodomain nt sequence, translated protein sequence) for a single
    input row. Returned dict augments the original row with these fields.

    Doing this once up front lets the duplicate filter and the row builder
    share the same translated protein sequence without re-deriving it.
    """
    strain = _strip_subtype_suffix(row["strain"], row["subtype"])
    subtype = row["subtype"].upper()

    if full_sequence:
        # Full HA CDS provided: locate the start codon and slice between
        # ectodomain_start and ectodomain_length (both measured from the
        # start codon).
        nt_seq = extract_ectodomain(
            nt_sequence=row["nt_sequence"],
            strain=strain,
            subtype=subtype,
            start_codon=start_codon,
            special_start_codons=special_start_codons,
            subtype_params=subtype_params,
        )
        # Per-subtype 5' offset trim happens after extraction.
        offset = subtype_offsets.get(subtype, 0)
        if offset:
            nt_seq = nt_seq[offset:]
    else:
        # Ectodomain region already provided (signal peptide already stripped
        # by the caller, so position 0 of the input is position 0 of the
        # ectodomain).
        #
        # Apply the 5' subtype offset first, then slice from position 0 for
        # a length of (ectodomain_length - ectodomain_start). This produces
        # the same final length as the true branch but treats the input as
        # already-trimmed ectodomain rather than full CDS.
        nt_seq = row["nt_sequence"]
        offset = subtype_offsets.get(subtype, 0)
        if offset:
            nt_seq = nt_seq[offset:]
        ectodomain_start = subtype_params[subtype]["ectodomain_start"]
        ectodomain_length = subtype_params[subtype]["ectodomain_length"]
        slice_length = ectodomain_length - ectodomain_start
        if len(nt_seq) < slice_length:
            raise ValueError(
                f"For {strain} ({subtype}): post-offset nt sequence is "
                f"{len(nt_seq)} nt, shorter than the expected ectodomain "
                f"slice length of {slice_length} nt "
                f"(ectodomain_length={ectodomain_length} - "
                f"ectodomain_start={ectodomain_start}; offset={offset}, "
                f"original input length={len(row['nt_sequence'])} nt). "
                f"Cannot slice the ectodomain."
            )
        nt_seq = nt_seq[:slice_length]

    protein_seq = _translate(nt_seq)

    return {
        **row,
        "_strain": strain,
        "_subtype": subtype,
        "_nt_ectodomain": nt_seq,
        "_protein_ectodomain": protein_seq,
    }


# ---------------------------------------------------------------------------
# Duplicate-sequence filtering
# ---------------------------------------------------------------------------

def filter_duplicate_sequences(
    prepared_rows: list[dict],
    past_protein_sequences: dict[str, list[str]],
) -> tuple[list[dict], list[dict]]:
    """
    Partition prepared rows into (kept, excluded) based on protein sequence dedup.

    A row is excluded if its protein_ectodomain either:
      - appears in `past_protein_sequences` (matched against existing constructs); or
      - has already been kept earlier in this same call (intra-order duplicate).

    For intra-order duplicates, the first occurrence is kept ONLY if it doesn't
    also match an existing past sequence; otherwise all occurrences are excluded.

    Each excluded entry is a dict with:
      - strain, subtype, genbank_accession, derived_haplotype, selection_file,
        protein_sequence_HA_ectodomain
      - reason: 'existing_sequence' or 'intra_order_duplicate'
      - matched_against: ';'-joined identifiers

    `derived_haplotype` and `selection_file` are taken from the input row if
    present; otherwise empty.
    """
    kept: list[dict] = []
    excluded: list[dict] = []
    seen_this_order: dict[str, dict] = {}  # protein_seq -> kept prepared row

    for prepared in prepared_rows:
        protein_seq = prepared["_protein_ectodomain"]
        strain = prepared["_strain"]
        subtype = prepared["_subtype"]
        accession = prepared.get("accession_w_aa_muts_added", "")
        derived_haplotype = prepared.get("derived_haplotype", "") or ""
        selection_file = prepared.get("selection_file", "") or ""

        if protein_seq in past_protein_sequences:
            matched_ids = past_protein_sequences[protein_seq]
            matched_against = ";".join(i for i in matched_ids if i) or "(unknown)"
            log.warning(
                f"Excluding {strain} ({subtype}, {accession}): protein sequence matches "
                f"existing construct(s) [{matched_against}]."
            )
            excluded.append({
                "strain": strain,
                "subtype": subtype,
                "genbank_accession": accession,
                "derived_haplotype": derived_haplotype,
                "selection_file": selection_file,
                "protein_sequence_HA_ectodomain": protein_seq,
                "reason": "existing_sequence",
                "matched_against": matched_against,
            })
            continue

        if protein_seq in seen_this_order:
            first = seen_this_order[protein_seq]
            matched_against = (
                f"{first['_strain']} ({first.get('accession_w_aa_muts_added', '')})"
            )
            log.warning(
                f"Excluding {strain} ({subtype}, {accession}): protein sequence duplicates "
                f"earlier row in this order [{matched_against}]."
            )
            excluded.append({
                "strain": strain,
                "subtype": subtype,
                "genbank_accession": accession,
                "derived_haplotype": derived_haplotype,
                "selection_file": selection_file,
                "protein_sequence_HA_ectodomain": protein_seq,
                "reason": "intra_order_duplicate",
                "matched_against": matched_against,
            })
            continue

        seen_this_order[protein_seq] = prepared
        kept.append(prepared)

    log.info(
        f"Sequence dedup: kept {len(kept)} of {len(prepared_rows)} input rows "
        f"({len(excluded)} excluded)."
    )
    return kept, excluded


def write_excluded_csv(excluded: list[dict], output_path: Path) -> None:
    """Write the excluded-strains CSV (always written, even when empty)."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_path,
        fieldnames=[
            "strain",
            "subtype",
            "genbank_accession",
            "derived_haplotype",
            "selection_file",
            "protein_sequence_HA_ectodomain",
            "reason",
            "matched_against",
        ],
        rows=excluded,
    )
    log.info(f"Excluded-strains CSV written ({len(excluded)} rows): {output_path}")


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def build_construct_rows(
    order: str,
    prepared_rows: list[dict],
    barcodes: list[str],
    n_barcodes: int,
    config: dict,
    annotations: dict[str, dict],
) -> list[dict]:
    """
    Build the enriched per-barcode rows shared by both output CSVs.

    Each prepared input row (one per kept strain, with pre-computed protein /
    ectodomain fields) expands to n_barcodes output rows, one per barcode.
    """
    order_config = config["orders"][order]
    library = order_config["library"]
    plasmid_start_idx = order_config["plasmid_start_idx"]

    # plasmid_start_idx can be an int (sequential numeric indexing) or a
    # plate-format string like 'plate14' / 'plate14-E9' (well-based indexing
    # that rolls over to the next plate after H12).
    if _is_plate_format(plasmid_start_idx):
        plate_iter = PlateIndex(plasmid_start_idx)
        plasmid_idx = None  # unused in plate mode
    else:
        plate_iter = None
        plasmid_idx = plasmid_start_idx

    subtype_params = {k.upper(): v for k, v in config["subtype_params"].items()}

    # Track per-subtype idx counters, initialised from config
    subtype_idx_counters = {
        subtype_key: order_config[f"shortname_{subtype_key.lower()}_idx"]
        for subtype_key in subtype_params
    }

    # Pull subtype-specific sequences once per subtype rather than per strain
    upstream_seq_by_subtype = {
        subtype_key: params["append_additional_upstream_sequence"]
        for subtype_key, params in subtype_params.items()
    }
    endodomain_seq_by_subtype = {
        subtype_key: params["endodomain_sequence"]
        for subtype_key, params in subtype_params.items()
    }

    rows = []
    for i, prepared in enumerate(prepared_rows):
        strain = prepared["_strain"]
        subtype = prepared["_subtype"]
        nt_sequence_HA_ectodomain = prepared["_nt_ectodomain"]
        protein_sequence_HA_ectodomain = prepared["_protein_ectodomain"]
        # derived_haplotype is optional in the input TSV; default to "" if absent.
        derived_haplotype = prepared.get("derived_haplotype", "") or ""

        idx = subtype_idx_counters[subtype]
        subtype_idx_counters[subtype] += 1

        upstream_seq = upstream_seq_by_subtype[subtype]
        endodomain_seq = endodomain_seq_by_subtype[subtype]

        strain_barcodes = barcodes[i * n_barcodes : (i + 1) * n_barcodes]
        for bc_i, barcode in enumerate(strain_barcodes, start=1):
            shortname = f"{library}_{subtype}_{idx}_bc{bc_i}"
            insert_sequence = upstream_seq + nt_sequence_HA_ectodomain + endodomain_seq + barcode
            if plate_iter is not None:
                current_idx, plate_num, well_id = plate_iter.next_idx()
            else:
                current_idx, plate_num, well_id = plasmid_idx, None, None
            plasmid_filename = build_plasmid_filename(
                idx=current_idx,
                vector=config["vector"],
                strain=strain,
                bc_idx=f"bc{bc_i}",
            )
            rows.append({
                "contributor": config["contributor"],
                "vector": config["vector"],
                "shortname": shortname,
                "insert_sequence": insert_sequence,
                "strain": strain,
                "subtype": subtype,
                "barcode": barcode,
                "nt_sequence_HA_ectodomain": nt_sequence_HA_ectodomain,
                "protein_sequence_HA_ectodomain": protein_sequence_HA_ectodomain,
                "genbank_accession": prepared["accession_w_aa_muts_added"],
                "vaccine_annotation": annotations.get(protein_sequence_HA_ectodomain, {}).get("vaccine_annotation", ""),
                "passage_history_annotation": annotations.get(protein_sequence_HA_ectodomain, {}).get("passage_history_annotation", ""),
                "bloom_lab_plasmid_log_id": plasmid_filename,
                "original_publication": "",
                "subclade": "",
                "derived_haplotype": derived_haplotype,
                "library": library,
                "equivalent_strains": "",
                # collection_date is the latest_sequence date propagated from the
                # curated selection TSV via the aggregated library_strains input.
                # Absent for some haplotypes (e.g. older_* additions); "" in that case.
                "collection_date": prepared.get("latest_sequence", "") or "",
                "plate": plate_num,
                "well": well_id,
            })
            if plate_iter is None:
                plasmid_idx += 1
    return rows


def write_order_csv(rows: list[dict], output_path: Path) -> None:
    """Write the order CSV for a single order."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_path,
        fieldnames=["shortname", "insert_sequence"],
        rows=[{"shortname": row["shortname"], "insert_sequence": row["insert_sequence"]} for row in rows],
    )
    log.info(f"Order CSV written: {output_path}")


def write_platemap_csv(rows: list[dict], output_path: Path) -> None:
    """
    Write the platemap CSV for a single order (plate-format orders only).

    Columns: shortname, insert_sequence, plate, well
    `plate` is the 1-based plate index relative to the order's starting plate.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_path,
        fieldnames=["shortname", "insert_sequence", "plate", "well"],
        rows=rows,
    )
    log.info(f"Platemap CSV written: {output_path}")


def write_log_csv(rows: list[dict], output_path: Path) -> None:
    """Write the log CSV for a single order."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_path,
        fieldnames=[
            'contributor', 'vector', 'shortname', 'strain', 'subtype', 'barcode',
            'nt_sequence_HA_ectodomain', 'protein_sequence_HA_ectodomain', 'genbank_accession',
            'vaccine_annotation', 'passage_history_annotation', 'bloom_lab_plasmid_log_id',
            'original_publication', 'subclade', 'derived_haplotype', 'library',
            'equivalent_strains', 'collection_date',
        ],
        rows=rows,
    )
    log.info(f"Log CSV written: {output_path}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _translate(nt_sequence: str) -> str:
    """Translate a nucleotide sequence to amino acids using Biopython."""
    return str(Seq(nt_sequence.upper()).translate())


def _strip_subtype_suffix(strain: str, subtype: str) -> str:
    """Strip _<subtype> suffix from strain name if present (case-insensitive)."""
    suffix = f"_{subtype}"
    if strain.lower().endswith(suffix.lower()):
        return strain[: -len(suffix)]
    return strain


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    """Write a list of dicts to a CSV file."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate CSVs for a single order.")
    parser.add_argument("--order",               required=True, help="Order name.")
    parser.add_argument("--config",              required=True, help="Path to YAML config file.")
    parser.add_argument("--output-order-csv",    required=True, help="Output path for order CSV.")
    parser.add_argument("--output-log-csv",      required=True, help="Output path for log CSV.")
    parser.add_argument("--output-excluded-csv", required=True, help="Output path for excluded CSV.")
    args = parser.parse_args()

    config = load_config(args.config)
    order_config = config["orders"][args.order]

    barcode_files = config["past_barcodes_to_avoid"]
    if isinstance(barcode_files, str):
        barcode_files = [barcode_files]
    random.seed(config["random_seed"])

    past_barcodes = load_past_barcodes(barcode_files)
    past_barcodes |= load_existing_order_barcodes()

    # ------------------------------------------------------------------
    # Protein-sequence dedup: load known sequences and filter input rows
    # before generating any barcodes (so dropped strains don't consume
    # barcodes or plate positions).
    # ------------------------------------------------------------------
    protein_seq_files = config.get("past_protein_sequences_to_avoid", []) or []
    if isinstance(protein_seq_files, str):
        protein_seq_files = [protein_seq_files]
    past_protein_sequences = load_past_protein_sequences(protein_seq_files)

    input_rows = read_input_tsv(order_config["input_file"])

    # Pre-compute derived fields once per input row (shared by dedup + row build)
    subtype_params = {k.upper(): v for k, v in config["subtype_params"].items()}
    subtype_offsets = {
        "H1N1": order_config.get("h1n1_offset", 0) or 0,
        "H3N2": order_config.get("h3n2_offset", 0) or 0,
    }
    start_codon = config["start_codon"]
    special_start_codons = config.get("special_start_codons") or {}
    if not isinstance(special_start_codons, dict):
        special_start_codons = {}
    full_sequence = order_config["full_sequence"]

    prepared_rows = [
        prepare_input_row(
            row=row,
            config=config,
            order_config=order_config,
            subtype_params=subtype_params,
            subtype_offsets=subtype_offsets,
            start_codon=start_codon,
            special_start_codons=special_start_codons,
            full_sequence=full_sequence,
        )
        for row in input_rows
    ]

    kept_rows, excluded_rows = filter_duplicate_sequences(
        prepared_rows=prepared_rows,
        past_protein_sequences=past_protein_sequences,
    )

    # Always write the excluded CSV, even when empty, so the artifact path
    # is predictable for downstream Snakemake rules / audit trails.
    order_csv_path = Path(args.output_order_csv)
    excluded_csv_path = Path(args.output_excluded_csv)
    write_excluded_csv(excluded_rows, excluded_csv_path)

    if not kept_rows:
        log.warning("All input rows were excluded; writing empty order/log CSVs.")

    n_barcodes = config["n_barcodes"]
    barcodes = generate_barcodes(
        n_barcodes=n_barcodes * len(kept_rows),
        past_barcodes=past_barcodes,
    )

    annotations = load_annotations(config["annotation_file"]) if config.get("annotation_file") else {}

    rows = build_construct_rows(
        order=args.order,
        prepared_rows=kept_rows,
        barcodes=barcodes,
        n_barcodes=n_barcodes,
        config=config,
        annotations=annotations,
    )

    write_order_csv(rows=rows, output_path=order_csv_path)
    write_log_csv(rows=rows, output_path=Path(args.output_log_csv))

    if _is_plate_format(order_config["plasmid_start_idx"]):
        platemap_path = order_csv_path.with_name(
            f"{order_csv_path.stem}_platemap{order_csv_path.suffix}"
        )
        write_platemap_csv(rows=rows, output_path=platemap_path)

    log.info("Done.")


if __name__ == "__main__":
    main()