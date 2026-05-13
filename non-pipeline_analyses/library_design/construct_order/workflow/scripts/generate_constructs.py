"""
generate_constructs.py

Processes a single order and writes two output CSVs.

Usage:
    python generate_constructs.py \
        --order            H3-testset \
        --config           config.yaml \
        --output-order-csv results/orders/H3-testset.csv \
        --output-log-csv   results/plasmid_log/2026-05-12_H3-testset.csv

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
MAX_BARCODE_ATTEMPTS = 100


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

def generate_barcodes(n_barcodes: int, past_barcodes: set[str]) -> list[str]:
    """
    Generate n_barcodes unique 16-nt barcodes, excluding any in past_barcodes.

    Barcodes are rejected if they:
      - start with 'GG'
      - already exist in past_barcodes

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
    assert start_position != -1, f"For {strain} - no start codon '{codon_to_use}' found in sequence"

    insert_start = start_position + ectodomain_start
    insert_end = start_position + ectodomain_length
    return nt_sequence[insert_start:insert_end]


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def build_construct_rows(
    order: str,
    input_rows: list[dict],
    barcodes: list[str],
    n_barcodes: int,
    config: dict,
    annotations: dict[str, dict],
) -> list[dict]:
    """
    Build the enriched per-barcode rows shared by both output CSVs.

    Each input row (one per strain) expands to n_barcodes output rows,
    one per barcode, with all derived fields pre-computed.
    """
    order_config = config["orders"][order]
    library = order_config["library"]

    start_codon = config["start_codon"]
    special_start_codons = config.get("special_start_codons")
    if not isinstance(special_start_codons, dict):
        special_start_codons = {}
    subtype_params = config["subtype_params"]

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
    for i, row in enumerate(input_rows):
        strain = _strip_subtype_suffix(row["strain"], row["subtype"])
        subtype = row["subtype"]

        nt_sequence_HA_ectodomain = extract_ectodomain(
            nt_sequence=row["nt_sequence"],
            strain=strain,
            subtype=subtype,
            start_codon=start_codon,
            special_start_codons=special_start_codons,
            subtype_params=subtype_params,
        )
        protein_sequence_HA_ectodomain = _translate(nt_sequence_HA_ectodomain)

        idx = subtype_idx_counters[subtype]
        subtype_idx_counters[subtype] += 1

        upstream_seq = upstream_seq_by_subtype[subtype]
        endodomain_seq = endodomain_seq_by_subtype[subtype]

        strain_barcodes = barcodes[i * n_barcodes : (i + 1) * n_barcodes]
        for bc_i, barcode in enumerate(strain_barcodes, start=1):
            shortname = f"{library}_{subtype}_{idx}_bc{bc_i}"
            insert_sequence = upstream_seq + nt_sequence_HA_ectodomain + endodomain_seq + barcode
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
                "genbank_accession": row["accession_w_aa_muts_added"],
                "vaccine_annotation": annotations.get(protein_sequence_HA_ectodomain, {}).get("vaccine_annotation", ""),
                "passage_history_annotation": annotations.get(protein_sequence_HA_ectodomain, {}).get("passage_history_annotation", ""),
                "bloom_lab_plasmid_log_id": "",
                "original_publication": "",
                "subclade": "",
                "derived_haplotype": "",
                "library": library,
                "equivalent_strains": "",
                "collection_date": "",
            })
    return rows


def write_order_csv(rows: list[dict], output_path: Path) -> None:
    """
    Write the order CSV for a single order.

    TODO: Add construct sequence column once surrounding sequence logic is built out.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_path,
        fieldnames=["shortname", "insert_sequence"],
        rows=[{"shortname": row["shortname"], "insert_sequence": row["insert_sequence"]} for row in rows],
    )
    log.info(f"Order CSV written: {output_path}")


def write_log_csv(rows: list[dict], output_path: Path) -> None:
    """
    Write the log CSV for a single order.

    TODO: Add remaining columns as logic is built out.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_path,
        fieldnames=['contributor', 'vector', 'shortname', 'strain', 'subtype', 'barcode', 'nt_sequence_HA_ectodomain', 'protein_sequence_HA_ectodomain', 'genbank_accession', 'vaccine_annotation', 'passage_history_annotation', 'bloom_lab_plasmid_log_id', 'original_publication', 'subclade', 'derived_haplotype', 'library', 'equivalent_strains', 'collection_date'],
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
    """Strip _<subtype> suffix from strain name if present."""
    suffix = f"_{subtype}"
    if strain.endswith(suffix):
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
    parser.add_argument("--order",            required=True, help="Order name.")
    parser.add_argument("--config",           required=True, help="Path to YAML config file.")
    parser.add_argument("--output-order-csv", required=True, help="Output path for order CSV.")
    parser.add_argument("--output-log-csv",   required=True, help="Output path for log CSV.")
    args = parser.parse_args()

    config = load_config(args.config)
    order_config = config["orders"][args.order]

    barcode_files = config["past_barcodes_to_avoid"]
    if isinstance(barcode_files, str):
        barcode_files = [barcode_files]
    random.seed(config["random_seed"])

    past_barcodes = load_past_barcodes(barcode_files)
    past_barcodes |= load_existing_order_barcodes()

    input_rows = read_input_tsv(order_config["input_file"])

    n_barcodes = config["n_barcodes"]
    barcodes = generate_barcodes(
        n_barcodes=n_barcodes * len(input_rows),
        past_barcodes=past_barcodes,
    )

    annotations = load_annotations(config["vaccine_annotation_file"]) if config.get("vaccine_annotation_file") else {}

    rows = build_construct_rows(
        order=args.order,
        input_rows=input_rows,
        barcodes=barcodes,
        n_barcodes=n_barcodes,
        config=config,
        annotations=annotations,
    )

    write_order_csv(rows=rows, output_path=Path(args.output_order_csv))
    write_log_csv(rows=rows, output_path=Path(args.output_log_csv))

    log.info("Done.")


if __name__ == "__main__":
    main()