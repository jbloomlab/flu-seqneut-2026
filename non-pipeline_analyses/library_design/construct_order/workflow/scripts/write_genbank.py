"""
write_genbank.py

Reads a log CSV and writes one GenBank plasmid map (.gb) per row.
Specific to cloning into the phw2000 backbone (Bloom lab plasmid #5764).

Usage:
    python write_genbank.py \
        --log-csv  results/plasmid_log/2026-05-12_H3-testset.csv \
        --outdir   results/genbank/2026-05-12_H3-testset

Intended for use in a Snakemake pipeline (one invocation per order).
"""

import argparse
import csv
import logging
import warnings
from pathlib import Path

import yaml

from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

warnings.simplefilter("ignore", BiopythonWarning)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Backbone and flanking sequences
# ---------------------------------------------------------------------------

BACKBONE_UPSTREAM = "acctgacgtcgatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatggtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagctctctggctaactagagaacccactgcttactggcttatcgaaattaatacgactcactatagggagacccaagctgttaacgctagcagttaaccggagtactggtcgacctccgaagttgggggggagcaaaagcaggggaaaataaaaacaaccaaa"
WSN_UPSTREAM_H3 = "atgaaggcaaaactactggtcctgttatatgcatttgtagctacagatgcagacaca"
WSN_UPSTREAM_H1 = "atgaaggcaaaactactggtcctgttatatgcatttgtagctacagatgcagacacaata"
BACKBONE_DOWNSTREAM = "agatcggaagagcgtcgtgtagggaaagagtgtgcggccgctatctactcaactgtcgccagttcactggtgctttaggtctccctgggggcaatcagtttctggatgtgttctaatgggtctttgcagtgcagaatatgcatctgagattaggatttcagaaatataaggaaaaacacccttgtttctactaataacccggcggcccaaaatgccgactcggagcgaaagatatacctcccccggggccgggaggtcgcgtcaccgaccacgccgccggcccaggcgacgcgcgacacggacacctgtccccaaaaacgccaccatcgcagccacacacggagcgcccggggccctctggtcaaccccaggacacacgcgggagcagcgccgggccggggacgccctcccggcggtcacctaaatgctagagctcgctgatcagcctcgactgtgccttctagttgccagccatctgttgtttgcccctcccccgtgccttccttgaccctggaaggtgccactcccactgtcctttcctaataaaatgaggaaattgcatcgcattgtctgagtaggtgtcattctattctggggggtggggtggggcaggacagcaagggggaggattgggaagacaatagcaggcatgctggggatgcggtgggctctatggcttctgaggcggaaagaaccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgcc"
CONSENSUS_H3_ENDODOMAIN = "atcaagggagttgagctgaagtcaggatacaaagattggatcctatggatttcctttgccatgtcttgcttcctactgtgcgtagcactactaggctttattatgtgggcgtgtcagaaa"
WSN_ENDODOMAIN = "aaattggaatcaatgggagtgtatcagattctggcgatatattctacagtggcaagctccttagtactgctagtttctttaggagcgattagcttttggatgtgctccaacggctccctacaatgtcggatttgtatttaatag"
WSN_RECODED_CT = "ggctccctacaatgtcggatttgtatttaatag"


# ---------------------------------------------------------------------------
# Input reading
# ---------------------------------------------------------------------------

def read_log_csv(log_csv: str) -> list[dict]:
    """Read the log CSV and return a list of rows as dicts."""
    with open(log_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    log.info(f"Read {len(rows)} rows from {log_csv}.")
    return rows


# ---------------------------------------------------------------------------
# GenBank record construction
# ---------------------------------------------------------------------------

def build_genbank_record(row: dict, date: str) -> SeqRecord:
    """
    Build a BioPython SeqRecord for a single construct row.
    Supports H3N2 and H1N1 subtypes.
    """
    strain = row["strain"]
    subtype = row["subtype"]
    shortname = row["shortname"]
    accession = row["genbank_accession"]
    barcode = row["barcode"]
    ectodomain = row["nt_sequence_HA_ectodomain"]
    library = row["library"]

    if subtype not in {"H3N2", "H1N1"}:
        raise ValueError(f"Unsupported subtype '{subtype}' for strain '{strain}'. Expected 'H3N2' or 'H1N1'.")

    # Build full plasmid sequence
    if subtype == "H3N2":
        plasmid_sequence = BACKBONE_UPSTREAM + WSN_UPSTREAM_H3 + ectodomain + CONSENSUS_H3_ENDODOMAIN + WSN_RECODED_CT + barcode + BACKBONE_DOWNSTREAM
    else:
        plasmid_sequence = BACKBONE_UPSTREAM + WSN_UPSTREAM_H1 + ectodomain + WSN_ENDODOMAIN + barcode + BACKBONE_DOWNSTREAM

    definition = (
        f"This pHW plasmid contains the HA ectodomain sequence for a {subtype} variant {strain}. "
        f"Signal peptide and 3'NCR from WSN, ectodomain from {strain} HA with GenBank accession {accession}, "
        f"and last 46 aa recoded WSN transmembrane and c-terminal domain. "
        f"With duplicated 5' packaging signals from WSN and the 16-nucleotide barcode {barcode}. "
        f"Library: {library}."
    )

    # Build features
    f1 = SeqFeature(FeatureLocation(14, 360, +1),   type="misc_feature", qualifiers={"label": "pCMV", "note": "HCMV-promoter-NdeI deletion"})
    f2 = SeqFeature(FeatureLocation(438, 473, +1),  type="misc_feature", qualifiers={"label": "tI", "note": "terminator for pol1"})
    f3 = SeqFeature(FeatureLocation(473, 505, +1),  type="misc_feature", qualifiers={"label": "3' NCR"})
    f4 = SeqFeature(FeatureLocation(473, 485, +1),  type="misc_feature", qualifiers={"label": "U12"})

    f5_start = 505
    if subtype == "H3N2":
        f5_end = f5_start + 57
        f5 = SeqFeature(FeatureLocation(f5_start, f5_end, +1), type="misc_feature", qualifiers={"label": "WSN first 19 aa"})
    else:
        f5_end = f5_start + 60
        f5 = SeqFeature(FeatureLocation(f5_start, f5_end, +1), type="misc_feature", qualifiers={"label": "WSN first 20 aa"})

    f6_end = f5_end + len(ectodomain)
    f6 = SeqFeature(FeatureLocation(f5_end, f6_end, +1), type="misc_feature", qualifiers={"label": f"HA ectodomain from {strain}"})

    if subtype == "H3N2":
        f7_end = f6_end + len(CONSENSUS_H3_ENDODOMAIN)
        f7 = SeqFeature(FeatureLocation(f6_end, f7_end, +1), type="misc_feature", qualifiers={"label": "consensus H3 endodomain"})
    else:
        f7_end = f6_end + len(WSN_ENDODOMAIN)
        f7 = SeqFeature(FeatureLocation(f6_end, f7_end, +1), type="misc_feature", qualifiers={"label": "WSN endodomain"})

    f8_end  = f7_end  + len(WSN_RECODED_CT)
    f9_end  = f8_end  + 16
    f10_end = f9_end  + 33
    f11_start = f10_end + 9
    f11_end = f11_start + 105
    f12_end = f11_end + 45
    f13_start = f11_end + 33
    f13_end = f12_end
    f14_end = f13_end + 212
    f15_start = f14_end + 25
    f15_end = f15_start + 232
    f16_start = f15_end + 222
    f16_end = f16_start + 561
    f17_start = f16_end + 249
    f17_end = f17_start + 861

    f8  = SeqFeature(FeatureLocation(f7_end,    f8_end,  +1), type="misc_feature", qualifiers={"label": "WSN recoded CT"})
    f9  = SeqFeature(FeatureLocation(f8_end,    f9_end,  +1), type="misc_feature", qualifiers={"label": "barcode"})
    f10 = SeqFeature(FeatureLocation(f9_end,    f10_end, +1), type="misc_feature", qualifiers={"label": "Illumina Read1"})
    f11 = SeqFeature(FeatureLocation(f11_start, f11_end, +1), type="misc_feature", qualifiers={"label": "WSN packaging signal"})
    f12 = SeqFeature(FeatureLocation(f11_end,   f12_end, -1), type="misc_feature", qualifiers={"label": "5' NCR"})
    f13 = SeqFeature(FeatureLocation(f13_start, f13_end, -1), type="misc_feature", qualifiers={"label": "U13"})
    f14 = SeqFeature(FeatureLocation(f13_end,   f14_end, -1), type="misc_feature", qualifiers={"label": "Human PolI promoter"})
    f15 = SeqFeature(FeatureLocation(f15_start, f15_end, +1), type="misc_feature", qualifiers={"label": "aBGH", "note": "polyadenylation site"})
    f16 = SeqFeature(FeatureLocation(f16_start, f16_end, +1), type="misc_feature", qualifiers={"label": "ori", "note": "colEI-origin of replication"})
    f17 = SeqFeature(FeatureLocation(f17_start, f17_end, +1), type="misc_feature", qualifiers={"label": "bla", "note": "beta-lactamase (ampR)"})

    return SeqRecord(
        Seq(plasmid_sequence),
        id=".",
        name=shortname,
        description=definition,
        features=[f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17],
        annotations={
            "source": "synthetic DNA construct",
            "organism": "synthetic DNA construct",
            "molecule_type": "ds-DNA",
            "topology": "circular",
            "data_file_division": "SYN",
            "date": date,
        },
    )


# ---------------------------------------------------------------------------
# Filename generation
# ---------------------------------------------------------------------------

def build_filename(idx: int, vector: str, strain: str, bc_idx: str) -> str:
    """
    Build the GenBank filename for a single construct.
    Format: {idx}_{vector}_{strain}_{bc_idx}.gb
    where strain has '/' replaced with '_'.
    """
    sanitized_strain = strain.replace("/", "_")
    return f"{idx}_{vector}_{sanitized_strain}_{bc_idx}.gb"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Write GenBank plasmid maps from a log CSV.")
    parser.add_argument("--log-csv", required=True, help="Path to input log CSV.")
    parser.add_argument("--outdir",  required=True, help="Output directory for .gb files.")
    parser.add_argument("--date",    required=True, help="Order date, passed from Snakemake wildcard.")
    parser.add_argument("--order",   required=True, help="Order name, passed from Snakemake wildcard.")
    parser.add_argument("--config",  required=True, help="Path to YAML config file.")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    plasmid_idx = config["orders"][args.order]["plasmid_start_idx"]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_log_csv(args.log_csv)

    for row in rows:
        bc_idx = row["shortname"].rsplit("_", 1)[-1]  # extract bc{n} from end of shortname
        filename = build_filename(
            idx=plasmid_idx,
            vector=row["vector"],
            strain=row["strain"],
            bc_idx=bc_idx,
        )
        record = build_genbank_record(row, date=args.date)
        with open(outdir / filename, "w") as f:
            SeqIO.write(record, f, "genbank")
        plasmid_idx += 1

    log.info(f"Done. Wrote {len(rows)} GenBank files to {outdir}.")


if __name__ == "__main__":
    main()