"""
write_genbank.py

Generate plasmid maps for all the plasmids specified in the log CSV file. Intended to be run after
the generate_constructs.py script as part of a Snakemake pipeline. This is specific to cloning into 
the pHW2000 backbone (Bloom lab plasmid #5764). 

Usage:
    python write_genbank.py \
        --log-csv  results/plasmid_log/2026-05-12_H3-testset.csv \
        --outdir   results/genbank/2026-05-12_H3-testset

"""

import argparse
import csv
import logging
import warnings
from pathlib import Path

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

# Fixed sequences within BACKBONE_DOWNSTREAM (sliced by position, verified against the string above).
# These are used to build GenBank feature annotations for the downstream backbone region.
# Layout (0-indexed within BACKBONE_DOWNSTREAM):
#   [  0: 33] Illumina Read1
#   [ 33: 42] gap (gcggccgct)
#   [ 42:147] WSN packaging signal
#   [147:192] 5' NCR  (minus strand; U13 overlaps last 12 nt: [180:192])
#   [192:404] Human PolI promoter  (minus strand)
#   [404:429] gap
#   [429:661] aBGH polyadenylation signal
#   [661:883] gap
#   [883:1444] ori (colE1)
#   [1444:1693] gap
#   [1693:2554] bla (ampR)
#   [2554:2680] gap (trailing backbone)
_BD = BACKBONE_DOWNSTREAM
ILLUMINA_READ1       = _BD[0:33]
WSN_PACKAGING_SIGNAL = _BD[42:147]
NCR_5PRIME           = _BD[147:192]
U13                  = _BD[180:192]
HUMAN_POL1_PROMOTER  = _BD[192:404]
ABGH                 = _BD[429:661]
ORI                  = _BD[883:1444]
BLA                  = _BD[1693:2554]

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
    contributor = row["contributor"]

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
        f"With duplicated 5' packaging signals from WSN with a single stop codon in the duplicated "
        f"packaging signal, and the 16-nucleotide barcode {barcode}. The plasmid was generated for "
        f"the {library} library. It was designed and logged by {contributor} and cloned and "
        f"sequence confirmed by GenScript"
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

    # Build post-endodomain features using a running position cursor.
    # Backbone feature positions are verified against the named sequence constants
    # (ILLUMINA_READ1, WSN_PACKAGING_SIGNAL, etc.) sliced from BACKBONE_DOWNSTREAM.
    pos = f7_end
    downstream_features = []

    def _add(label, length, strand, **kw):
        nonlocal pos
        feat = SeqFeature(
            FeatureLocation(pos, pos + length, strand),
            type="misc_feature",
            qualifiers={"label": label, **kw},
        )
        downstream_features.append(feat)
        pos += length

    # For H3N2, WSN_RECODED_CT follows the endodomain as a separate sequence.
    # For H1N1, it is embedded at the end of WSN_ENDODOMAIN, so annotate it as an
    # overlapping sub-feature spanning the last 33 nt of f7.
    if subtype == "H3N2":
        _add("WSN recoded CT", len(WSN_RECODED_CT), +1)
    else:
        ct_start = f7_end - len(WSN_RECODED_CT)
        downstream_features.append(SeqFeature(
            FeatureLocation(ct_start, f7_end, +1),
            type="misc_feature",
            qualifiers={"label": "WSN recoded CT"},
        ))

    _add("barcode", len(barcode), +1)

    # pos is now at the start of BACKBONE_DOWNSTREAM; verify against the sequence
    bd_start = pos
    assert plasmid_sequence[bd_start:bd_start + len(BACKBONE_DOWNSTREAM)] == BACKBONE_DOWNSTREAM, (
        f"Backbone downstream not found at expected position {bd_start}"
    )

    # Downstream backbone features — positions relative to start of BACKBONE_DOWNSTREAM:
    #   [  0: 33] Illumina Read1         (+1)
    #   [ 42:147] WSN packaging signal   (+1)
    #   [147:192] 5' NCR                 (-1)  (U13 overlaps last 12 nt: [180:192])
    #   [180:192] U13                    (-1)
    #   [192:404] Human PolI promoter    (-1)
    #   [429:661] aBGH                   (+1)
    #   [883:1444] ori (colEI)           (+1)
    #   [1693:2554] bla (ampR)           (+1)
    def _bd(label, seq, strand, **kw):
        offset = BACKBONE_DOWNSTREAM.index(seq)
        feat = SeqFeature(
            FeatureLocation(bd_start + offset, bd_start + offset + len(seq), strand),
            type="misc_feature",
            qualifiers={"label": label, **kw},
        )
        downstream_features.append(feat)

    _bd("Illumina Read1",      ILLUMINA_READ1,       +1)
    _bd("WSN packaging signal", WSN_PACKAGING_SIGNAL, +1)
    _bd("5' NCR",              NCR_5PRIME,            -1)
    _bd("U13",                 U13,                   -1)
    _bd("Human PolI promoter", HUMAN_POL1_PROMOTER,   -1)
    _bd("aBGH",                ABGH,                  +1, note="polyadenylation site")
    _bd("ori",                 ORI,                   +1, note="colEI-origin of replication")
    _bd("bla",                 BLA,                   +1, note="beta-lactamase (ampR)")

    return SeqRecord(
        Seq(plasmid_sequence),
        id=".",
        name=shortname,
        description=definition,
        features=[f1, f2, f3, f4, f5, f6, f7] + downstream_features,
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
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Write GenBank plasmid maps from a log CSV.")
    parser.add_argument("--log-csv", required=True, help="Path to input log CSV.")
    parser.add_argument("--outdir",  required=True, help="Output directory for .gb files.")
    parser.add_argument("--date",    required=True, help="Order date, passed from Snakemake wildcard.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_log_csv(args.log_csv)

    for row in rows:
        filename = row["bloom_lab_plasmid_log_id"]
        record = build_genbank_record(row, date=args.date)
        with open(outdir / filename, "w") as f:
            SeqIO.write(record, f, "genbank")

    log.info(f"Done. Wrote {len(rows)} GenBank files to {outdir}.")


if __name__ == "__main__":
    main()