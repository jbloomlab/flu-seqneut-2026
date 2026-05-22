"""
Prepare input files for nextstrain-prot-titers-tree pipeline.

This script combines selected and nonselected haplotypes, validates sequences,
and creates alignment and metadata files for building nextstrain trees.
"""

import pandas as pd
from pathlib import Path
import sys


def validate_sequences(df, seq_col="representative_strain_ha1_sequence"):
    """Validate that all sequences are aligned and contain valid amino acids."""
    sequences = df[seq_col].tolist()

    # Check all sequences exist and are non-empty
    if any(pd.isna(sequences)) or any(len(seq) == 0 for seq in sequences):
        raise ValueError("Some sequences are missing or empty")

    # Check all sequences have the same length (aligned)
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) != 1:
        raise ValueError(
            f"Sequences are not aligned - found lengths: {sorted(set(seq_lengths))}"
        )

    # Check all sequences contain only valid amino acids
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*X-")  # Include stop (*), unknown (X), and gap (-)
    for i, seq in enumerate(sequences):
        invalid_chars = set(seq.upper()) - valid_amino_acids
        if invalid_chars:
            raise ValueError(
                f"Sequence for {df.iloc[i]['representative_strain']} contains "
                f"invalid amino acid characters: {invalid_chars}"
            )

    print(f"Validated {len(sequences)} sequences of length {seq_lengths[0]}")


def create_date_column(df):
    """Create date column from latest_sequence if available, else use 2025.12."""
    if "latest_sequence" in df.columns:
        # Parse dates and convert to decimal year
        df["date"] = pd.to_datetime(df["latest_sequence"]).apply(
            lambda x: x.year + (x.dayofyear - 1) / 365.25
        )
        print(f"Created dates from 'latest_sequence' column (range: {df['date'].min():.2f} - {df['date'].max():.2f})")
    else:
        df["date"] = 2025.12
        print("No 'latest_sequence' column found, using default date of 2025.12")

    return df

def get_ha1_coords(gff_path):
    """Return (start, end) of the HA1 gene from a GFF3 file (1-based, inclusive)."""
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            attrs = dict(kv.split("=") for kv in fields[8].split(";") if "=" in kv)
            attrs = {k: v.strip('"') for k, v in attrs.items()}
            if attrs.get("gene_name") == "HA1":
                return int(fields[3]), int(fields[4])
    raise ValueError('gene_name="HA1" not found in GFF file')

def read_fasta(fasta_path):
    """Return (header, sequence) from a single-record FASTA file."""
    with open(fasta_path) as f:
        header = next(f).strip()
        sequence = "".join(line.strip() for line in f)
    return header, sequence

def extract_ha1_region(fasta_path, ha1_coords, subtype):
    """Return (header, sequence) truncated to ha1_coords (1-based, inclusive, nucleotide).
    Adjusts start by subtype: H1N1 = 6aa before HA1 start, H3N2 = 5aa before."""
    start, end = ha1_coords
    header, sequence = read_fasta(fasta_path)
    aa_start = (start - 1) // 3
    aa_end = end // 3
 
    offsets = {"H1N1": 6, "H3N2": 5}
    if subtype not in offsets:
        raise ValueError(f"Unknown subtype '{subtype}', expected one of {list(offsets)}")
    aa_start -= offsets[subtype]
    aa_end -= offsets[subtype]
 
    return header, sequence[aa_start : aa_end]
 
 
def write_fasta(output_path, header, sequence):
    """Write a sequence to a FASTA file on a single line."""
    with open(output_path, "w") as f:
        f.write(header + "\n")
        f.write(sequence + "\n")

def main(snakemake):
    """Main function called by Snakemake."""

    sys.stderr = sys.stdout = open(snakemake.log[0], "w")

    # Read input files
    print(f"Reading selected haplotypes from {snakemake.input.selected}")
    selected_df = pd.read_csv(snakemake.input.selected, sep="\t")
    assert selected_df["selected_haplotype"].all()

    print(f"Reading nonselected haplotypes from {snakemake.input.nonselected}")
    nonselected_df = pd.read_csv(snakemake.input.nonselected, sep="\t")
    assert not nonselected_df["selected_haplotype"].all()

    # Combine dataframes
    combined_df = pd.concat([selected_df, nonselected_df], ignore_index=True)
    print(f"Combined {len(selected_df)} selected + {len(nonselected_df)} nonselected = {len(combined_df)} total haplotypes")

    # Validate sequences
    validate_sequences(combined_df)

    # Create date column
    combined_df = create_date_column(combined_df)

    # Check that all color_by_metadata columns exist
    color_by_metadata = snakemake.params.color_by_metadata
    missing_columns = []
    for col in color_by_metadata.keys():
        if col not in combined_df.columns:
            missing_columns.append(col)

    if missing_columns:
        raise ValueError(
            f"The following color_by_metadata columns are missing from the data: {missing_columns}"
        )

    print(f"Verified all {len(color_by_metadata)} color_by_metadata columns exist")

    # Write alignment FASTA
    alignment_path = Path(snakemake.output.alignment)
    alignment_path.parent.mkdir(parents=True, exist_ok=True)

    with open(alignment_path, "w") as f:
        for _, row in combined_df.iterrows():
            f.write(f">{row['representative_strain']}\n")
            f.write(f"{row['representative_strain_ha1_sequence']}\n")

    print(f"Wrote alignment to {alignment_path}")

    # Write metadata TSV
    metadata_cols = ["representative_strain"] + ["date"] + list(color_by_metadata.keys())
    metadata_df = combined_df[metadata_cols].copy()
    metadata_df.rename(columns={"representative_strain": "strain"}, inplace=True)

    metadata_path = Path(snakemake.output.metadata)
    metadata_df.to_csv(metadata_path, sep="\t", index=False)
    print(f"Wrote metadata with {len(metadata_df)} strains and {len(metadata_cols)} columns to {metadata_path}")

    # Truncate outgroup sequences to just HA1
    ha1_coords = get_ha1_coords(snakemake.input.gff)
    header, ha1_sequence = extract_ha1_region(snakemake.input.outgroup_fasta, ha1_coords, snakemake.params.subtype)
    write_fasta(snakemake.output.outgroup_ha1, f"{header}_HA1_{ha1_coords[0]}_{ha1_coords[1]}", ha1_sequence)
    print(f"Extracted {len(ha1_sequence)} (positions {ha1_coords[0]}-{ha1_coords[1]}) -> {snakemake.output.outgroup_ha1}")

    print("Done!")


if __name__ == "__main__":
    main(snakemake)