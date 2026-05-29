"""
Prepare input files for nextstrain-prot-titers-tree pipeline.

This script combines selected and nonselected haplotypes, validates sequences,
and creates alignment and metadata files for building nextstrain trees.
"""

import pandas as pd
from pathlib import Path
import sys


def validate_sequences(df, seq_col="representative_strain_ha1_ha2_sequence"):
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
    # Intentionally broader than the 20-aa VALID_AA used upstream: outgroup sequences
    # and tree alignments may contain stop codons (*), unknown residues (X), or gaps (-).
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*X-")
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


def get_ha1_ha2_coords(gff_path):
    """
    Return (ha1_start, ha1_end, ha2_start, ha2_end) from a GFF3 file (1-based, inclusive).

    Both HA1 and HA2 features must be present with gene_name attributes.
    """
    ha1_coords = None
    ha2_coords = None

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            attrs = dict(kv.split("=") for kv in fields[8].split(";") if "=" in kv)
            attrs = {k: v.strip('"') for k, v in attrs.items()}
            gene_name = attrs.get("gene_name")
            if gene_name == "HA1":
                ha1_coords = (int(fields[3]), int(fields[4]))
            elif gene_name == "HA2":
                ha2_coords = (int(fields[3]), int(fields[4]))

    if ha1_coords is None:
        raise ValueError('gene_name="HA1" not found in GFF file')
    if ha2_coords is None:
        raise ValueError('gene_name="HA2" not found in GFF file')

    return ha1_coords[0], ha1_coords[1], ha2_coords[0], ha2_coords[1]


def read_fasta(fasta_path):
    """Return (header, sequence) from a single-record FASTA file."""
    with open(fasta_path) as f:
        header = next(f).strip()
        sequence = "".join(line.strip() for line in f)
    return header, sequence


def extract_ha1_ha2_region(fasta_path, ha1_start, ha1_end, ha2_start, ha2_end, subtype):
    """
    Extract and concatenate HA1 and HA2 regions from a full-length FASTA.

    Applies a subtype-specific amino acid offset to the HA1 start. HA2 begins
    immediately where the shifted HA1 ends, so the same offset is applied
    implicitly by concatenating from that point.

    Parameters
    ----------
    fasta_path : str
        Path to full-length FASTA file
    ha1_start : int
        1-based nucleotide start of HA1 in GFF
    ha1_end : int
        1-based nucleotide end of HA1 in GFF
    ha2_start : int
        1-based nucleotide start of HA2 in GFF
    ha2_end : int
        1-based nucleotide end of HA2 in GFF
    subtype : str
        Subtype name used to look up the amino acid offset

    Returns
    -------
    tuple of (str, str, str)
        FASTA header, HA1 sequence, HA2 sequence
    """
    offsets = {"H1N1": 6, "H3N2": 5}
    if subtype not in offsets:
        raise ValueError(f"Unknown subtype '{subtype}', expected one of {list(offsets)}")
    offset = offsets[subtype]

    header, full_sequence = read_fasta(fasta_path)

    # Convert GFF nucleotide coords to 0-indexed amino acid positions
    ha1_aa_start = (ha1_start - 1) // 3 - offset
    ha1_aa_end   = ha1_end // 3 - offset

    ha2_aa_start = (ha2_start - 1) // 3 - offset
    ha2_aa_end   = ha2_end // 3 - offset

    ha1_sequence = full_sequence[ha1_aa_start:ha1_aa_end]
    ha2_sequence = full_sequence[ha2_aa_start:ha2_aa_end]

    return header, ha1_sequence, ha2_sequence


def write_fasta(output_path, header, sequence):
    """Write a sequence to a FASTA file on a single line."""
    with open(output_path, "w") as f:
        f.write(header + "\n")
        f.write(sequence + "\n")


def main(snakemake):
    """Main function called by Snakemake."""

    sys.stdout = sys.stderr = open(snakemake.log[0], "w")

    # Read input files
    print(f"Reading selected haplotypes from {snakemake.input.selected}")
    selected_df = pd.read_csv(snakemake.input.selected, sep="\t")
    if not selected_df["selected_haplotype"].all():
        raise ValueError(f"Selected haplotypes file '{snakemake.input.selected}' contains rows with selected_haplotype=False")

    print(f"Reading nonselected haplotypes from {snakemake.input.nonselected}")
    nonselected_df = pd.read_csv(snakemake.input.nonselected, sep="\t")
    if nonselected_df["selected_haplotype"].any():
        raise ValueError(f"Nonselected haplotypes file '{snakemake.input.nonselected}' contains rows with selected_haplotype=True")

    # Combine dataframes
    combined_df = pd.concat([selected_df, nonselected_df], ignore_index=True)
    print(f"Combined {len(selected_df)} selected + {len(nonselected_df)} nonselected = {len(combined_df)} total haplotypes")

    # Validate combined HA1+HA2 sequences
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

    # Write alignment FASTA using combined HA1+HA2 sequences
    alignment_path = Path(snakemake.output.alignment)
    alignment_path.parent.mkdir(parents=True, exist_ok=True)

    with open(alignment_path, "w") as f:
        for _, row in combined_df.iterrows():
            f.write(f">{row['representative_strain']}\n")
            f.write(f"{row['representative_strain_ha1_ha2_sequence']}\n")

    print(f"Wrote alignment to {alignment_path}")

    # Write metadata TSV
    metadata_cols = ["representative_strain"] + ["date"] + list(color_by_metadata.keys())
    metadata_df = combined_df[metadata_cols].copy()
    metadata_df.rename(columns={"representative_strain": "strain"}, inplace=True)

    metadata_path = Path(snakemake.output.metadata)
    metadata_df.to_csv(metadata_path, sep="\t", index=False)
    print(f"Wrote metadata with {len(metadata_df)} strains and {len(metadata_cols)} columns to {metadata_path}")

    # Extract and concatenate HA1+HA2 from outgroup FASTA
    ha1_start, ha1_end, ha2_start, ha2_end = get_ha1_ha2_coords(snakemake.input.gff)

    header, ha1_sequence, ha2_sequence = extract_ha1_ha2_region(
        snakemake.input.outgroup_fasta,
        ha1_start, ha1_end,
        ha2_start, ha2_end,
        snakemake.params.subtype,
    )

    combined_outgroup_sequence = ha1_sequence + ha2_sequence
    outgroup_header = f"{header}_HA1_HA2_{ha1_start}_{ha1_end}_{ha2_start}_{ha2_end}"
    write_fasta(snakemake.output.outgroup_ha1ha2, outgroup_header, combined_outgroup_sequence)

    print(
        f"Extracted HA1 ({len(ha1_sequence)} aa, nt {ha1_start}-{ha1_end}) + "
        f"HA2 ({len(ha2_sequence)} aa, nt {ha2_start}-{ha2_end}) = "
        f"{len(combined_outgroup_sequence)} aa combined -> {snakemake.output.outgroup_ha1ha2}"
    )

    print("Done!")


if __name__ == "__main__":
    main(snakemake)