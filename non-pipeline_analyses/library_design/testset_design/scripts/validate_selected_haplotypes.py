"""
Validate HA protein sequences in haplotype data.

This script:
1. Ensures all HA sequences are of equal length
2. Checks which sequences contain only valid amino acids
"""

import sys
import os
import pandas as pd

# Valid amino acid characters (20 standard amino acids)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def has_valid_amino_acids(seq):
    """Check if sequence contains only valid amino acids."""
    return all(aa in VALID_AA for aa in seq)


def validate_sequences(tsv_path, report_path, valid_path, invalid_path, library_fasta_path):
    """
    Validate HA sequences and generate reports.

    Parameters
    ----------
    tsv_path : str
        Path to input TSV file with haplotype data
    report_path : str
        Path to output validation report
    valid_path : str
        Path to output TSV of valid haplotypes
    invalid_path : str
        Path to output TSV of invalid haplotypes
    library_fasta_path : str
        Path to output FASTA file with library (selected) haplotype sequences
    """
    # Read data
    df = pd.read_csv(tsv_path, sep='\t')

    # Get the input TSV filename (basename only)
    selection_file = os.path.basename(tsv_path)

    if 'representative_strain' not in df.columns:
        raise ValueError("Required column 'representative_strain' not found in input data")

    # Check for duplicate values in representative_strain
    duplicates = df['representative_strain'].duplicated(keep=False)
    if duplicates.any():
        dup_values = df[duplicates]['representative_strain'].value_counts()
        examples = []
        for val, count in dup_values.head(5).items():
            dup_rows = df[df['representative_strain'] == val].index.tolist()
            examples.append(f"    '{val}' appears {count} times (rows: {dup_rows[:5]})")
        raise ValueError(
            "Input data contains duplicate values in 'representative_strain':\n" +
            "\n".join(examples) +
            (f"\n    ... and {len(dup_values) - 5} more" if len(dup_values) > 5 else "")
        )

    for col in ('valid_sequence', 'selected_haplotype', 'selection_file'):
        if col in df.columns:
            raise ValueError(
                f"Input data already contains a '{col}' column. "
                "This column will be created by this script and must not already exist."
            )

    # Extract sequences
    sequences = df['representative_strain_ha_sequence']

    # Check sequence lengths
    seq_lengths = sequences.str.len()
    unique_lengths = seq_lengths.unique()

    if len(unique_lengths) != 1:
        raise ValueError(
            f"Sequences have different lengths: {sorted(unique_lengths)}\n"
            f"All sequences must be the same length for proper alignment."
        )

    seq_length = unique_lengths[0]

    # Check amino acid validity
    df['valid_sequence'] = sequences.apply(has_valid_amino_acids)

    # Create selected_haplotype column
    # All specified haplotypes are selected here
    df['selected_haplotype'] = True

    # Add selection_file column with the input TSV filename
    df['selection_file'] = selection_file

    # Validate that all selected haplotypes have valid sequences
    selected_invalid = df[df['selected_haplotype'] & ~df['valid_sequence']]
    if len(selected_invalid) > 0:
        examples = []
        for _, row in selected_invalid.head().iterrows():
            examples.append(
                f"  - {row['representative_strain']}\n"
                f"    Sequence: {row['representative_strain_ha_sequence']}"
            )
        raise ValueError(
            f"Found {len(selected_invalid)} selected haplotype(s) with invalid sequences.\n"
            f"All selected haplotypes must have valid amino acid sequences.\n"
            f"Examples:\n" +
            "\n".join(examples)
        )

    # Count valid vs invalid
    n_valid = df['valid_sequence'].sum()
    n_invalid = (~df['valid_sequence']).sum()
    n_selected = df['selected_haplotype'].sum()

    # Identify invalid characters if any
    invalid_chars = set()
    if n_invalid > 0:
        for seq in sequences[~df['valid_sequence']]:
            invalid_chars.update(set(seq) - VALID_AA)

    # Save valid and invalid haplotypes to separate files
    valid_df = df[df['valid_sequence']]
    invalid_df = df[~df['valid_sequence']]

    valid_df.to_csv(valid_path, sep='\t', index=False)
    invalid_df.to_csv(invalid_path, sep='\t', index=False)

    # Write library FASTA file (selected haplotypes only)
    with open(library_fasta_path, 'w') as f:
        for _, row in valid_df.iterrows():
            f.write(f">{row['representative_strain']}\n")
            f.write(f"{row['representative_strain_ha_sequence']}\n")

    # Write validation report
    with open(report_path, 'w') as f:
        f.write("HA Sequence Validation Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Input file: {tsv_path}\n")
        f.write(f"Selection file (basename): {selection_file}\n\n")

        f.write("Sequence Length Check:\n")
        f.write(f"  All sequences have length: {seq_length}\n")
        f.write(f"  ✓ All sequences are the same length (properly aligned)\n\n")

        f.write("Amino Acid Validity Check:\n")
        f.write(f"  Total haplotypes: {len(df)}\n")
        f.write(f"  Haplotypes with valid amino acids only: {n_valid} ({100*n_valid/len(df):.1f}%)\n")
        f.write(f"  Haplotypes with invalid characters: {n_invalid} ({100*n_invalid/len(df):.1f}%)\n")

        if n_invalid > 0:
            f.write(f"\n  Invalid characters found: {sorted(invalid_chars)}\n")
            f.write("\n  Examples of haplotypes with invalid characters:\n")
            invalid_examples = df[~df['valid_sequence']].head(5)
            for idx, row in invalid_examples.iterrows():
                invalid_in_seq = set(row['representative_strain_ha_sequence']) - VALID_AA
                f.write(f"    {row['representative_strain']}: contains {sorted(invalid_in_seq)}\n")
        else:
            f.write(f"  ✓ All sequences contain only valid amino acids\n")


    print(f"Validation complete:")
    print(f"  Sequence length: {seq_length}")
    print(f"  Valid haplotypes: {n_valid}/{len(df)} ({100*n_valid/len(df):.1f}%)")
    print(f"  Invalid haplotypes: {n_invalid}/{len(df)} ({100*n_invalid/len(df):.1f}%)")
    print(f"  Selected haplotypes: {n_selected}/{len(df)} ({100*n_selected/len(df):.1f}%)")
    print(f"  Report: {report_path}")
    print(f"  Valid haplotypes TSV: {valid_path}")
    print(f"  Invalid haplotypes TSV: {invalid_path}")


if __name__ == '__main__':
    # Snakemake provides these variables
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    validate_sequences(
        tsv_path=snakemake.input.tsv,
        report_path=snakemake.output.report,
        valid_path=snakemake.output.valid_haplotypes,
        invalid_path=snakemake.output.invalid_haplotypes,
        library_fasta_path=snakemake.output.library_fasta_path
    )