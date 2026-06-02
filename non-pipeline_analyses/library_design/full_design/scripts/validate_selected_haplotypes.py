"""
Validate HA protein sequences in haplotype data.

This script:
1. Ensures all HA1 and HA2 sequences are of equal length (checked independently)
2. Checks which sequences contain only valid amino acids (checked independently per segment)
3. Marks a haplotype as valid only if both HA1 and HA2 sequences pass
4. Missing HA2 sequences are treated as invalid
5. Generates a histogram of haplotype counts colored by sequence validity
6. Records the source selection file(s) for each selected haplotype in a `selection_file`
   column (comma-separated basenames if a haplotype was named by multiple selection files).
   Override-excluded rows keep their `selection_file` populated as an audit trail even
   though their `selected_haplotype` flag is set to False.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

# Valid amino acid characters (20 standard amino acids)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def check_sequence_lengths(sequences, col_name):
    """
    Ensure all non-missing sequences in a Series are the same length.

    NaN values are ignored for the length check — missing sequences are
    handled separately as invalid.

    Parameters
    ----------
    sequences : pd.Series
        Series of sequence strings (may contain NaN)
    col_name : str
        Column name, used in error messages

    Returns
    -------
    int
        The uniform sequence length among present sequences,
        or None if all values are missing

    Raises
    ------
    ValueError
        If non-missing sequences have differing lengths
    """
    present = sequences.dropna()
    if len(present) == 0:
        return None

    seq_lengths = present.str.len()
    unique_lengths = seq_lengths.unique()
    if len(unique_lengths) != 1:
        raise ValueError(
            f"Column '{col_name}' has sequences of different lengths: {sorted(unique_lengths)}\n"
            f"All sequences must be the same length for proper alignment."
        )
    return int(unique_lengths[0])


def collect_invalid_chars(sequences, valid_mask):
    """
    Collect the set of invalid characters found in sequences that failed validation.

    Parameters
    ----------
    sequences : pd.Series
        Series of sequence strings
    valid_mask : pd.Series of bool
        Boolean mask where True = valid sequence

    Returns
    -------
    set
        Set of invalid character strings
    """
    invalid_chars = set()
    for seq in sequences[~valid_mask].dropna():
        invalid_chars.update(set(seq) - VALID_AA)
    return invalid_chars


def has_valid_amino_acids(seq):
    """Check if sequence contains only valid amino acids."""
    return all(aa in VALID_AA for aa in seq)


def validate_sequences(tsv_path, report_path, scatter_path, valid_path, invalid_path, select_files, override_haplotypes):
    """
    Validate HA1 and HA2 sequences and generate reports.

    Parameters
    ----------
    tsv_path : str
        Path to input TSV file with haplotype data
    report_path : str
        Path to output validation report
    scatter_path : str
        Path to output scatter plot PDF
    valid_path : str
        Path to output TSV of valid haplotypes
    invalid_path : str
        Path to output TSV of invalid haplotypes
    select_files : list of str
        List of paths to files containing haplotype names to select (one name per line).
        A haplotype is selected if it appears in at least one file.
    override_haplotypes : list
        List of haplotype names (from derived_haplotype column) to exclude from selection
    """
    # Read data
    df = pd.read_csv(tsv_path, sep='\t')

    # Check for duplicate values in key columns
    duplicate_errors = []

    for col_name in ['representative_strain', 'derived_haplotype']:
        if col_name not in df.columns:
            raise ValueError(f"Required column '{col_name}' not found in input data")

        duplicates = df[col_name].duplicated(keep=False)
        if duplicates.any():
            dup_values = df[duplicates][col_name].value_counts()
            examples = []
            for val, count in dup_values.head(5).items():
                dup_rows = df[df[col_name] == val].index.tolist()
                examples.append(f"    '{val}' appears {count} times (rows: {dup_rows[:5]})")

            duplicate_errors.append(
                f"  Column '{col_name}' contains {len(dup_values)} duplicate value(s):\n" +
                "\n".join(examples) +
                (f"\n    ... and {len(dup_values) - 5} more" if len(dup_values) > 5 else "")
            )

    if duplicate_errors:
        raise ValueError(
            "Input data contains duplicate values in columns that must be unique:\n" +
            "\n\n".join(duplicate_errors)
        )

    # Check that output columns don't already exist
    for col_name in ['valid_ha1_sequence', 'valid_ha2_sequence', 'valid_sequence', 'selected_haplotype', 'selection_file']:
        if col_name in df.columns:
            raise ValueError(
                f"Input data already contains a '{col_name}' column. "
                f"This column will be created by this script and must not already exist."
            )

    # Convert count column to numeric, handling comma separators
    if df['count'].dtype == 'object':
        df['count'] = df['count'].str.replace(',', '').astype(int)
    else:
        df['count'] = df['count'].astype(int)

    # Verify all counts are valid
    if df['count'].isna().any():
        raise ValueError(
            f"Some 'count' values could not be converted to numeric.\n"
            f"Rows with invalid counts: {df[df['count'].isna()][['derived_haplotype']].head()}"
        )

    # Read haplotype names from selection files, tracking which file(s) named each haplotype.
    # Build {haplotype_name: [basename, basename, ...]} so we can record the source file(s)
    # in the selection_file column. Order is preserved from the order of select_files.
    print(f"\nReading selection files...")
    selected_haplotype_names = set()
    haplotype_to_files = {}  # derived_haplotype -> list of selection file basenames (in input order)
    existing_haplotypes = set(df['derived_haplotype'])

    for path in select_files:
        basename = os.path.basename(path)
        with open(path) as f:
            names = {line.strip() for line in f if line.strip()}

        if not names:
            raise ValueError(
                f"Selection file '{path}' is empty. "
                f"Each selection file must contain at least one haplotype name."
            )

        print(f"  {path}: {len(names)} haplotype name(s)")

        missing_from_tsv = names - existing_haplotypes
        if missing_from_tsv:
            missing_list = "\n".join(f"    - {h}" for h in sorted(missing_from_tsv))
            print(
                f"  Warning: {len(missing_from_tsv)} haplotype(s) in '{path}' were not found "
                f"in the input TSV and will be ignored:\n{missing_list}"
            )

        selected_haplotype_names.update(names)

        # Record this file as a source for each haplotype it names
        for name in names:
            haplotype_to_files.setdefault(name, []).append(basename)

    print(f"  Total unique haplotype names across all selection files: {len(selected_haplotype_names)}")

    found_names = selected_haplotype_names & existing_haplotypes
    missing_names = selected_haplotype_names - existing_haplotypes
    if missing_names:
        print(
            f"  Warning: {len(missing_names)} total haplotype name(s) from selection files "
            f"were not found in the input TSV and will be ignored."
        )
    print(f"  Haplotype names matched in TSV: {len(found_names)}")

    # Extract sequences
    sequences_ha1 = df['representative_strain_ha1_sequence']
    sequences_ha2 = df['representative_ha2_sequence']

    # Count missing HA2 sequences before any other checks
    n_ha2_missing = sequences_ha2.isna().sum()
    if n_ha2_missing > 0:
        print(f"\nWarning: {n_ha2_missing} haplotype(s) have no HA2 sequence and will be marked invalid.")

    # Check sequence lengths independently per segment (NaNs excluded from length check)
    ha1_seq_length = check_sequence_lengths(sequences_ha1, 'representative_strain_ha1_sequence')
    ha2_seq_length = check_sequence_lengths(sequences_ha2, 'representative_ha2_sequence')

    # Check amino acid validity independently per epitope.
    # Missing (NaN) HA2 sequences are treated as invalid — has_valid_amino_acids
    # is only called on present values; NaN rows are set to False directly.
    df['valid_ha1_sequence'] = sequences_ha1.apply(has_valid_amino_acids)
    df['valid_ha2_sequence'] = sequences_ha2.where(sequences_ha2.notna(), other=None).apply(
        lambda seq: has_valid_amino_acids(seq) if seq is not None else False
    )

    # A haplotype is only valid if both epitopes pass
    df['valid_sequence'] = df['valid_ha1_sequence'] & df['valid_ha2_sequence']

    # Create selected_haplotype column
    df['selected_haplotype'] = df['derived_haplotype'].isin(selected_haplotype_names)

    # Create selection_file column immediately after selected_haplotype.
    # A row's value is a comma-separated string of selection-file basenames that named
    # its derived_haplotype, or pd.NA if no file named it. Note: override exclusions are
    # applied below but do NOT clear this column — selection_file is an audit trail of
    # which files named this haplotype, independent of the final selected_haplotype flag.
    df['selection_file'] = df['derived_haplotype'].map(
        lambda h: ",".join(haplotype_to_files[h]) if h in haplotype_to_files else pd.NA
    )

    # Validate that all override haplotypes exist in the data
    if override_haplotypes:
        missing_overrides = [h for h in override_haplotypes if h not in existing_haplotypes]
        if missing_overrides:
            raise ValueError(
                f"The following override haplotypes were not found in the input data:\n" +
                "\n".join(f"  - {h}" for h in missing_overrides) +
                f"\n\nAvailable haplotypes: {sorted(existing_haplotypes)[:10]}..."
            )

    # Apply override exclusions
    if override_haplotypes:
        initially_selected = df['selected_haplotype'].copy()
        excluded_mask = df['derived_haplotype'].isin(override_haplotypes)
        excluded_haplotypes = df[initially_selected & excluded_mask]

        if len(excluded_haplotypes) > 0:
            print(f"\nExcluding {len(excluded_haplotypes)} haplotype(s) via override_select_recent_haplotypes:")
            for _, row in excluded_haplotypes.iterrows():
                print(f"  - {row['derived_haplotype']} (strain: {row['representative_strain']}, "
                      f"count: {row['count']}, valid_sequence: {row['valid_sequence']}, "
                      f"selection_file: {row['selection_file']})")

        df.loc[excluded_mask, 'selected_haplotype'] = False
        # selection_file is intentionally NOT cleared for override-excluded rows.

    # Validate that all selected haplotypes have valid sequences (both epitopes)
    selected_invalid = df[df['selected_haplotype'] & ~df['valid_sequence']]
    if len(selected_invalid) > 0:
        examples = []
        for _, row in selected_invalid.head().iterrows():
            ha1_ok = row['valid_ha1_sequence']
            ha2_ok = row['valid_ha2_sequence']
            failed = []
            if not ha1_ok:
                failed.append(f"HA1: {row['representative_strain_ha1_sequence']}")
            if not ha2_ok:
                ha2_val = row['representative_ha2_sequence']
                failed.append(f"HA2: {'(missing)' if pd.isna(ha2_val) else ha2_val}")
            examples.append(
                f"  - {row['derived_haplotype']} (strain: {row['representative_strain']})\n"
                + "\n".join(f"    {f}" for f in failed)
            )
        raise ValueError(
            f"Found {len(selected_invalid)} selected haplotype(s) with invalid sequences.\n"
            f"All selected haplotypes must have valid amino acid sequences in both HA1 and HA2.\n"
            f"Examples:\n" +
            "\n".join(examples)
        )

    # Reorder columns so selection_file sits immediately after selected_haplotype.
    # Both columns were created by this script and are appended at the end by default;
    # we move selection_file to right after selected_haplotype for readability.
    cols = list(df.columns)
    cols.remove('selection_file')
    insert_at = cols.index('selected_haplotype') + 1
    cols.insert(insert_at, 'selection_file')
    df = df[cols]

    # Compute per-epitope validity counts
    n_ha1_valid = df['valid_ha1_sequence'].sum()
    n_ha1_invalid = (~df['valid_ha1_sequence']).sum()
    n_ha2_valid = df['valid_ha2_sequence'].sum()
    n_ha2_invalid = (~df['valid_ha2_sequence']).sum()

    # Both valid / one invalid / both invalid breakdowns
    both_valid_mask = df['valid_sequence']
    only_ha1_invalid_mask = ~df['valid_ha1_sequence'] & df['valid_ha2_sequence']
    only_ha2_invalid_mask = df['valid_ha1_sequence'] & ~df['valid_ha2_sequence']
    both_invalid_mask = ~df['valid_ha1_sequence'] & ~df['valid_ha2_sequence']

    n_both_valid = both_valid_mask.sum()
    n_only_ha1_invalid = only_ha1_invalid_mask.sum()
    n_only_ha2_invalid = only_ha2_invalid_mask.sum()
    n_both_invalid = both_invalid_mask.sum()
    n_any_invalid = (~both_valid_mask).sum()
    n_selected = df['selected_haplotype'].sum()

    # Sequence count totals
    total_sequences_valid = df[both_valid_mask]['count'].sum()
    total_sequences_invalid = df[~both_valid_mask]['count'].sum()

    # Collect invalid characters per epitope (NaNs skipped in collect_invalid_chars)
    ha1_invalid_chars = collect_invalid_chars(sequences_ha1, df['valid_ha1_sequence'])
    ha2_invalid_chars = collect_invalid_chars(sequences_ha2, df['valid_ha2_sequence'])

    # Save valid and invalid haplotypes to separate files
    df[both_valid_mask].to_csv(valid_path, sep='\t', index=False)
    df[~both_valid_mask].to_csv(invalid_path, sep='\t', index=False)

    # Write validation report
    with open(report_path, 'w') as f:
        f.write("HA Sequence Validation Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Input file: {tsv_path}\n\n")

        f.write("Sequence Length Check:\n")
        f.write(f"  HA1 sequences: all have length {ha1_seq_length} ✓\n")
        if ha2_seq_length is not None:
            f.write(f"  HA2 sequences: all present sequences have length {ha2_seq_length} ✓\n")
            if n_ha2_missing > 0:
                f.write(f"  HA2 sequences: {n_ha2_missing} missing (treated as invalid)\n")
        else:
            f.write(f"  HA2 sequences: all {n_ha2_missing} values are missing\n")
        f.write("\n")

        f.write("Amino Acid Validity Check:\n")
        f.write(f"  Total haplotypes: {len(df)}\n\n")

        f.write(f"  HA1 (representative_strain_ha1_sequence):\n")
        f.write(f"    Valid:   {n_ha1_valid} ({100*n_ha1_valid/len(df):.1f}%)\n")
        f.write(f"    Invalid: {n_ha1_invalid} ({100*n_ha1_invalid/len(df):.1f}%)\n")
        if ha1_invalid_chars:
            f.write(f"    Invalid characters: {sorted(ha1_invalid_chars)}\n")
            f.write(f"    Examples:\n")
            for _, row in df[~df['valid_ha1_sequence']].head(5).iterrows():
                bad = sorted(set(row['representative_strain_ha1_sequence']) - VALID_AA)
                f.write(f"      {row['derived_haplotype']}: contains {bad}\n")
        else:
            f.write(f"    ✓ All HA1 sequences contain only valid amino acids\n")

        f.write(f"\n  HA2 (representative_ha2_sequence):\n")
        f.write(f"    Valid:   {n_ha2_valid} ({100*n_ha2_valid/len(df):.1f}%)\n")
        f.write(f"    Invalid: {n_ha2_invalid} ({100*n_ha2_invalid/len(df):.1f}%)")
        if n_ha2_missing > 0:
            f.write(f" (includes {n_ha2_missing} missing)\n")
        else:
            f.write("\n")
        if ha2_invalid_chars:
            f.write(f"    Invalid characters: {sorted(ha2_invalid_chars)}\n")
            f.write(f"    Examples:\n")
            for _, row in df[~df['valid_ha2_sequence'] & sequences_ha2.notna()].head(5).iterrows():
                bad = sorted(set(row['representative_ha2_sequence']) - VALID_AA)
                f.write(f"      {row['derived_haplotype']}: contains {bad}\n")
        elif n_ha2_missing == 0:
            f.write(f"    ✓ All HA2 sequences contain only valid amino acids\n")

        if n_ha2_missing > 0:
            f.write(f"    Haplotypes with missing HA2 sequence:\n")
            for _, row in df[sequences_ha2.isna()].head(5).iterrows():
                f.write(f"      {row['derived_haplotype']} (strain: {row['representative_strain']})\n")
            if n_ha2_missing > 5:
                f.write(f"      ... and {n_ha2_missing - 5} more\n")

        f.write(f"\n  Combined validity (both segments must pass):\n")
        f.write(f"    Both valid:       {n_both_valid} ({100*n_both_valid/len(df):.1f}%)\n")
        f.write(f"    HA1 invalid only: {n_only_ha1_invalid} ({100*n_only_ha1_invalid/len(df):.1f}%)\n")
        f.write(f"    HA2 invalid only: {n_only_ha2_invalid} ({100*n_only_ha2_invalid/len(df):.1f}%)\n")
        f.write(f"    Both invalid:     {n_both_invalid} ({100*n_both_invalid/len(df):.1f}%)\n")
        f.write(f"\n  Total sequences (valid haplotypes):   {total_sequences_valid:,}\n")
        f.write(f"  Total sequences (invalid haplotypes): {total_sequences_invalid:,}\n")
        f.write(f"  Total sequences (all):                {total_sequences_valid + total_sequences_invalid:,}\n")

        f.write("\nHaplotype Selection:\n")
        f.write(f"  Selection files: {select_files}\n")
        f.write(f"  Total unique haplotype names in selection files: {len(selected_haplotype_names)}\n")
        if missing_names:
            f.write(f"  Warning: {len(missing_names)} name(s) in selection files not found in TSV:\n")
            for name in sorted(missing_names):
                f.write(f"    - {name}\n")
        if override_haplotypes:
            f.write(f"  Override exclusions: {len(override_haplotypes)} haplotype(s) specified\n")
            excluded_mask = df['derived_haplotype'].isin(override_haplotypes)
            n_excluded = excluded_mask.sum()
            if n_excluded > 0:
                f.write(f"  Excluded via override: {n_excluded} haplotype(s)\n")
                for haplotype in override_haplotypes:
                    if haplotype in df['derived_haplotype'].values:
                        row = df[df['derived_haplotype'] == haplotype].iloc[0]
                        f.write(
                            f"    - {haplotype} (strain: {row['representative_strain']}, "
                            f"count: {row['count']}, selection_file: {row['selection_file']})\n"
                        )
        f.write(f"  Selected haplotypes: {n_selected} ({100*n_selected/len(df):.1f}%)\n")
        total_sequences_selected = df[df['selected_haplotype']]['count'].sum()
        f.write(f"  Total sequences (selected haplotypes): {total_sequences_selected:,}\n")
        f.write(f"  ✓ All selected haplotypes have valid HA1 and HA2 sequences\n")

        # selection_file summary
        f.write(f"\n  selection_file value counts (non-null):\n")
        sf_counts = df['selection_file'].dropna().value_counts()
        if len(sf_counts) == 0:
            f.write(f"    (no haplotypes were named by any selection file)\n")
        else:
            for val, count in sf_counts.items():
                f.write(f"    {val}: {count}\n")

        f.write("\nSequence Count Distribution:\n")
        f.write(f"  Total sequences represented: {df['count'].sum()}\n")
        f.write(f"  Mean sequences per haplotype: {df['count'].mean():.1f}\n")
        f.write(f"  Median sequences per haplotype: {df['count'].median():.0f}\n")
        f.write(f"  Max sequences for a haplotype: {df['count'].max()}\n")

    # Create scatter plot
    # Four validity categories: both valid, HA1 invalid only, HA2 invalid only, both invalid
    fig, ax = plt.subplots(figsize=(12, 8))

    plot_groups = [
        (both_valid_mask,      'steelblue',    'darkblue',  'o',  60,  f'Both valid (n={n_both_valid}, seqs={total_sequences_valid:,})'),
        (only_ha1_invalid_mask,'coral',        'darkred',   '^', 100,  f'HA1 invalid only (n={n_only_ha1_invalid})'),
        (only_ha2_invalid_mask,'mediumpurple', 'indigo',    's', 100,  f'HA2 invalid only / missing (n={n_only_ha2_invalid})'),
        (both_invalid_mask,    'gold',         'goldenrod', 'D', 100,  f'Both invalid (n={n_both_invalid})'),
    ]

    for mask, color, edge_color, marker, size, label in plot_groups:
        group_counts = df[mask].groupby('count').size()
        if len(group_counts) > 0:
            ax.scatter(group_counts.index, group_counts.values,
                       alpha=0.7, s=size, marker=marker,
                       color=color, edgecolors=edge_color, linewidth=0.8,
                       label=label, zorder=2)

    ax.set_xlabel('Number of sequences with haplotype (count)', fontsize=12)
    ax.set_ylabel('Number of haplotypes', fontsize=12)
    ax.set_title(f'Distribution of Haplotype Counts by Sequence Validity\n{tsv_path.split("/")[-1]}', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(axis='both', alpha=0.3, linestyle='--', zorder=1)

    if df['count'].max() > 100:
        ax.set_xscale('log')
        ax.set_xlabel('Number of sequences with haplotype (log scale)', fontsize=12)

    plt.tight_layout()
    plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
    plt.close()

    total_sequences_selected = df[df['selected_haplotype']]['count'].sum()
    print(f"Validation complete:")
    print(f"  HA1 sequence length: {ha1_seq_length}")
    print(f"  HA2 sequence length: {ha2_seq_length} ({n_ha2_missing} missing)")
    print(f"  HA1 valid: {n_ha1_valid}/{len(df)} ({100*n_ha1_valid/len(df):.1f}%)")
    print(f"  HA2 valid: {n_ha2_valid}/{len(df)} ({100*n_ha2_valid/len(df):.1f}%)")
    print(f"  Both valid: {n_both_valid}/{len(df)} ({100*n_both_valid/len(df):.1f}%)")
    print(f"  Any invalid: {n_any_invalid}/{len(df)} ({100*n_any_invalid/len(df):.1f}%)")
    print(f"  Selected haplotypes: {n_selected}/{len(df)} ({100*n_selected/len(df):.1f}%)")
    print(f"  Total sequences (valid): {total_sequences_valid:,}")
    print(f"  Total sequences (selected): {total_sequences_selected:,}")
    print(f"  selection_file value counts (non-null):")
    sf_counts = df['selection_file'].dropna().value_counts()
    if len(sf_counts) == 0:
        print(f"    (no haplotypes were named by any selection file)")
    else:
        for val, count in sf_counts.items():
            print(f"    {val}: {count}")
    print(f"  Report: {report_path}")
    print(f"  Scatter plot: {scatter_path}")
    print(f"  Valid haplotypes TSV: {valid_path}")
    print(f"  Invalid haplotypes TSV: {invalid_path}")


if __name__ == '__main__':
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    validate_sequences(
        tsv_path=snakemake.input.tsv,
        report_path=snakemake.output.report,
        scatter_path=snakemake.output.scatter,
        valid_path=snakemake.output.valid_haplotypes,
        invalid_path=snakemake.output.invalid_haplotypes,
        select_files=snakemake.input.select_files,
        override_haplotypes=snakemake.params.override_haplotypes
    )