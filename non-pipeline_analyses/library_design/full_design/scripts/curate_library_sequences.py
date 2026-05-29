"""
Curate library sequences by combining validated recent haplotypes and additional haplotypes.

This script:
1. Reads validated recent haplotypes (with selected_haplotype already defined)
2. Reads any additional haplotypes in one of two formats:
   a) Split format: HA1 and HA2 columns already separated
   b) Combined format: Full HA amino acid sequence that needs to be split by subtype
3. For combined format, parses GFF (with amino acid coordinates) to extract HA1/HA2 regions
   and slices the amino acid sequence accordingly
4. Combines into a single dataframe. The `selection_file` column is preserved from the
   upstream validation script for validated rows; additional-haplotype rows are tagged
   with the basename of their source file in the `selection_file` column.
5. Creates combined HA1+HA2 sequence column early so it flows through all downstream steps
6. Validates that representative_strain is unique
7. Validates that all HA1, HA2, and combined sequences are the same length and contain only valid amino acids
8. Computes Hamming distance to nearest library strain independently for HA1 and HA1+HA2
9. Computes epitope-specific distances to nearest HA1 library strain for each HA1 epitope set
10. Outputs two TSVs: one for selected haplotypes (sorted by HA1+HA2 distance ascending,
    includes `selection_file` column at the end), one for non-selected haplotypes (sorted by
    HA1+HA2 distance descending, does NOT include `selection_file` column)
11. Outputs a FASTA file with combined HA1+HA2 sequences for selected haplotypes
"""
import sys
import os
import pandas as pd
import re

# Valid amino acid characters (20 standard amino acids)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')


def has_valid_amino_acids(seq):
    """Check if sequence contains only valid amino acids."""
    return all(aa in VALID_AA for aa in seq)


def parse_gff_for_regions(gff_path, gene_names=['HA1', 'HA2']):
    """
    Parse GFF file and extract nucleotide coordinates for specified genes.
    
    Parameters
    ----------
    gff_path : str
        Path to GFF file (with 1-indexed nucleotide coordinates)
    gene_names : list of str
        Gene names to extract (default: ['HA1', 'HA2'])
    
    Returns
    -------
    dict
        Mapping of gene_name → (start_nt, end_nt) in 1-indexed GFF format
    
    Raises
    ------
    ValueError
        If required genes are not found in GFF
    """
    regions = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('##') or line.startswith('#'):
                continue
            if not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
            
            try:
                start_nt = int(fields[3])
                end_nt = int(fields[4])
            except ValueError:
                continue
            
            attributes = fields[8]
            
            # Parse gene_name from attributes (format: gene_name="HA1")
            match = re.search(r'gene_name="([^"]+)"', attributes)
            if match:
                gene_name = match.group(1)
                if gene_name in gene_names:
                    regions[gene_name] = (start_nt, end_nt)
    
    # Validate all required genes were found
    missing_genes = set(gene_names) - set(regions.keys())
    if missing_genes:
        raise ValueError(
            f"GFF file {gff_path} does not contain genes: {missing_genes}\n"
            f"Found genes: {list(regions.keys())}"
        )
    
    return regions


def extract_ha_regions(combined_seq, ha1_coords_nt, ha2_coords_nt, subtype):
    """
    Extract HA1 and HA2 amino acid sequences from combined amino acid sequence.
    
    Converts GFF nucleotide coordinates to amino acid coordinates, applying a
    subtype-specific offset to account for signal peptide removal.
    
    Parameters
    ----------
    combined_seq : str
        Combined HA amino acid sequence
    ha1_coords_nt : tuple of (int, int)
        1-indexed nucleotide coordinates (start, end) for HA1 from GFF
    ha2_coords_nt : tuple of (int, int)
        1-indexed nucleotide coordinates (start, end) for HA2 from GFF
    subtype : str
        Flu subtype ('H1N1' or 'H3N2') to determine amino acid offset
    
    Returns
    -------
    tuple of (str, str)
        HA1 and HA2 amino acid sequences
    
    Raises
    ------
    ValueError
        If sequence length is incompatible with coordinates or unknown subtype
    """
    # Subtype-specific amino acid offsets (account for signal peptide removal)
    offsets = {"H1N1": 6, "H3N2": 5}
    if subtype not in offsets:
        raise ValueError(f"Unknown subtype '{subtype}', expected one of {list(offsets.keys())}")
    offset = offsets[subtype]
    
    # Convert GFF nucleotide coords (1-indexed) to amino acid positions (0-indexed)
    ha1_start, ha1_end = ha1_coords_nt
    ha2_start, ha2_end = ha2_coords_nt
    
    ha1_aa_start = (ha1_start - 1) // 3 - offset
    ha1_aa_end = ha1_end // 3 - offset
    
    ha2_aa_start = (ha2_start - 1) // 3 - offset
    ha2_aa_end = ha2_end // 3 - offset
    
    # Validate sequence is long enough
    required_length = max(ha1_aa_end, ha2_aa_end)
    if len(combined_seq) < required_length:
        raise ValueError(
            f"Combined sequence length ({len(combined_seq)}) is shorter than "
            f"required ({required_length}) based on GFF coordinates and subtype offset"
        )
    
    ha1_aa = combined_seq[ha1_aa_start:ha1_aa_end]
    ha2_aa = combined_seq[ha2_aa_start:ha2_aa_end]
    
    return ha1_aa, ha2_aa


def detect_input_format(df):
    """
    Detect which format the input dataframe is in.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe to check
    
    Returns
    -------
    str
        Either 'split' (has representative_strain_ha1_sequence and representative_ha2_sequence)
        or 'combined' (has representative_strain_ha_sequence)
    
    Raises
    ------
    ValueError
        If neither format is detected
    """
    has_split = (
        'representative_strain_ha1_sequence' in df.columns and
        'representative_ha2_sequence' in df.columns
    )
    has_combined = 'representative_strain_ha_sequence' in df.columns
    
    if has_split and has_combined:
        raise ValueError(
            "Input file has both split (HA1/HA2) and combined (HA) columns. "
            "Please provide only one format."
        )
    
    if has_split:
        return 'split'
    elif has_combined:
        return 'combined'
    else:
        raise ValueError(
            "Input file has neither split format (representative_strain_ha1_sequence + "
            "representative_ha2_sequence) nor combined format (representative_strain_ha_sequence). "
            f"Available columns: {list(df.columns)}"
        )


def process_additional_haplotypes(df_additional, gff_path, input_format, subtype):
    """
    Process additional haplotypes, converting combined format to split if needed.
    
    Parameters
    ----------
    df_additional : pd.DataFrame
        Additional haplotypes dataframe
    gff_path : str
        Path to GFF file (required for combined format)
    input_format : str
        Either 'split' or 'combined'
    subtype : str
        Flu subtype ('H1N1' or 'H3N2') for coordinate conversion
    
    Returns
    -------
    pd.DataFrame
        Dataframe with standard columns (representative_strain_ha1_sequence, representative_ha2_sequence)
    """
    df = df_additional.copy()
    
    if input_format == 'split':
        # Already in the right format, nothing to do
        return df
    
    elif input_format == 'combined':
        # Parse GFF to get HA1/HA2 nucleotide coordinates
        regions = parse_gff_for_regions(gff_path, gene_names=['HA1', 'HA2'])
        ha1_coords = regions['HA1']
        ha2_coords = regions['HA2']
        
        # Extract HA1 and HA2 for each row
        ha1_sequences = []
        ha2_sequences = []
        
        for idx, row in df.iterrows():
            combined_seq = row['representative_strain_ha_sequence']
            try:
                ha1_aa, ha2_aa = extract_ha_regions(combined_seq, ha1_coords, ha2_coords, subtype)
                ha1_sequences.append(ha1_aa)
                ha2_sequences.append(ha2_aa)
            except Exception as e:
                raise ValueError(
                    f"Failed to extract HA regions for {row['representative_strain']}: {e}"
                )
        
        df['representative_strain_ha1_sequence'] = ha1_sequences
        df['representative_ha2_sequence'] = ha2_sequences
        
        # Drop the combined column as it's no longer needed
        df = df.drop(columns=['representative_strain_ha_sequence'])
        
        return df
    
    else:
        raise ValueError(f"Unknown input format: {input_format}")


def hamming_distance(seq1, seq2):
    """
    Calculate Hamming distance between two sequences of equal length.

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    int
        Number of positions where sequences differ
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be same length: {len(seq1)} vs {len(seq2)}")
    return sum(a != b for a, b in zip(seq1, seq2))


def check_sequence_lengths(df, col_name):
    """
    Ensure all sequences in a column are the same length.

    Parameters
    ----------
    df : pd.DataFrame
    col_name : str

    Returns
    -------
    int
        The uniform sequence length

    Raises
    ------
    ValueError
        If sequences have differing lengths
    """
    sequences = df[col_name]
    seq_lengths = sequences.str.len()
    unique_lengths = seq_lengths.unique()

    if len(unique_lengths) != 1:
        length_examples = []
        for length in sorted(unique_lengths)[:5]:
            strains = df[seq_lengths == length]['representative_strain'].head(3).tolist()
            length_examples.append(f"    Length {length}: {strains}")

        raise ValueError(
            f"Column '{col_name}' has sequences of different lengths: {sorted(unique_lengths)}\n"
            f"All sequences must be the same length for proper alignment.\n"
            f"Examples:\n" +
            "\n".join(length_examples) +
            (f"\n    ... and {len(unique_lengths) - 5} more lengths" if len(unique_lengths) > 5 else "")
        )

    return int(unique_lengths[0])


def check_amino_acid_validity(df, col_name):
    """
    Ensure all sequences in a column contain only valid amino acids.

    Parameters
    ----------
    df : pd.DataFrame
    col_name : str

    Raises
    ------
    ValueError
        If any sequences contain invalid characters
    """
    sequences = df[col_name]
    invalid_mask = ~sequences.apply(has_valid_amino_acids)

    if invalid_mask.any():
        invalid_examples = []
        for idx in df[invalid_mask].index[:5]:
            row = df.loc[idx]
            seq = row[col_name]
            invalid_chars = sorted(set(seq) - VALID_AA)
            invalid_examples.append(
                f"    {row['representative_strain']}: contains invalid character(s) {invalid_chars}"
            )

        n_invalid = invalid_mask.sum()
        raise ValueError(
            f"Column '{col_name}': found {n_invalid} sequence(s) with invalid amino acids.\n"
            f"All sequences must contain only the 20 standard amino acids: {sorted(VALID_AA)}\n"
            f"Examples:\n" +
            "\n".join(invalid_examples) +
            (f"\n    ... and {n_invalid - 5} more" if n_invalid > 5 else "")
        )


def find_nearest_library_strain(sequence, current_strain_name, library_sequences, library_strains):
    """
    Find the nearest selected (library) strain to a given sequence.

    Parameters
    ----------
    sequence : str
        The sequence to compare
    current_strain_name : str
        The name of the current strain (excluded from comparison)
    library_sequences : array-like of str
        Sequences of all library strains
    library_strains : array-like of str
        Names of all library strains

    Returns
    -------
    tuple of (int, str)
        Hamming distance and strain name of nearest library strain
    """
    distances = [
        (hamming_distance(sequence, lib_seq), lib_strain)
        for lib_seq, lib_strain in zip(library_sequences, library_strains)
        if lib_strain != current_strain_name
    ]

    if not distances:
        raise ValueError(
            f"Cannot find nearest library strain for {current_strain_name}: "
            f"no other library strains available for comparison"
        )

    return min(distances)


def curate_library_sequences(valid_haplotypes_path, additional_haplotypes_paths,
                              gff_path, site_annotations_path, selected_output_path,
                              nonselected_output_path, library_fasta_path, subtype):
    """
    Curate library sequences by combining validated and additional haplotypes.

    Parameters
    ----------
    valid_haplotypes_path : str
        Path to TSV file with validated recent haplotypes (from validate_recent_haplotypes)
    additional_haplotypes_paths : list of str
        List of paths to TSV files with additional haplotypes to include
    gff_path : str
        Path to GFF file with HA1/HA2 region definitions (nucleotide coordinates)
    site_annotations_path : str
        Path to TSV file with site annotations including epitope definitions
    selected_output_path : str
        Path to output TSV of selected haplotypes
    nonselected_output_path : str
        Path to output TSV of non-selected haplotypes
    library_fasta_path : str
        Path to output FASTA file with library (selected) haplotype combined HA1+HA2 sequences
    subtype : str
        Flu subtype ('H1N1' or 'H3N2') for coordinate conversion
    """
    # Read validated recent haplotypes
    df_validated = pd.read_csv(valid_haplotypes_path, sep='\t')

    # Check that required columns exist
    required_cols = [
        'representative_strain',
        'representative_strain_ha1_sequence',
        'representative_ha2_sequence',
        'selected_haplotype',
    ]
    for col in required_cols:
        if col not in df_validated.columns:
            raise ValueError(
                f"Input file '{valid_haplotypes_path}' does not contain required column '{col}'.\n"
                f"Available columns: {list(df_validated.columns)}"
            )

    # The `selection_file` column is set upstream by validate_recent_haplotypes and flows
    # through unchanged here. Verify it's present so a missing-column failure surfaces early
    # with a clear message rather than later as a KeyError during concat.
    if 'selection_file' not in df_validated.columns:
        raise ValueError(
            f"Input file '{valid_haplotypes_path}' does not contain required column "
            f"'selection_file'. This column should be set by the upstream "
            f"validate_recent_haplotypes script.\n"
            f"Available columns: {list(df_validated.columns)}"
        )

    print(f"Read {len(df_validated)} validated haplotypes from {valid_haplotypes_path}")
    print(f"  Selected: {df_validated['selected_haplotype'].sum()}")
    print(f"  Non-selected: {(~df_validated['selected_haplotype']).sum()}")
    print(f"  selection_file value counts (non-null):")
    sf_counts = df_validated['selection_file'].dropna().value_counts()
    if len(sf_counts) == 0:
        print(f"    (none)")
    else:
        for val, count in sf_counts.items():
            print(f"    {val}: {count}")

    # Read additional haplotypes if any
    dfs_to_combine = [df_validated]

    if additional_haplotypes_paths:
        for path in additional_haplotypes_paths:
            df_additional = pd.read_csv(path, sep='\t')

            # Detect input format
            input_format = detect_input_format(df_additional)
            print(f"\nRead {len(df_additional)} additional haplotypes from {path}")
            print(f"  Input format: {input_format}")

            # Process based on format
            if input_format == 'combined':
                print(f"  Extracting HA1/HA2 from combined sequence using GFF: {gff_path}")
                df_additional = process_additional_haplotypes(df_additional, gff_path, input_format, subtype)
                print(f"  Successfully extracted HA1 and HA2 regions")

            # Ensure required columns exist after processing
            for col in ['representative_strain', 'representative_strain_ha1_sequence', 'representative_ha2_sequence']:
                if col not in df_additional.columns:
                    raise ValueError(
                        f"Additional haplotypes file '{path}' does not contain required column '{col}'.\n"
                        f"Available columns: {list(df_additional.columns)}"
                    )

            df_additional['selected_haplotype'] = True

            # Tag with basename of this additional file
            df_additional['selection_file'] = os.path.basename(path)

            print(f"  All additional haplotypes are marked as selected")
            print(f"  selection_file: {os.path.basename(path)}")

            dfs_to_combine.append(df_additional)

    # Combine all dataframes
    df_combined = pd.concat(dfs_to_combine, ignore_index=True)

    print(f"\nCombined dataframe contains {len(df_combined)} total haplotypes")
    print(f"  Selected: {df_combined['selected_haplotype'].sum()}")
    print(f"  Non-selected: {(~df_combined['selected_haplotype']).sum()}")

    # Create combined HA1+HA2 sequence column early so it flows through all downstream steps.
    # Both individual sequence columns must be present and non-null at this point (the upstream
    # validation script should have filtered any rows with missing HA2 to the invalid TSV).
    if df_combined['representative_ha2_sequence'].isna().any():
        n_missing = df_combined['representative_ha2_sequence'].isna().sum()
        raise ValueError(
            f"{n_missing} haplotype(s) have a missing HA2 sequence. "
            f"Only haplotypes with valid HA1 and HA2 sequences should reach this script. "
            f"Check that the input comes from the valid_haplotypes output of validate_haplotypes."
        )

    df_combined['representative_strain_ha1_ha2_sequence'] = (
        df_combined['representative_strain_ha1_sequence'] +
        df_combined['representative_ha2_sequence']
    )

    print(f"\nCreated combined HA1+HA2 sequence column 'representative_strain_ha1_ha2_sequence'")

    # Check for duplicate representative_strain values
    duplicates = df_combined['representative_strain'].duplicated(keep=False)
    if duplicates.any():
        dup_values = df_combined[duplicates]['representative_strain'].value_counts()
        examples = []
        for val, count in dup_values.head(10).items():
            dup_rows = df_combined[df_combined['representative_strain'] == val].index.tolist()
            examples.append(f"    '{val}' appears {count} times (rows: {dup_rows})")

        raise ValueError(
            f"Combined data contains {len(dup_values)} duplicate representative_strain value(s).\n"
            f"Each strain name must be unique.\n"
            f"Examples:\n" +
            "\n".join(examples) +
            (f"\n    ... and {len(dup_values) - 10} more" if len(dup_values) > 10 else "")
        )

    print("  ✓ All representative_strain values are unique")

    # Validate sequence lengths and amino acid validity for all three sequence columns
    for col in [
        'representative_strain_ha1_sequence',
        'representative_ha2_sequence',
        'representative_strain_ha1_ha2_sequence',
    ]:
        seq_length = check_sequence_lengths(df_combined, col)
        print(f"  ✓ All '{col}' sequences have length {seq_length}")
        check_amino_acid_validity(df_combined, col)
        print(f"  ✓ All '{col}' sequences contain only valid amino acids")

    # Split into selected and non-selected
    df_selected = df_combined[df_combined['selected_haplotype']].copy()

    # Get library sequences and strain names for nearest-strain searches
    ha1_library_sequences = df_selected['representative_strain_ha1_sequence'].values
    ha1_ha2_library_sequences = df_selected['representative_strain_ha1_ha2_sequence'].values
    library_strains = df_selected['representative_strain'].values

    # Compute HA1 nearest library strain (independent search)
    print(f"\nComputing HA1 Hamming distances to nearest library strain...")
    ha1_distances = []
    ha1_nearest_strains = []

    for _, row in df_combined.iterrows():
        dist, nearest = find_nearest_library_strain(
            row['representative_strain_ha1_sequence'],
            row['representative_strain'],
            ha1_library_sequences,
            library_strains,
        )
        ha1_distances.append(dist)
        ha1_nearest_strains.append(nearest)

    df_combined['hamming_distance_nearest_library_strain_ha1'] = ha1_distances
    df_combined['nearest_library_strain_ha1'] = ha1_nearest_strains
    print(f"  Computed HA1 Hamming distances for {len(df_combined)} haplotypes")

    # Compute HA1+HA2 nearest library strain (independent search)
    print(f"\nComputing HA1+HA2 Hamming distances to nearest library strain...")
    ha1_ha2_distances = []
    ha1_ha2_nearest_strains = []

    for _, row in df_combined.iterrows():
        dist, nearest = find_nearest_library_strain(
            row['representative_strain_ha1_ha2_sequence'],
            row['representative_strain'],
            ha1_ha2_library_sequences,
            library_strains,
        )
        ha1_ha2_distances.append(dist)
        ha1_ha2_nearest_strains.append(nearest)

    df_combined['hamming_distance_nearest_library_strain_ha1_ha2'] = ha1_ha2_distances
    df_combined['nearest_library_strain_ha1_ha2'] = ha1_ha2_nearest_strains
    print(f"  Computed HA1+HA2 Hamming distances for {len(df_combined)} haplotypes")

    # Update selected for distance reporting
    df_selected = df_combined[df_combined['selected_haplotype']].copy()
    df_nonselected = df_combined[~df_combined['selected_haplotype']].copy()

    for dist_col in [
        'hamming_distance_nearest_library_strain_ha1',
        'hamming_distance_nearest_library_strain_ha1_ha2',
    ]:
        print(f"  Selected {dist_col}: range {df_selected[dist_col].min()}-{df_selected[dist_col].max()}")
        if len(df_nonselected) > 0:
            print(f"  Non-selected {dist_col}: range {df_nonselected[dist_col].min()}-{df_nonselected[dist_col].max()}")

    # Read site annotations to identify HA1 epitope sites
    print(f"\nReading site annotations from {site_annotations_path}")
    df_sites = pd.read_csv(site_annotations_path, sep='\t')

    # Filter to HA1 rows only for epitope position extraction
    if 'protein' not in df_sites.columns:
        raise ValueError(
            f"Site annotations file {site_annotations_path} does not contain a 'protein' column.\n"
            f"Available columns: {df_sites.columns.tolist()}"
        )

    df_sites_ha1 = df_sites[df_sites['protein'] == 'HA1']

    # Identify epitope columns (columns ending with "_epitope")
    epitope_cols = [col for col in df_sites_ha1.columns if col.endswith('_epitope')]

    if not epitope_cols:
        raise ValueError(
            f"No epitope columns found in site annotations file {site_annotations_path}.\n"
            f"Available columns: {df_sites.columns.tolist()}"
        )

    print(f"  Found {len(epitope_cols)} HA1 epitope sets: {epitope_cols}")

    # Extract 0-indexed positions for each epitope from HA1 rows only
    epitope_positions = {}
    for epitope_col in epitope_cols:
        positions = (df_sites_ha1[df_sites_ha1[epitope_col] == 1]['sequential_site'] - 1).tolist()
        epitope_positions[epitope_col] = positions
        print(f"    {epitope_col}: {len(positions)} HA1 sites")

    # Compute epitope-specific distances using HA1 sequences and HA1 nearest library strain
    print(f"\nComputing HA1 epitope-specific distances to nearest HA1 library strain...")

    # Build a lookup from strain name to HA1 sequence for the nearest-strain dereference
    ha1_seq_by_strain = dict(zip(
        df_selected['representative_strain'],
        df_selected['representative_strain_ha1_sequence'],
    ))

    epitope_distance_cols = {}
    for epitope_col in epitope_cols:
        distance_col_name = epitope_col.replace('_epitope', '_epitope_distance')
        epitope_distance_cols[epitope_col] = distance_col_name

        positions = epitope_positions[epitope_col]

        if not positions:
            print(f"  Warning: {epitope_col} has no HA1 epitope sites, setting all distances to 0")
            df_combined[distance_col_name] = 0
            continue

        epitope_distances = []
        for _, row in df_combined.iterrows():
            seq = row['representative_strain_ha1_sequence']
            nearest_seq = ha1_seq_by_strain[row['nearest_library_strain_ha1']]
            dist = sum(seq[pos] != nearest_seq[pos] for pos in positions)
            epitope_distances.append(dist)

        df_combined[distance_col_name] = epitope_distances
        print(f"  {distance_col_name}: computed for {len(df_combined)} haplotypes")

    # Update selected and non-selected after all columns are added
    df_selected = df_combined[df_combined['selected_haplotype']].copy()
    df_nonselected = df_combined[~df_combined['selected_haplotype']].copy()

    for epitope_col in epitope_cols:
        distance_col_name = epitope_distance_cols[epitope_col]
        print(f"  Selected {distance_col_name}: range {df_selected[distance_col_name].min()}-{df_selected[distance_col_name].max()}")
        if len(df_nonselected) > 0:
            print(f"  Non-selected {distance_col_name}: range {df_nonselected[distance_col_name].min()}-{df_nonselected[distance_col_name].max()}")

    # Reorder columns: identifiers, distances, sequences, then remaining.
    # `selection_file` is intentionally pinned to the very end of the selected TSV
    # (and dropped entirely from the non-selected TSV).
    epitope_distance_col_names = sorted(epitope_distance_cols.values())

    priority_cols = [
        'representative_strain',
        'derived_haplotype',
        'hamming_distance_nearest_library_strain_ha1',
        'nearest_library_strain_ha1',
        'hamming_distance_nearest_library_strain_ha1_ha2',
        'nearest_library_strain_ha1_ha2',
    ] + epitope_distance_col_names + [
        'representative_strain_ha1_sequence',
        'representative_ha2_sequence',
        'representative_strain_ha1_ha2_sequence',
    ]

    missing_priority_cols = [col for col in priority_cols if col not in df_combined.columns]
    if missing_priority_cols:
        raise ValueError(
            f"Expected columns are missing from combined dataframe: {missing_priority_cols}\n"
            f"Available columns: {list(df_combined.columns)}"
        )

    other_cols = [
        col for col in df_combined.columns
        if col not in priority_cols and col != 'selection_file'
    ]
    # selected TSV: priority cols, remaining cols, then selection_file pinned at the end
    selected_ordered_cols = priority_cols + other_cols + ['selection_file']
    # non-selected TSV: same order but selection_file dropped entirely
    nonselected_ordered_cols = priority_cols + other_cols

    df_selected = df_selected[selected_ordered_cols]
    df_nonselected = df_nonselected[nonselected_ordered_cols]

    # Sort: selected ascending by HA1+HA2 distance, non-selected descending
    df_selected = df_selected.sort_values('hamming_distance_nearest_library_strain_ha1_ha2', ascending=True)
    df_nonselected = df_nonselected.sort_values('hamming_distance_nearest_library_strain_ha1_ha2', ascending=False)

    # Save outputs
    df_selected.to_csv(selected_output_path, sep='\t', index=False)
    df_nonselected.to_csv(nonselected_output_path, sep='\t', index=False)

    # Write library FASTA using combined HA1+HA2 sequences
    with open(library_fasta_path, 'w') as f:
        for _, row in df_selected.iterrows():
            f.write(f">{row['representative_strain']}\n")
            f.write(f"{row['representative_strain_ha1_ha2_sequence']}\n")

    print(f"\nOutputs:")
    print(f"  Selected haplotypes: {selected_output_path}")
    print(f"    {len(df_selected)} haplotypes")
    print(f"    selection_file value counts:")
    for val, count in df_selected['selection_file'].value_counts().items():
        print(f"      {val}: {count}")
    print(f"  Non-selected haplotypes: {nonselected_output_path}")
    print(f"    {len(df_nonselected)} haplotypes")
    print(f"  Library FASTA: {library_fasta_path}")
    print(f"    {len(df_selected)} sequences (combined HA1+HA2)")

    print(f"\nCuration complete ✓")


if __name__ == '__main__':
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    curate_library_sequences(
        valid_haplotypes_path=snakemake.input.valid_haplotypes,
        additional_haplotypes_paths=snakemake.input.additional_haplotypes,
        gff_path=snakemake.input.gff,
        site_annotations_path=snakemake.input.site_annotations,
        selected_output_path=snakemake.output.selected_haplotypes,
        nonselected_output_path=snakemake.output.nonselected_haplotypes,
        library_fasta_path=snakemake.output.library_fasta,
        subtype=snakemake.wildcards.subtype
    )