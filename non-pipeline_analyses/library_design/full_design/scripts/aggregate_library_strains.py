"""
Aggregate library strains from multiple subtypes into a single TSV.

This script combines the match_genbank results from multiple subtypes (H1N1, H3N2)
into a single output file with validation checks. For each subtype, it also joins
against the corresponding curation TSV to bring in the `selection_file` and
`derived_haplotype` columns, matching on protein sequence and strain name (curation
`representative_strain_ha1_ha2_sequence` + `representative_strain` ==
match_genbank `prot_sequence` + `strain`).
"""

import pandas as pd
import sys


def build_curation_lookup(curated_tsv_path, subtype):
    """
    Build a {(prot_sequence, strain): (selection_file, derived_haplotype)} lookup from a curated TSV.

    The lookup uses (protein_sequence, strain_name) as a composite key so that if the
    same protein sequence appears under multiple strain names (e.g., one from a selection
    file and one from an additional haplotype file), each strain maps to its own
    `selection_file` and `derived_haplotype` values.

    Parameters
    ----------
    curated_tsv_path : str
        Path to a curation selected-haplotypes TSV.
    subtype : str
        Subtype label, used only in error messages.

    Returns
    -------
    dict
        Mapping from (representative_strain_ha1_ha2_sequence, representative_strain) to
        (selection_file, derived_haplotype).
    """
    df = pd.read_csv(curated_tsv_path, sep='\t')

    required = ['representative_strain_ha1_ha2_sequence', 'representative_strain', 'selection_file', 'derived_haplotype']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Curated TSV '{curated_tsv_path}' (subtype {subtype}) is missing required "
            f"column(s): {missing}\nAvailable columns: {list(df.columns)}"
        )

    lookup = {}

    for _, row in df.iterrows():
        seq = row['representative_strain_ha1_ha2_sequence']
        strain = row['representative_strain']
        selection_file = row['selection_file']
        derived_haplotype = row['derived_haplotype']
        
        # Use (sequence, strain) as the composite key
        lookup[(seq, strain)] = (selection_file, derived_haplotype)

    return lookup


def aggregate_library_strains(input_csvs, curated_tsvs, output_tsv, subtypes):
    """
    Aggregate library strains from multiple subtype-specific CSVs, joining each
    against its corresponding curation TSV to attach the `selection_file` and
    `derived_haplotype` columns.

    Parameters
    ----------
    input_csvs : list of str
        List of paths to input CSV files from match_genbank rule.
    curated_tsvs : list of str
        List of paths to curation selected-haplotypes TSVs, parallel to input_csvs
        and subtypes (one per subtype).
    output_tsv : str
        Path to output TSV file.
    subtypes : list of str
        List of subtypes corresponding to each input CSV.
    """
    if not (len(input_csvs) == len(curated_tsvs) == len(subtypes)):
        raise ValueError(
            f"input_csvs ({len(input_csvs)}), curated_tsvs ({len(curated_tsvs)}), "
            f"and subtypes ({len(subtypes)}) must all be the same length"
        )

    # Read and process each subtype
    dfs = []
    for csv_path, curated_path, subtype in zip(input_csvs, curated_tsvs, subtypes):
        df = pd.read_csv(csv_path)

        # Add subtype column
        df['subtype'] = subtype

        # Build per-subtype lookup and attach selection_file and derived_haplotype
        # via (prot_sequence, strain) join
        lookup = build_curation_lookup(curated_path, subtype)
        
        # Apply the lookup to each row, returning both columns
        join_results = df.apply(
            lambda row: lookup.get((row['prot_sequence'], row['strain']), (pd.NA, pd.NA)),
            axis=1,
            result_type='expand'
        )
        df['selection_file'] = join_results[0]
        df['derived_haplotype'] = join_results[1]

        # Every strain must have a matching selection_file and derived_haplotype
        unmatched_mask = df['selection_file'].isna()
        if unmatched_mask.any():
            unmatched = df.loc[unmatched_mask, ['strain', 'prot_sequence']].copy()
            unmatched['subtype'] = subtype
            msg_lines = [
                f"ERROR: {unmatched_mask.sum()} strain(s) in '{csv_path}' (subtype {subtype}) "
                f"could not be matched to any curated row via (prot_sequence, strain) pair."
            ]
            for idx, row in unmatched.head(20).iterrows():
                msg_lines.append(f"  - {row['subtype']} / {row['strain']} / seq_len={len(row['prot_sequence'])}")
            if len(unmatched) > 20:
                msg_lines.append(f"  ... and {len(unmatched) - 20} more")
            raise ValueError("\n".join(msg_lines))

        # Select only required columns with derived_haplotype after strain
        df = df[[
            'subtype',
            'strain',
            'derived_haplotype',
            'accession_w_aa_muts_added',
            'prot_sequence',
            'nt_sequence',
            'selection_file',
        ]]

        dfs.append(df)

    # Combine all dataframes
    combined_df = pd.concat(dfs, ignore_index=True)

    # Validation 1: Check for null values
    null_counts = combined_df.isnull().sum()
    if null_counts.any():
        null_info = null_counts[null_counts > 0]
        error_msg = "ERROR: Found null values in the following columns:\n"
        for col, count in null_info.items():
            error_msg += f"  - {col}: {count} null values\n"

        # Show rows with null values for debugging
        error_msg += "\nRows with null values:\n"
        null_rows = combined_df[combined_df.isnull().any(axis=1)]
        error_msg += null_rows.to_string()

        raise ValueError(error_msg)

    # Validation 2: Check for duplicate strain names
    duplicate_strains = combined_df['strain'].duplicated()
    if duplicate_strains.any():
        duplicate_names = combined_df.loc[duplicate_strains, 'strain'].unique()
        error_msg = f"ERROR: Found {len(duplicate_names)} duplicate strain names:\n"
        for name in duplicate_names:
            count = (combined_df['strain'] == name).sum()
            error_msg += f"  - {name} (appears {count} times)\n"

        raise ValueError(error_msg)

    # Write output
    combined_df.to_csv(output_tsv, sep='\t', index=False)

    # Print summary statistics
    print(f"Successfully aggregated library strains:")
    print(f"  - Total strains: {len(combined_df)}")
    for subtype in subtypes:
        count = (combined_df['subtype'] == subtype).sum()
        print(f"  - {subtype}: {count} strains")
    print(f"  - selection_file value counts:")
    for val, count in combined_df['selection_file'].value_counts().items():
        print(f"      {val}: {count}")
    print(f"  - Output written to: {output_tsv}")


if __name__ == "__main__":
    # Redirect stdout and stderr to log file
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")

    # Snakemake workflow will provide these variables
    aggregate_library_strains(
        input_csvs=snakemake.input.csvs,
        curated_tsvs=snakemake.input.curated_tsvs,
        output_tsv=snakemake.output.tsv,
        subtypes=snakemake.params.subtypes
    )