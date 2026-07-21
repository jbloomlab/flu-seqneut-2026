"""
Aggregate library strains from multiple subtypes into a single TSV.
This script combines the match_genbank results from multiple subtypes (H1N1, H3N2)
into a single output file with validation checks, and enriches with selection_file
information from the valid haplotypes file.
"""
import pandas as pd
import sys

def aggregate_library_strains(input_csvs, input_csvs_valid_haplotypes, output_tsv, subtypes):
    """
    Aggregate library strains from multiple subtype-specific CSVs.
    
    Parameters
    ----------
    input_csvs : list of str
        List of paths to input CSV files from match_genbank rule
    input_csvs_valid_haplotypes : list of str
        List of paths to valid haplotypes TSV files containing selection_file
    output_tsv : str
        Path to output TSV file
    subtypes : list of str
        List of subtypes corresponding to each input CSV
    """
    # Read and combine valid haplotypes files to get selection_file mapping
    haplotype_dfs = []
    for hap_path in input_csvs_valid_haplotypes:
        hap_df = pd.read_csv(hap_path, sep='\t')
        # Keep only the columns needed for matching plus the metadata we carry
        # through (selection_file, derived_haplotype, and latest_sequence).
        hap_df = hap_df[['representative_strain', 'representative_strain_ha_sequence', 'selection_file', 'derived_haplotype', 'latest_sequence']]
        haplotype_dfs.append(hap_df)
    
    # Combine all haplotype data
    combined_haplotypes = pd.concat(haplotype_dfs, ignore_index=True)
    
    # Read and process each subtype
    dfs = []
    for csv_path, subtype in zip(input_csvs, subtypes):
        df = pd.read_csv(csv_path)
        
        # Add subtype column
        df['subtype'] = subtype
        
        # Merge with haplotypes to get selection_file
        # Match on strain (from input csv) with representative_strain (from haplotypes)
        # and prot_sequence with representative_strain_ha_sequence
        df = df.merge(
            combined_haplotypes,
            left_on=['strain', 'prot_sequence'],
            right_on=['representative_strain', 'representative_strain_ha_sequence'],
            how='left'
        )
        
        # Every strain must have a matching selection_file
        unmatched = df[df['selection_file'].isnull()]
        if len(unmatched) > 0:
            msg_lines = [f"ERROR: {len(unmatched)} strain(s) in {subtype} could not be matched to any haplotype row:"]
            for _, row in unmatched.head(20).iterrows():
                msg_lines.append(f"  - {row['strain']} / seq_len={len(row['prot_sequence'])}")
            if len(unmatched) > 20:
                msg_lines.append(f"  ... and {len(unmatched) - 20} more")
            raise ValueError("\n".join(msg_lines))
        
        # Select required columns, with derived_haplotype after strain and
        # latest_sequence last, matching the full_design aggregated schema.
        df = df[['subtype', 'strain', 'derived_haplotype', 'accession_w_aa_muts_added', 'prot_sequence', 'nt_sequence', 'selection_file', 'latest_sequence']]
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
    print(f"  - Output written to: {output_tsv}")


if __name__ == "__main__":
    # Redirect stdout and stderr to log file
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    
    # Snakemake workflow will provide these variables
    aggregate_library_strains(
        input_csvs=snakemake.input.csvs,
        input_csvs_valid_haplotypes=snakemake.input.csvs_valid_haplotypes,
        output_tsv=snakemake.output.tsv,
        subtypes=snakemake.params.subtypes
    )