"""Aggregate and validate sera metadata from multiple cohorts.

Reads sera metadata CSVs from multiple cohorts, validates required columns and
data formats, standardizes values, and outputs a single aggregated CSV.
"""

import re
import sys

import pandas as pd

sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

# Get inputs and outputs from snakemake
input_csvs = snakemake.input.csvs
output_csv = snakemake.output.csv

# Required columns that must exist in each input file
REQUIRED_COLUMNS = [
    "bloom_lab_id",
    "cohort",
    "species",
    "age",
    "sex",
    "collection_date",
]

# Collection dates are recorded to month precision only
COLLECTION_DATE_REGEX = r"\d{4}-(0[1-9]|1[0-2])"


def parse_age_to_numeric(age_str):
    """Parse age string to numeric midpoint in years.

    Handles formats:
    - Numeric (e.g., "40"): returns float
    - Range with 'y' suffix (e.g., "10-19y"): returns midpoint
    - Range without 'y' (e.g., "18-29"): returns midpoint
    - Open-ended (e.g., "75+"): returns lower bound

    Raises ValueError for unrecognized formats.
    """
    age_str = str(age_str).strip()

    # Try numeric first
    try:
        return float(age_str)
    except ValueError:
        pass

    # Remove trailing 'y' if present
    if age_str.lower().endswith("y"):
        age_str = age_str[:-1]

    # Handle open-ended ranges (e.g., "75+")
    if age_str.endswith("+"):
        try:
            return float(age_str[:-1])
        except ValueError:
            raise ValueError(f"Cannot parse open-ended age: {age_str}")

    # Handle ranges (e.g., "10-19" or "18-29")
    range_match = re.match(r"^(\d+)-(\d+)$", age_str)
    if range_match:
        low = float(range_match.group(1))
        high = float(range_match.group(2))
        assert high > low, f"{low=}, {high=}, {age_str=}"
        return (low + high) / 2

    raise ValueError(f"Cannot parse age: {age_str}")


def normalize_sex(sex_str):
    """Normalize sex string to 'Male', 'Female', or 'Unknown'.

    Handles various input formats (M/F, male/female, Male/Female, etc.).
    Raises ValueError for unrecognized values.
    """
    sex_str = str(sex_str).strip().lower()

    if sex_str in ("m", "male"):
        return "Male"
    elif sex_str in ("f", "female"):
        return "Female"
    elif sex_str in ("unknown", ""):
        return "Unknown"
    else:
        raise ValueError(f"Cannot normalize sex: {sex_str}")


# =============================================================================
# Load and validate each input file
# =============================================================================
print(f"Aggregating sera metadata from {len(input_csvs)} files...")
print()

all_dfs = []

for csv_path in input_csvs:
    print(f"Processing: {csv_path}")

    df = pd.read_csv(csv_path)
    print(f"  Loaded {len(df)} rows")

    # Check required columns exist
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in {csv_path}: {missing_cols}")

    # Check no existing 'serum' column (will rename bloom_lab_id)
    if "serum" in df.columns:
        raise ValueError(
            f"File {csv_path} already has 'serum' column; "
            "expected 'bloom_lab_id' which will be renamed"
        )

    # Strip whitespace from all string columns
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.strip()
            # Convert "nan" strings back to actual NaN
            df[col] = df[col].replace("nan", pd.NA)

    # Validate bloom_lab_id is non-null and unique within this file
    if df["bloom_lab_id"].isna().any():
        null_count = df["bloom_lab_id"].isna().sum()
        raise ValueError(f"File {csv_path} has {null_count} null bloom_lab_id values")

    duplicates = df["bloom_lab_id"].duplicated()
    if duplicates.any():
        dup_ids = df.loc[duplicates, "bloom_lab_id"].tolist()[:5]
        raise ValueError(
            f"File {csv_path} has duplicate bloom_lab_id values: {dup_ids}"
        )

    # Validate cohort is non-null
    if df["cohort"].isna().any():
        null_count = df["cohort"].isna().sum()
        raise ValueError(f"File {csv_path} has {null_count} null cohort values")

    # Validate species is non-null
    if df["species"].isna().any():
        null_count = df["species"].isna().sum()
        raise ValueError(f"File {csv_path} has {null_count} null species values")

    # Validate collection_date is YYYY-MM (null dates fail this too)
    bad_dates = df.loc[
        ~df["collection_date"].str.fullmatch(COLLECTION_DATE_REGEX, na=False),
        "collection_date",
    ]
    if len(bad_dates):
        raise ValueError(
            f"File {csv_path} has collection_date values not in YYYY-MM format: "
            f"{bad_dates.unique().tolist()}"
        )

    all_dfs.append(df)
    print("  Validated successfully")

print()

# =============================================================================
# Combine all dataframes
# =============================================================================
combined_df = pd.concat(all_dfs, ignore_index=True)
print(f"Combined {len(combined_df)} total rows")

# Check bloom_lab_id is globally unique
duplicates = combined_df["bloom_lab_id"].duplicated()
if duplicates.any():
    dup_ids = combined_df.loc[duplicates, "bloom_lab_id"].tolist()[:5]
    raise ValueError(f"Duplicate bloom_lab_id values across files: {dup_ids}")

print("All bloom_lab_id values are unique across all files")
print()

# =============================================================================
# Normalize and standardize values
# =============================================================================
print("Normalizing values...")

# Rename bloom_lab_id to serum
combined_df = combined_df.rename(columns={"bloom_lab_id": "serum"})

# Normalize sex
print("  Normalizing sex values...")
try:
    combined_df["sex"] = combined_df["sex"].apply(normalize_sex)
except ValueError as e:
    raise ValueError(f"Sex normalization failed: {e}")

# Parse age to numeric
print("  Parsing age to numeric...")
try:
    combined_df["age_numeric"] = combined_df["age"].apply(parse_age_to_numeric)
except ValueError as e:
    raise ValueError(f"Age parsing failed: {e}")

print()

# =============================================================================
# Summary statistics
# =============================================================================
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()

# Sera per cohort
print("Sera per cohort:")
cohort_counts = combined_df["cohort"].value_counts().sort_index()
for cohort, count in cohort_counts.items():
    print(f"  {cohort}: {count}")
print()
print(f"Total sera: {len(combined_df)}")
print()

# Age distribution
print("Age distribution:")
print(f"  Min: {combined_df['age_numeric'].min():.1f}")
print(f"  Max: {combined_df['age_numeric'].max():.1f}")
print(f"  Mean: {combined_df['age_numeric'].mean():.1f}")
print(f"  Median: {combined_df['age_numeric'].median():.1f}")
print()

# Sex distribution
print("Sex distribution:")
sex_counts = combined_df["sex"].value_counts()
for sex, count in sex_counts.items():
    print(f"  {sex}: {count}")
print()

# =============================================================================
# Output
# =============================================================================

# Reorder columns to put serum first, then required columns, then age_numeric,
# then any additional columns
output_cols = [
    "serum",
    "cohort",
    "species",
    "age",
    "sex",
    "collection_date",
    "age_numeric",
]
other_cols = [col for col in combined_df.columns if col not in output_cols]
output_cols.extend(other_cols)

combined_df = combined_df[output_cols]

# Write output
combined_df.to_csv(output_csv, index=False)
print(f"Output written to: {output_csv}")

# Flush log before reading it back
log.flush()
