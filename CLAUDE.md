# CLAUDE.md 

## Critical Scientific Coding Principles

**This is scientific research code.** Data integrity and reproducibility are paramount. Follow these principles:

### 1. Fail Fast - No Silent Errors
- **NEVER allow silent failures** or default to placeholder values
- All data processing should raise explicit exceptions when issues are encountered
- Validate inputs at entry points (file loading, configuration parsing)
- Use assertions for critical assumptions
- Log warnings for unexpected but non-fatal conditions

**Example - Good**:
```python
if barcode_counts.sum() < min_counts:
    raise ValueError(f"Barcode counts {barcode_counts.sum()} below minimum {min_counts}")
```

**Example - Bad**:
```python
if barcode_counts.sum() < min_counts:
    print("Warning: low counts")  # Silent failure - DO NOT DO THIS
    barcode_counts = None  # Might cause issues downstream
```

### 2. Single Source of Truth (DRY Principle)
- **Parameters should be specified in exactly ONE place** (typically `config.yml`)
- Never duplicate parameter values in code, documentation, or multiple config sections
- If a parameter exists in config, reference it - don't redefine it
- This prevents inconsistencies and improves maintainability

**Example - Good**:
```python
min_counts = config["qc_thresholds"]["min_counts"]
```

**Example - Bad**:
```python
min_counts = 500  # Duplicates value from config - DO NOT DO THIS
```

**Documentation Principle**:
- Both code documentation and README should **reference where** values are set (e.g., "see `config.yml`", "configured in `data/viral_libraries/`")
- Do NOT repeat current configuration values in Markdown text
- Describe WHAT parameters control and HOW to set them, not their current values
- This keeps documentation maintainable as configuration changes

**Example - Good**:
```markdown
QC thresholds are configured in `config.yml` under `default_serum_qc_thresholds`.
Key parameters include `min_replicates` and `max_fold_change_from_median`.
```

**Example - Bad**:
```markdown
The minimum replicates threshold is set to 1, and outliers are flagged at 3-fold change.
```

### 3. Explicit Over Implicit
- Be explicit about data transformations and filtering
- Document QC drops in YAML files (already implemented in pipeline)
- Avoid "magic numbers" - use named configuration parameters
- Type hints and docstrings for complex functions

### 4. Reproducibility
- All analysis controlled by `config.yml`
- Random seeds set where stochastic methods used
- Track exact versions (conda environment, submodule versions)
- Results committed to git for key QC files

### Non-Pipeline Analyses (GENERALLY IGNORE)

The `non-pipeline_analyses/` directory contains one-off analyses for library design and pooling optimization. These are **NOT part of the main neutralization assay pipeline** and are documented separately in that directory. Ignore unless the user specifically asks about them.

## Code Style and Quality Requirements

All code must pass these checks before committing:

### Python
```bash
ruff check .        # Linting (fast, comprehensive)
black .             # Code formatting (auto-fix)
```

### Snakemake
```bash
snakefmt .          # Snakemake formatting
snakemake --lint    # Snakemake validation
```

### Configuration
**File**: `config.yml`
- Single source of truth for all pipeline parameters
