# CLAUDE.md 

Also read [seqneut-pipeline/CLAUDE.md](seqneut-pipeline/CLAUDE.md), which covers the
pipeline included here as a submodule, and follow its principles as well.

## Overrides of the lab coding standards

The `age` column in `data/sera_metadata/*.csv` holds either an exact age or an age
range, depending on what the cohort released. Both forms are accepted and documented in
`README.md`, and `scripts/aggregate_sera_metadata.py` parses them into `age_numeric`.
This is a deliberate exception to the standards' "one column, one job" rule: do not flag
it, and do not propose splitting it into separate bound columns.

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

## Organization of Snakemake Rules

The core analysis is the `seqneut-pipeline` submodule; do not add rules specific to this
project to it. Analyses in `snakemake` that are outside that core pipeline are organized
like this:

- Each coherent chunk of analysis gets its own `.smk` file in `rules/`, which `Snakefile`
  includes after the pipeline include. `Snakefile` itself holds only the `configfile`,
  the `include` statements, and `rule all`.
- Each file in `rules/` ends by defining the list of final outputs it contributes (e.g.
  `trees_outputs`), which `rule all` collects. This mirrors `seqneut_pipeline_outputs`,
  and keeps `rule all` from having to know path patterns defined elsewhere.
- `snakemake` resolves the paths in `script:`, `conda:`, and `module: snakefile:`
  relative to the file that defines the rule, so from `rules/` these need a `../`
  prefix. Input, output, and log paths are relative to the working directory instead,
  and so do not.
- Included files all share one global namespace, so prefix helper variables with `_`.
- Defining a helper function alongside rules in the same `.smk` fails `snakemake --lint`
  ("Mixed rules and functions in same snakefile"). Derive a variable from the config
  instead, as `rules/library_qc.smk` does with `pool_plates`; a rule can then index it
  from a `lambda wc:` without repeating the lookup.
- The top-level config keys that select which library QC analyses to run (`analyze_pools`,
  `analyze_repools`, `analyze_single_well_infections`) are read with a default of empty, so
  a configuration that omits one, or leaves it blank, runs no jobs for it rather than
  failing. This is a deliberate exception to the fail-fast principle above, since not every
  project does every kind of QC. It applies **only** to those keys: everything *inside* one
  of those sections is still required, and a missing or unexpected key there must raise.
  `rules/library_qc.smk` applies the default in one place, at the top, and reads the
  resulting variables everywhere else.

### Testing a rule's script without running the pipeline

A `script:` rule reads a global `snakemake` object, so its script can be run directly
against a stub of one. This exercises the real script rather than a copy of it, and reaches
failure cases the committed data does not contain, such as a well with no counts or a
malformed configuration value:

```python
snakemake = types.SimpleNamespace(
    wildcards=types.SimpleNamespace(pool="test_pool"),
    params=types.SimpleNamespace(...),  # as the rule's `params` would be
    input=types.SimpleNamespace(...),   # point at copies of the real inputs, edited to
    output=types.SimpleNamespace(...),  # fabricate the failure being tested
    log=[str(out / "log.txt")],
)
builtins.snakemake = snakemake
globals_ = runpy.run_path("scripts/<script>.py", run_name="__main__")
```

Two things to get right. These scripts redirect `sys.stdout` and `sys.stderr` to their log
file, so save and restore both around the call or the harness loses its own output. And
that log file is left open and buffered, so flush it through the globals `runpy` returns
(`globals_["log"].flush()`) before reading it, or it will look empty.

**IMPORTANT: rules outside the core pipeline must be kept concise and non-redundant.**
They are glue, and everything below follows from that:

- A rule should be a thin wrapper that names its inputs and outputs, passes config values
  through `params`, and delegates the analysis itself to a script in `scripts/` or a
  notebook in `notebooks/`. Do not put analysis logic in a `.smk` file.
- Never write the same path pattern twice. Refer to another rule's output as
  `rules.<rule>.output.<name>` rather than repeating the string, and build target lists
  with `expand` over `rules.<rule>.output.<name>`.
- Consume an existing rule's output rather than recomputing the same quantity, and add a
  wildcard to an existing rule rather than copying it into a near-duplicate rule.
- Do not add a config key or a rule whose only purpose is to switch an analysis off; let
  the script or notebook report that there is nothing to show.

## Keeping README.md Current and Proportionate

`README.md` must be updated along with any change that adds, removes, or redirects an
analysis, an input data format, or a result file. It describes the key points only, and
its structure carries meaning:

- The early sections describe the underlying study: the assay, the viral library, the
  sera, and the titers. They are about what the data *are*, not about the analyses run on
  them; do not add an analysis description there.
- Each analysis outside the core pipeline gets one short subsection under "Additional
  analyses" near the end. Keep these evenly weighted and of the same shape: what the
  analysis does, which `.smk` file runs it, where its outputs go, and where it is
  configured. A new analysis that needs more room than its neighbors is a sign it is
  being described in too much detail, not that it deserves a longer section.
- Say where something is configured, never what the current configuration is (see the
  Documentation Principle above). Config values, thresholds, and file lists live in
  `config.yml` and would go stale here.
- Describe a thing in exactly one place and cross-reference it from anywhere else that
  needs it, rather than restating it.

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
