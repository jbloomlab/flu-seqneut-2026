# Analyze single virus-per well infections to look for contamination

This directory analyzes single-virus-per-well infections to look for contamination.

[notebook.py](notebook.py) classifies each barcode read in each well as coming from that well's own designated strain, an adjacent well, a distant well, or not in the plate at all — flagging likely contamination. Runs across one or more plates (see `plate_configs`).

**Inputs** (per plate, in `plate_configs`): `countsdir` (counts/fates CSVs), `samplesfile` (well metadata), `viral_library_file` (barcode-to-strain map). `neut_standard_file` is shared across plates.

**Outputs** (to `results/`): `fates_chart.html`, `count_by_barcode_status.html`.

**Notes:**
- Plate layout assumed: 96-well, column-major (A1, B1, ... H1, A2, ...).
- Well names are matched to strains via two conventions (numbered wells via library `shortname`, strain-named wells via library `strain`); unmatched wells are dropped with a printed warning.
- Compare results to actual titers of each strain in the pool: if a well's own strain failed to rescue (zero titer), its contaminating counts/fractions can look misleadingly large.