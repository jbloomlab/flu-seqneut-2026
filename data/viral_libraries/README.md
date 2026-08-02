# Viral libraries data

The *designed* library lists every barcode designed into the library.
It is built by [create-designed-library-csv.py](create-designed-library-csv.py) from the final library of the [library design analysis](../../non-pipeline_analyses/library_design/construct_order).

Some strains and barcodes were then determined experimentally not to be usable, having failed rescue or given too low a titer.
These are listed with the reason for each in [drop_strains.csv](drop_strains.csv).
The *actual* library is the designed library with those rows removed, and is built by [create-actual-library-csv.py](create-actual-library-csv.py).
It is the library used for the titer measurements.

Both scripts are run by hand from this directory when their inputs change; the pipeline reads only the CSVs.
Which library each analysis uses is set in [config.yml](../../config.yml).

## Files

- [create-designed-library-csv.py](create-designed-library-csv.py) — builds the designed library
- [flu-seqneut-2026-barcode-to-strain-designed.csv](flu-seqneut-2026-barcode-to-strain-designed.csv) — designed library, every barcode
- [drop_strains.csv](drop_strains.csv) — dropped strain/barcode pairs, with a note on why
- [create-actual-library-csv.py](create-actual-library-csv.py) — builds the actual library
- [flu-seqneut-2026-barcode-to-strain-actual.csv](flu-seqneut-2026-barcode-to-strain-actual.csv) — actual library, used for titers
