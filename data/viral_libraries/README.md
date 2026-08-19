# Viral libraries data

The *designed* library lists every barcode designed into the library.
It is built by [create-designed-library-csv.py](create-designed-library-csv.py) from the final library of the [library design analysis](../../non-pipeline_analyses/library_design/construct_order).

That analysis supplies *subclade* and *derived_haplotype* from hand-curated haplotype tables, which left some strains labelled with only a bare subclade.
[assign-haplotypes-nextclade.py](assign-haplotypes-nextclade.py) therefore re-assigns both columns in the designed library from `nt_sequence_HA_ectodomain` using `nextclade`, and prints what it changed.
It always uses the latest `nextclade` dataset, recording which one in [assign-haplotypes-nextclade_provenance.txt](assign-haplotypes-nextclade_provenance.txt).

Some strains and barcodes were then determined experimentally not to be usable, having failed rescue or given too low a titer.
These are listed with the reason for each in [drop_strains.csv](drop_strains.csv).
The *actual* library is the designed library with those rows removed, and is built by [create-actual-library-csv.py](create-actual-library-csv.py).
It is the library used for the titer measurements.

The scripts are run by hand from this directory when their inputs change, in the order they are listed below; the pipeline reads only the CSVs.
`assign-haplotypes-nextclade.py` needs `nextclade`, so it is run in the environment defined by [assign-haplotypes-nextclade_environment.yml](assign-haplotypes-nextclade_environment.yml); the other two need only the pipeline environment.
Which library each analysis uses is set in [config.yml](../../config.yml).

## Files

- [create-designed-library-csv.py](create-designed-library-csv.py) — builds the designed library
- [assign-haplotypes-nextclade.py](assign-haplotypes-nextclade.py) — re-assigns *subclade* and *derived_haplotype* in the designed library
- [assign-haplotypes-nextclade_environment.yml](assign-haplotypes-nextclade_environment.yml) — conda environment for that script
- [assign-haplotypes-nextclade_provenance.txt](assign-haplotypes-nextclade_provenance.txt) — `nextclade` version and dataset tags used for the current assignment
- [flu-seqneut-2026-barcode-to-strain-designed.csv](flu-seqneut-2026-barcode-to-strain-designed.csv) — designed library, every barcode
- [drop_strains.csv](drop_strains.csv) — dropped strain/barcode pairs, with a note on why
- [create-actual-library-csv.py](create-actual-library-csv.py) — builds the actual library
- [flu-seqneut-2026-barcode-to-strain-actual.csv](flu-seqneut-2026-barcode-to-strain-actual.csv) — actual library, used for titers
