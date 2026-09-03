# Key H1 HA mutations on the protein structure

Configuration for one [`prot-struct-viz`](https://github.com/jbloomlab/prot-struct-viz)
page: influenza A/Victoria/2570/2019 (H1N1)pdm09 hemagglutinin, uncleaved, drawn as the
deposited trimer, in three views:

- the sites that differ between subclades D.3.1 and D.3.1.1, with sites 155 and 157 of
  antigenic region Sa also picked out;
- the same differing sites alone;
- the classically defined antigenic regions of H1 HA.

`spec.yaml` is the whole input. The CSVs it names color and label the residues, and the
`.md` files are the captions shown beneath each view. Paths inside the spec resolve
relative to the spec file, so this directory can be moved or copied intact. It has no
`out` key: `rule prot_struct_viz_page` passes `--out`, because where the page is written
belongs to the workflow that owns `results/`.

The coordinates are not stored here. `spec.yaml` names PDB entry `9GSP`, which
`prot-struct-viz` fetches from `https://files.rcsb.org/download/9GSP.cif` at render time,
so the entry ID in the spec is the record of what is drawn. That entry is an HA0 whose
protease loop is disordered, which is why those residues are absent from every view.

## Regenerating the CSVs

The CSVs are generated, committed, and read as-is by `spec.yaml`; no rule runs the
generator. Regenerate them by hand from the repository root when its inputs change:

```bash
python data/prot-struct-viz_config/H1/make_coloring_csv.py
```

It reads `data/nextstrain-prot-titers-tree_data/H1N1_site_numbering_map.tsv`,
`data/viral_libraries/flu-seqneut-2026-barcode-to-strain-actual.csv`, `config.yml` and PDB
9GSP, and needs `gemmi`. Its own docstring says what each CSV holds and where the
antigenic-region definitions come from.

The differing sites are not transcribed: the generator picks out the strains whose
`derived_haplotype` is exactly `D.3.1` and `D.3.1.1`, and computes them by comparing the
library's `protein_sequence_HA_ectodomain` column. Correcting a sequence in the library
and re-running is therefore all it takes to correct the page.

The two sites picked out on top of them are the exception: which sites are of interest is
a choice, so they are named in the generator. What the page then says about them is
derived the same way, from every strain the library carries on the D.3.1.1 background. A
named site no such strain mutates, or one outside the antigenic region whose color they
share, fails the generator.
