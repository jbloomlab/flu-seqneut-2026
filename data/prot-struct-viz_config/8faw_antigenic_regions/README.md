# H3 HA antigenic regions and recent substitutions

Configuration for one [`prot-struct-viz`](https://github.com/jbloomlab/prot-struct-viz)
page: influenza A/Victoria/22/2020 (H3N2) hemagglutinin with the LSTc receptor analogue
bound, drawn as the biological trimer, in five views.

`spec.yaml` is the whole input. The CSVs it names color and label the residues, and the
`.md` files are the captions shown beneath each view. Paths inside the spec resolve
relative to the spec file, so this directory can be moved or copied intact. It has no
`out` key: `rule prot_struct_viz_page` passes `--out`, because where the page is written
belongs to the workflow that owns `results/`.

## Provenance

Copied from `examples/8faw_antigenic_regions` of `prot-struct-viz` v0.2.0 (commit
`30ba448be8474d8452bb7b4576adb9537e049b0f`) on 2026-09-02, with three changes: the
`out` key was dropped from `spec.yaml`; `make_coloring_csv.py` reads this repository's
numbering map from disk instead of fetching it over HTTP; and four implicit string
concatenations in that script were parenthesized, which this repository's `ruff`
configuration requires and upstream's does not. The CSVs and captions are verbatim.

The coordinates are not stored here. `spec.yaml` names PDB entry `8FAW`, which
`prot-struct-viz` fetches from `https://files.rcsb.org/download/8FAW.cif` at render time,
so the entry ID in the spec is the record of what is drawn.

## Regenerating the CSVs

The five CSVs are generated, committed, and read as-is by `spec.yaml`; no rule runs the
generator. Regenerate them by hand from the repository root when its inputs change:

```bash
python data/prot-struct-viz_config/8faw_antigenic_regions/make_coloring_csv.py
```

It reads `data/nextstrain-prot-titers-tree_data/H3N2_site_numbering_map.tsv` and PDB
8FAW, and needs `gemmi`. Its own docstring says what each CSV holds and where the
antigenic-region definitions come from.

Its two mutation lists are transcribed constants rather than being recomputed from
`data/viral_libraries/`, inherited from the upstream example where those libraries were
not available. Deriving them from the libraries instead would remove the transcription
step; until then the script's assertions on list length, site-to-residue mapping, and
whether each site is modeled in 8FAW are what guard them.
