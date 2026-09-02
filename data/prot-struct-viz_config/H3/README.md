# Key H3 HA mutations on the protein structure

Configuration for one [`prot-struct-viz`](https://github.com/jbloomlab/prot-struct-viz)
page: influenza A/Victoria/22/2020 (H3N2) hemagglutinin with the LSTc receptor analogue
bound, drawn as the biological trimer, in three views:

- the sites that differ between the 2025-2026 and 2026-2027 vaccine strains, with sites
  222 and 223 of antigenic region D also picked out;
- the same differing sites alone;
- the classically defined antigenic regions of H3 HA.

`spec.yaml` is the whole input. The CSVs it names color and label the residues, and the
`.md` files are the captions shown beneath each view. Paths inside the spec resolve
relative to the spec file, so this directory can be moved or copied intact. It has no
`out` key: `rule prot_struct_viz_page` passes `--out`, because where the page is written
belongs to the workflow that owns `results/`.

The coordinates are not stored here. `spec.yaml` names PDB entry `8FAW`, which
`prot-struct-viz` fetches from `https://files.rcsb.org/download/8FAW.cif` at render time,
so the entry ID in the spec is the record of what is drawn.

## Regenerating the CSVs

The three CSVs are generated, committed, and read as-is by `spec.yaml`; no rule runs the
generator. Regenerate them by hand from the repository root when its inputs change:

```bash
python data/prot-struct-viz_config/H3/make_coloring_csv.py
```

It reads `data/nextstrain-prot-titers-tree_data/H3N2_site_numbering_map.tsv` and PDB
8FAW, and needs `gemmi`. Its own docstring says what each CSV holds and where the
antigenic-region definitions come from.

Its list of sites differing between the two vaccine strains is a transcribed constant
rather than being recomputed from `data/viral_libraries/`. Deriving it from the libraries
would remove the transcription step; until then the script's assertions on list length,
site-to-residue mapping, and whether each site is modeled in 8FAW are what guard it.
