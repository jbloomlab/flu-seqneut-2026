# library_pooling

Marimo notebooks for library-pooling analyses. Each notebook runs two ways from
the same file:

- **Interactively** with `marimo edit`
- **Non-interactively** via Snakemake (locally or submitted to SLURM), which
  exports each notebook to a self-contained HTML

This directory is a standalone analysis (see `non-pipeline_analyses/`), not part
of the main seqneut-pipeline.

## Environment

All commands run in the `library_pooling` conda environment. Build and activate
it once:

```bash
conda env create -f environment.yml
conda activate library_pooling
```

The environment is pinned in `environment.yml`; add analysis dependencies there
(the single source of truth) rather than installing ad hoc. Snakemake also uses
this file to provision the per-rule environment when run with `--use-conda`.

## Layout

```
library_pooling/
├── config.yaml                 # notebooks + their inputs/outputs/params (single source of truth)
├── environment.yml             # conda environment (name: library_pooling)
├── run_snakemake.bash          # SLURM controller (sbatch this)
├── data/                       # analysis-specific input data (committed)
├── notebooks/                  # marimo notebooks, one .py per analysis
│   └── example_analysis.py     # template — copy to start a new analysis
├── profiles/config.yaml        # Snakemake SLURM profile
├── results/                    # exported HTML + declared outputs
└── workflow/
    ├── Snakefile               # one rule per notebook in config
    └── scripts/
        └── run_marimo_w_context_pickle.py   # driver: pickles context, exports HTML
```

## Running

From this directory:

```bash
# Interactive editing (opens the notebook in a browser)
marimo edit notebooks/example_analysis.py

# Run all notebooks locally
snakemake -s workflow/Snakefile --use-conda --cores 1

# Run a single notebook
snakemake -s workflow/Snakefile --use-conda --cores 1 results/example_analysis.html

# Submit to SLURM (dispatches each notebook as its own job)
sbatch run_snakemake.bash
```

Each run writes `results/<notebook>.html` (the rendered notebook) plus any files
the notebook declares under `output` in `config.yaml`.

## Adding a notebook

1. Copy `notebooks/example_analysis.py` to `notebooks/<your_analysis>.py`.
2. Add an entry under `notebooks` in `config.yaml` with its `input`, `output`,
   and `params`. Reference files by these **keys** inside the notebook (via the
   context), not by hard-coded paths, so the interactive and Snakemake run modes
   stay in sync.
3. Run it with one of the commands above.

## How the two run modes share one file

The Snakemake driver (`workflow/scripts/run_marimo_w_context_pickle.py`) pickles
the rule's inputs/outputs/params into `results/<notebook>_context.pickle` and
runs `marimo export html ... -- --context-pickle <pickle>`. The notebook's
context cell detects the `--context-pickle` argument:

- **present** → running under Snakemake; load that pickle and verify the working
  directory.
- **absent** → running in `marimo edit`; fall back to an optional development
  pickle (set `context_pickle_path` in the notebook's first cell to a real pickle
  after one Snakemake run) or a stub context so the notebook still opens.

The config's `input` / `output` / `params` mappings are surfaced as
`context["input"]` / `context["output"]` / `context["params"]`, so notebooks read
e.g. `context["input"]["designed_library"]` regardless of run mode (matching the
seqneut-pipeline notebooks).
