# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.23.0"
app = marimo.App(width="full")


@app.cell
def _():
    # Load context from pickled file.
    #
    # This cell supports multiple ways to provide context:
    # 1. Via command-line: marimo export html notebook.py -- --context-pickle path/to/context.pickle
    # 2. Via saved pickle: Manually save a context pickle to results/context_dev.pickle
    # 3. Stub context: If no pickle available, creates minimal empty context for exploration
    #
    # For interactive development with `marimo edit`, you can:
    # - Run the pipeline once to generate a real context pickle, then copy it to context_dev.pickle
    # - Or work with the stub context (downstream cells will show warnings/empty data)

    import argparse
    import os
    import pathlib
    import pickle
    import sys

    import marimo as mo

    # Check if context-pickle argument is provided (run by driver script)
    from_cmdline = "--context-pickle" in sys.argv

    if from_cmdline:
        # Running via driver script - parse args
        print("Loading context from command-line argument")
        p = argparse.ArgumentParser()
        p.add_argument("--context-pickle", required=True)
        args = p.parse_args()
        context_pickle_path = pathlib.Path(args.context_pickle)
    else:
        # Running in marimo edit - try to use development pickle
        print("Running in marimo edit mode")
        # if running in edit mode, set `context_pickle_path` to valid pickle
        context_pickle_path = None
        # context_pickle_path = pathlib.Path("results/example_analysis_context.pickle")

    # Load context if pickle path exists and is valid
    if context_pickle_path and context_pickle_path.exists():
        print(f"Reading context from {context_pickle_path}")
        with open(context_pickle_path, "rb") as f_context:
            context = pickle.load(f_context)

        # Handle working directory
        context_workdir = context["workdir"]
        current_workdir = os.getcwd()

        if from_cmdline:
            # Running via snakemake - verify workdir matches
            if context_workdir != current_workdir:
                raise RuntimeError(
                    f"Context workdir mismatch!\n"
                    f"  Context was created in: {context_workdir}\n"
                    f"  Currently running in:   {current_workdir}\n"
                    f"This should not happen when running via Snakemake."
                )
            print(f"Verified working directory: {current_workdir}")
        else:
            # Running in marimo edit - change to context workdir
            if context_workdir and context_workdir != current_workdir:
                print(f"Changing directory from {current_workdir} to {context_workdir}")
                os.chdir(context_workdir)
            elif context_workdir:
                print(f"Already in correct working directory: {context_workdir}")
    else:
        # Create a minimal stub context for interactive development
        print("Creating minimal stub context that you need to complete")
        context = {
            "input": {},
            "output": {},
            "params": {},
            "wildcards": {},
            "threads": 1,
            "resources": {},
        }
    return context, mo


@app.cell
def _(context, mo):
    # Extract variables from context - raises KeyError if required keys missing.
    # The Snakefile wires the config `input`/`output`/`params` mappings into the
    # rule's native input/output/params, so they are read by key here.
    stub_context = not context["input"]

    if not stub_context:
        input_designed_library = context["input"]["designed_library"]
        output_summary_csv = context["output"]["summary_csv"]
        min_barcodes_per_strain = context["params"]["min_barcodes_per_strain"]
    else:
        # Interactive stub: fill in a real designed-library path to explore.
        input_designed_library = None
        output_summary_csv = None
        min_barcodes_per_strain = 2

    # Show informative message about context mode
    if stub_context:
        mo.output.append(
            mo.callout(
                mo.md(
                    "**⚠️ Running in interactive mode with stub context**\n\n"
                    "To run with real data:\n"
                    "1. Run the pipeline to generate a context pickle\n"
                    "2. Copy it to `results/context_dev.pickle` and point "
                    "`context_pickle_path` at it\n"
                    "3. Or run: `marimo export html notebook.py -- --context-pickle path/to/context.pickle`\n\n"
                    "Meanwhile, set `input_designed_library` below to a real CSV to explore."
                ),
                kind="warn",
            )
        )
    return (
        input_designed_library,
        min_barcodes_per_strain,
        output_summary_csv,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Example library_pooling analysis
    Template marimo notebook for the `library_pooling` analyses. Runs both
    interactively (`marimo edit`) and via Snakemake (locally or on SLURM). Copy it
    to start a new analysis, adding an entry under `notebooks` in `config.yaml` and
    referencing files by their config keys.
    """)
    return


@app.cell
def _():
    # Setup and read data

    import altair as alt

    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()
    return alt, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load the designed library
    Read the barcode-to-strain designed library CSV named in the config.
    """)
    return


@app.cell
def _(input_designed_library, mo, pd):
    if input_designed_library is None:
        mo.stop(
            True,
            mo.callout(
                mo.md(
                    "No `input_designed_library` set; complete the stub context above."
                ),
                kind="warn",
            ),
        )

    df = pd.read_csv(input_designed_library)
    mo.output.append(
        mo.md(f"Loaded **{len(df)}** barcode rows from `{input_designed_library}`.")
    )
    df
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summarize barcodes and strains per subtype
    Count strains and barcodes per subtype, flagging strains below the configured
    minimum barcode count.
    """)
    return


@app.cell
def _(df, min_barcodes_per_strain, mo, output_summary_csv):
    import pathlib as _pathlib

    per_strain = (
        df.groupby(["subtype", "strain"]).size().rename("n_barcodes").reset_index()
    )
    summary = (
        per_strain.groupby("subtype")
        .agg(
            n_strains=("strain", "nunique"),
            n_barcodes=("n_barcodes", "sum"),
            n_strains_below_min=(
                "n_barcodes",
                lambda s: int((s < min_barcodes_per_strain).sum()),
            ),
        )
        .reset_index()
    )

    if output_summary_csv:
        out_path = _pathlib.Path(output_summary_csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_path, index=False)
        mo.output.append(mo.md(f"Wrote summary to `{out_path}`."))

    mo.output.append(
        mo.md(
            f"Strains with fewer than **{min_barcodes_per_strain}** barcodes are flagged."
        )
    )
    summary
    return (per_strain,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plot the distribution of barcodes per strain
    """)
    return


@app.cell
def _(alt, mo, per_strain):
    chart = (
        alt.Chart(per_strain)
        .mark_bar()
        .encode(
            alt.X(
                "n_barcodes:Q",
                bin=alt.Bin(maxbins=20),
                title="barcodes per strain",
            ),
            alt.Y("count()", title="number of strains"),
            color=alt.Color("subtype:N", title="subtype"),
            tooltip=["subtype:N", alt.Tooltip("count()", title="strains")],
        )
        .properties(
            width=500,
            title="Distribution of barcodes per strain",
        )
    )
    mo.output.append(mo.ui.altair_chart(chart))
    return


if __name__ == "__main__":
    app.run()
