# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full")


@app.cell(hide_code=True)
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

    from pathlib import Path
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
    return Path, context, mo, os


@app.cell(hide_code=True)
def _(context, mo):
    # Extract variables from context - raises KeyError if required keys missing.
    # The Snakefile wires the config `input`/`output`/`params` mappings into the
    # rule's native input/output/params, so they are read by key here.
    stub_context = not context["input"]

    if not stub_context:
        viral_library_csv = context["input"]["viral_library_csv"]
        neut_standard_set_csv = context["input"]["neut_standard_set"]
        samplesfile = context["input"]["samplesfile"]
        platedir = context["input"]["platedir"]
        repooling_math = context["input"]["repooling_math"]
    else:
        # Interactive stub: fill in parameters.
        viral_library_csv = '../../data/viral_libraries/flu-seqneut-2026-barcode-to-strain-actual.csv'
        neut_standard_set_csv = '../../data/neut_standard_sets/loes2023_neut_standards.csv'
        samplesfile = '../../data/miscellaneous_plates/2026-07-31_repool.csv'
        platedir = '../../results/miscellaneous_plates/20260731_repool/'
        repooling_math = './../results/pooling_math/2026-07-15_repooling_math.csv'

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
                    "Meanwhile, set variables above to explore."
                ),
                kind="warn",
            )
        )
    return (
        neut_standard_set_csv,
        platedir,
        repooling_math,
        samplesfile,
        viral_library_csv,
    )


@app.cell
def _():
    import altair as alt

    import numpy as np
    import pandas as pd

    _ = alt.data_transformers.disable_max_rows()

    # Basic color palette
    color_palette = [
        '#345995', #blue
        '#03cea4', #teal
        '#ca1551', #red
        '#eac435', #yellow
                   ]
    return alt, color_palette, np, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Add input data locations
    Some of these files are defined as data, and some of these files are generated by running the specified library pooling data as `miscellaneous_plates` through the `seqneut-pipeline`. For details on how these files are generated, see the `README.md' in [https://github.com/jbloomlab/seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline)
    """)
    return


@app.cell
def _(os, platedir):
    # Identify all counts and fates CSVs
    file_list = os.listdir(platedir)
    count_csvs = [os.path.join(platedir, f) for f in file_list if "_counts" in f]
    fate_csvs = [os.path.join(platedir, f) for f in file_list if "_fates" in f]
    return count_csvs, fate_csvs


@app.cell
def _(pd, samplesfile):
    # Define a samples dataframe using the samples file
    samples_df = pd.read_csv(samplesfile)
    samples_df = samples_df.drop(columns=['fastq'])
    samples_df['sample'] = samples_df.apply(
        lambda x: '-'.join(x.astype(str)), axis=1
    )

    samples = samples_df["sample"].unique().tolist()
    print(f"There are {len(samples)} barcode runs.")

    samples_df
    return samples, samples_df


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Statistics on barcode-parsing for each sample
    Make interactive chart of the "fates" of the sequencing reads parsed for each sample on the plate.

    If most sequencing reads are not "valid barcodes", this could potentially indicate some problem in the sequencing or barcode set you are parsing.

    Potential fates are:

    - **valid barcode**: barcode that matches a known virus or neutralization standard, we hope most reads are this.
    -  **invalid barcode**: a barcode with proper flanking sequences, but does not match a known virus or neutralization standard. If you  have a lot of reads of this type, it is probably a good idea to look at the invalid barcode CSVs (in the `./results/barcode_invalid/` subdirectory created by the pipeline) to see what these invalid barcodes are.
    -  **unparseable barcode**: could not parse a barcode from this read as there was not a sequence of the correct length with the appropriate flanking sequence.
    - **low quality barcode**: low-quality or `N` nucleotides in barcode, could indicate problem with sequencing.
    - **failed chastity filter**: reads that failed the Illumina chastity filter, if these are reported in the FASTQ (they may not be).

    Also, if the number of reads per sample is very uneven, that could indicate that you did not do a good job of balancing the different samples in the Illumina sequencing.
    """)
    return


@app.cell
def _(Path, alt, fate_csvs, pd, samples, samples_df):
    assert len(fate_csvs) == len(samples)

    fates = (
        pd.concat([
            pd.read_csv(f).assign(well=Path(f).stem.removesuffix('_fates'))
            for f in fate_csvs
        ])
        .merge(samples_df, validate="many_to_one", on="well")
        .assign(
            fate_counts=lambda x: x.groupby("fate")["count"].transform("sum"),
            sample_well=lambda x: x["sample"] + " (" + x["well"] + ")",
        )
        .query("fate_counts > 0")[  # only keep fates with at least one count
            ["fate", "count", "well", "sample_well", "dilution_factor"]
        ]
    )

    assert len(fates) == len(fates.drop_duplicates())

    sample_wells = list(
        fates.sort_values(["dilution_factor"])["sample_well"]
    )

    fates_chart = (
        alt.Chart(fates)
        .encode(
            alt.X("count", scale=alt.Scale(nice=False, padding=3)),
            alt.Y(
                "sample_well",
                title=None,
                sort=sample_wells,
            ),
            alt.Color("fate", sort=sorted(fates["fate"].unique(), reverse=True)),
            alt.Order("fate", sort="descending"),
            tooltip=fates.columns.tolist(),
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=200,
            title=f"Barcode parsing for initial titering plate",
        )
        .configure_axis(grid=False)
    )

    fates_chart
    return (sample_wells,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Read barcode counts
    Read the counts per barcode:
    """)
    return


@app.cell
def _(
    Path,
    count_csvs,
    neut_standard_set_csv,
    pd,
    sample_wells,
    samples,
    samples_df,
    viral_library_csv,
):
    # get barcode counts
    assert len(count_csvs) == len(samples)
    counts = (
        pd.concat([
            pd.read_csv(c).assign(well=Path(c).stem.removesuffix('_counts'))
            for c in count_csvs
        ])
        .merge(samples_df, validate="many_to_one", on="well")
        .drop(columns=["replicate"])
        .assign(sample_well=lambda x: x["sample"] + " (" + x["well"] + ")")
    )


    # classify barcodes as viral or neut standard
    barcode_class = pd.concat(
        [
            pd.read_csv(viral_library_csv)[["barcode", "strain"]].assign(
                neut_standard=False,
            ),
            pd.read_csv(neut_standard_set_csv)[["barcode"]].assign(
                neut_standard=True,
                strain=pd.NA,
            ),
        ],
        ignore_index=True,
    )

    # merge counts and classification of barcodes
    assert set(counts["barcode"]) == set(barcode_class["barcode"])
    counts = counts.merge(barcode_class, on="barcode", validate="many_to_one")
    assert set(sample_wells) == set(counts["sample_well"])
    return (counts,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Average counts per barcode in each well

    Plot average counts per barcode.
    If a sample has inadequate barcode counts, it may not have good enough statistics for accurate analysis, and a QC-threshold is applied:
    """)
    return


@app.cell
def _(alt, counts, pd, sample_wells):
    MIN_AVG_BARCODE_COUNT = 500
    avg_barcode_counts = (
        counts.groupby(
            ["well", "sample_well"],
            dropna=False,
            as_index=False,
        )
        .aggregate(avg_count=pd.NamedAgg("count", "mean"))
        .assign(
            fails_qc=lambda x: (
                x["avg_count"] < MIN_AVG_BARCODE_COUNT
            ),
        )
    )

    avg_barcode_counts_chart = (
        alt.Chart(avg_barcode_counts)
        .encode(
            alt.X(
                "avg_count",
                title="average barcode counts per well",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Color(
                "fails_qc",
                title=f"fails {MIN_AVG_BARCODE_COUNT=}",
                legend=alt.Legend(titleLimit=500),
            ),
            tooltip=[
                alt.Tooltip(c, format=".3g") if pd.api.types.is_float_dtype(avg_barcode_counts[c]) else c
                for c in avg_barcode_counts.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=250,
            title=f"Average barcode counts per well for titering plate",
        )
        .configure_axis(grid=False)
    )

    avg_barcode_counts_chart
    return (avg_barcode_counts,)


@app.cell
def _(avg_barcode_counts):
    # wells failing QC (checked against the wells chosen for balancing, below)
    avg_barcode_counts_per_well_drops = list(avg_barcode_counts.query("fails_qc")["well"])
    print(f"wells failing {avg_barcode_counts_per_well_drops=}")
    return (avg_barcode_counts_per_well_drops,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Fraction of counts from neutralization standard
    Determine the fraction of counts from the neutralization standard in each sample, and make sure this fraction passess the QC threshold.
    """)
    return


@app.cell
def _(alt, counts, np, pd, sample_wells):
    MIN_NEUT_STANDARD_FRAC = 0.001
    neut_standard_fracs = (
        counts.assign(
            neut_standard_count=lambda x: x["count"] * x["neut_standard"].astype(int)
        )
        .groupby(
            ["well", "sample_well", 'dilution_factor'],
            dropna=False,
            as_index=False,
        )
        .aggregate(
            total_count=pd.NamedAgg("count", "sum"),
            neut_standard_count=pd.NamedAgg("neut_standard_count", "sum"),
        )
        .assign(
            neut_standard_frac=lambda x: x["neut_standard_count"] / x["total_count"],
            fails_qc=lambda x: (
                x["neut_standard_frac"] < MIN_NEUT_STANDARD_FRAC
            ),
            neut_standard_logit=lambda x: np.log(
                x["neut_standard_frac"] / (1 - x["neut_standard_frac"])
            )
        )
    )

    neut_standard_fracs_chart = (
        alt.Chart(neut_standard_fracs)
        .encode(
            alt.X(
                "neut_standard_frac",
                title="frac counts from neutralization standard per well",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("sample_well", sort=sample_wells),
            alt.Color(
                "fails_qc",
                title=f"fails {MIN_NEUT_STANDARD_FRAC=}",
                legend=alt.Legend(titleLimit=500),
            ),
            tooltip=[
                alt.Tooltip(c, format=".3g") if neut_standard_fracs[c].dtype == float else c
                for c in neut_standard_fracs.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=250,
            title=f"Neutralization-standard fracs per well for titering plate, initial pool",
        )
        .configure_axis(grid=False)
        .configure_legend(titleLimit=1000)
    )

    neut_standard_fracs_chart
    return (neut_standard_fracs,)


@app.cell
def _(neut_standard_fracs):
    # wells failing QC (checked against the wells chosen for balancing, below)
    min_neut_standard_frac_per_well_drops = list(
        neut_standard_fracs.query("fails_qc")["well"]
    )
    print(f"wells failing {min_neut_standard_frac_per_well_drops=}")
    neut_standard_fracs
    return (min_neut_standard_frac_per_well_drops,)


@app.cell
def _(alt, neut_standard_fracs):
    # Scatterplot of the same data as above, plotted by dilution factor
    alt.Chart(neut_standard_fracs).mark_circle(size=60).encode(
        alt.X('dilution_factor:Q',
              scale=alt.Scale(type='log'),
              title='library pool reciprocal dilution factor'),
        alt.Y('neut_standard_frac:Q',
              title='fraction of reads = neutralization standard'),
        color='fails_qc',
        tooltip=['well', 'dilution_factor', 'neut_standard_frac', 'total_count']
    ).interactive()
    return


@app.cell
def _(alt, neut_standard_fracs):
    # Scatterplot with the logit instead, plotted by dilution factor
    alt.Chart(neut_standard_fracs).mark_circle(size=60).encode(
        alt.X('dilution_factor:Q', 
              scale=alt.Scale(type='log'),
              title='library pool reciprocal dilution factor'),
        alt.Y('neut_standard_logit:Q', 
              title='logit of fraction of reads = neutralization standard'),
        color='fails_qc',
        tooltip=['well', 'dilution_factor', 'neut_standard_frac', 'total_count']
    ).interactive()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Assess balancing of strains contained in the library
    """)
    return


@app.cell
def _(counts, np):
    # Strain balancing concerns the viral strains only, so drop the neutralization
    # standard explicitly rather than relying on its NaN `strain` being discarded
    # by the groupby's default dropna=True (the other groupbys in this notebook
    # pass dropna=False, so that default is easy to flip by accident).
    viral_counts = counts[~counts['neut_standard']]
    assert viral_counts['strain'].notna().all(), "viral barcodes with missing strain"
    assert counts.loc[counts['neut_standard'], 'strain'].isna().all(), (
        "neut-standard barcodes unexpectedly carry a strain"
    )

    # Get summed barcode counts for all strains across all wells
    straincounts_allbarcodes = (
        viral_counts.groupby(
            ['sample', 'sample_well', 'strain', 'dilution_factor', 'serum', 'well'],
            dropna=False,
        )['count']
        .sum()
        .reset_index()
    )

    # Get sum of all virus/barcode counts per well
    sumperwell = (
        straincounts_allbarcodes.groupby(
            ['sample', 'dilution_factor', 'serum', 'well'], dropna=False
        )['count']
        .sum()
        .reset_index()
        .rename(columns={'count': 'counts_perwell'})
    )

    # Merge dataframes and calculate fraction of each well devoted to each strain
    merged_df = straincounts_allbarcodes.merge(sumperwell, on=['sample','dilution_factor','serum','well'])
    # Count how many barcodes were used for each strain (this can vary, e.g. 1-4)
    barcodes_per_strain = (
        viral_counts[['strain', 'barcode']]
        .drop_duplicates()
        .groupby('strain')
        .size()
        .rename('n_barcodes')
        .reset_index()
    )
    print(barcodes_per_strain['n_barcodes'].value_counts())
    merged_df = merged_df.merge(barcodes_per_strain, on='strain', validate='many_to_one')
    # Fraction of the well's viral counts belonging to this strain. `count` is
    # already summed over the strain's barcodes, so this must NOT be divided by
    # n_barcodes -- doing so rescales strains unequally (n_barcodes varies 1-3)
    # and makes the total sum to <1, which breaks the comparison against the
    # absolute 1/num_strains expectation used to call over-representation.
    merged_df['fraction_strain'] = merged_df['count'] / merged_df['counts_perwell']

    # Fractions are over viral strains only, so they must sum to 1 per well.
    _frac_sums = merged_df.groupby('well')['fraction_strain'].sum()
    assert np.allclose(_frac_sums, 1), f"fraction_strain does not sum to 1:\n{_frac_sums}"

    merged_df
    return (merged_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We now have this fraction of reads devoted to all strains calculated for all wells. However, ideally we should just focus on those wells containing dilutions that we would use for actual neutralization assays. We should choose a set of replicate wells where the fraction of neutralization standard reads begins to increase linearly with the increasing reciprocal dilution factor. See plots above for choosing these wells.
    """)
    return


@app.cell
def _(
    avg_barcode_counts_per_well_drops,
    merged_df,
    min_neut_standard_frac_per_well_drops,
):
    # Choose a pair of replicate wells near the beginning of the linear range
    chosen_wells = ['A5', 'B5']

    # Match the well column exactly rather than substring-matching the composite
    # `sample` string, where e.g. 'A1' would also match 'A10'/'A11'/'A12'.
    missing_wells = set(chosen_wells) - set(merged_df['well'])
    assert not missing_wells, f"chosen wells not in data: {sorted(missing_wells)}"

    # The wells above are picked by eye off the QC plots, so enforce that they
    # actually passed the QC thresholds applied earlier in the notebook.
    failed_wells = set(chosen_wells) & (
        set(avg_barcode_counts_per_well_drops) | set(min_neut_standard_frac_per_well_drops)
    )
    assert not failed_wells, f"chosen wells failed QC: {sorted(failed_wells)}"

    single_well = merged_df.loc[merged_df['well'].isin(chosen_wells)]
    return (chosen_wells, single_well)


@app.cell
def _(np, single_well):
    # Calculate mean fraction strain across both wells
    mean_df = (
        single_well.groupby('strain')['fraction_strain']
        .mean()
        .rename('mean_fraction_strains')
        .reset_index()
    )
    mean_single_well = single_well.merge(mean_df, on='strain', how='left')

    # calcualte ratios to add for equal pool
    num_strains = mean_single_well['strain'].nunique()
    mean_single_well['ratio_to_add'] = (1/num_strains)/mean_single_well['fraction_strain']
    mean_single_well['mean_ratio_to_add'] = (1/num_strains)/mean_single_well['mean_fraction_strains']

    # Unlike the initial equal-volume pool, every strain in this repool is present
    # at non-zero frequency, so there are no infinite ratios. Assert that rather
    # than carrying the equal-volume notebook's zero-titer filtering, which would
    # silently do nothing here.
    assert np.isfinite(mean_single_well['mean_ratio_to_add']).all(), (
        "strains absent from the pool:\n"
        f"{mean_single_well.loc[~np.isfinite(mean_single_well['mean_ratio_to_add']), 'strain'].unique()}"
    )

    print(f'this library has {num_strains} total strains')
    print(mean_single_well[['mean_ratio_to_add']].describe())
    return mean_single_well, num_strains


@app.cell
def _(mean_single_well):
    mean_single_well.sort_values("mean_fraction_strains", ascending = False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Visualize barcode- and strain-level balancing in the current pool
    """)
    return


@app.cell
def _(alt, mean_single_well, num_strains, pd):
    # One row per strain for the bars. `mean_single_well` carries a row per
    # (strain, well), and a bar mark with two rows in the same y-category stacks
    # them -- which would draw every strain at twice its true fraction, against a
    # reference line drawn at the single-well expectation.
    strain_means = mean_single_well[['strain', 'mean_fraction_strains']].drop_duplicates()
    assert len(strain_means) == num_strains

    # Mean fraction across the chosen replicate wells.
    strains_chart = (
        alt.Chart(strain_means)
        .mark_bar(height={"band": 0.85})
        .encode(
            alt.X(
                "mean_fraction_strains",
                title="fraction of pool",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("strain"),
            tooltip=['strain', 'mean_fraction_strains'],
        )
    )

    # Individual wells overlaid, so replicate disagreement stays visible rather
    # than being averaged away.
    replicate_points = (
        alt.Chart(mean_single_well)
        .mark_point(size=15, filled=True, opacity=0.9, color="black")
        .encode(
            alt.X("fraction_strain"),
            alt.Y("strain"),
            tooltip=['strain', 'well', 'fraction_strain'],
        )
    )

    # add vertical line where we would expect equal representation of all strains in pool
    expected_line = alt.Chart(
        pd.DataFrame({'x': [1/num_strains]})
    ).mark_rule(strokeDash=[2, 2], strokeWidth=2).encode(x='x')

    (
        (strains_chart + replicate_points + expected_line)
        .properties(
            height=alt.Step(10),
            width=200,
            title="Strain representation, repool (bars = mean, points = individual wells)",
        )
        .configure_axis(grid=False)
    )
    return


@app.cell
def _(alt, chosen_wells, color_palette, counts):
    # Each barcode fraction across strains. The dropna() drops the neutralization
    # standard, whose `strain` is NA, leaving viral barcodes only.
    all_barcode_counts = counts[['strain', 'barcode', 'count', 'well']].dropna()
    single_well_all_barcode_counts = all_barcode_counts[all_barcode_counts['well'].isin(chosen_wells)]

    # Get tidy single well means
    tidy_single_well = single_well_all_barcode_counts[['strain','barcode','count']].groupby(['strain', 'barcode']).mean().reset_index()
    # Get sums for each strain. Select `count` explicitly: an unrestricted .sum()
    # also "sums" the barcode column, concatenating the sequence strings.
    strain_sums_df = (
        tidy_single_well.groupby('strain')['count']
        .sum()
        .rename('strain_count_sum')
        .reset_index()
    )
    # Merge and calculate per strain the fraction represented by each barcode
    tidy_single_well = tidy_single_well.merge(
        strain_sums_df,
        on='strain',
        validate="many_to_one",
    )
    tidy_single_well['per_strain_fraction_barcode'] = tidy_single_well['count'] / tidy_single_well['strain_count_sum']
    # Purely positional index used to colour the segments of each strain's bar.
    # The rows are ordered by barcode sequence (groupby sorts its keys), so this
    # does NOT correspond to the library's own _bc1/_bc2 designation -- it just
    # separates one strain's barcodes from each other visually. The true barcode
    # is in the tooltip.
    tidy_single_well['barcode_index'] = tidy_single_well.groupby('strain').cumcount()

    # Plot as colored bar chart
    bar_chart = alt.Chart(tidy_single_well).mark_bar(height={"band": 0.85}).encode(
        x = 'per_strain_fraction_barcode',
        y = 'strain',
        color=alt.Color('barcode_index:N', legend=None).scale(range=color_palette),
        tooltip = ['strain', 'per_strain_fraction_barcode', 'barcode'],
    ).configure_axis(grid=False).properties(
            height = alt.Step(10),
            width = 200,
            title = "Barcode fraction for each strain, repool")

    bar_chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Repooling for over-represented strains

    Some strains are over-represented here. It seems like these are systematically the strains that had the highest titer in the initial equal volume pool. We now want to attempt a second repool to better balance these strains. The following are strains where the observed fraction exceeded the expected fraction (`1/num_strains`) by the factor set in `OVER_REP_FACTOR` below.

    Note that a handful of strains are over-represented by a wide margin rather than marginally: the most abundant strain sits at roughly 12x its expected fraction, and the top few strains together account for a substantial share of the whole pool. These dominate any rebalancing regardless of which subpool they fall in.
    """)
    return


@app.cell
def _(mean_single_well, num_strains):
    mean_single_well_sorted = mean_single_well.sort_values('mean_fraction_strains', ascending=False)
    # Call a strain over-represented when its fraction exceeds the equal-pool
    # expectation (1/num_strains) by this factor.
    OVER_REP_FACTOR = 1.25
    over_rep_threshold = OVER_REP_FACTOR / num_strains
    over_rep_strains = mean_single_well_sorted[mean_single_well_sorted['mean_fraction_strains'] > over_rep_threshold]
    over_rep_strains.drop_duplicates(subset='strain')[["strain", "mean_fraction_strains"]].reset_index(drop=True)
    return mean_single_well_sorted, over_rep_strains


@app.cell
def _(mean_single_well_sorted, pd, viral_library_csv):
    # Collapse to one row per strain on both sides before merging. `mean_single_well_sorted`
    # has a row per (strain, well) and the library has a row per (strain, barcode), so
    # merging them directly cross-products to ~4 rows per strain and forces every
    # downstream step to re-deduplicate.
    strain_fracs = (
        mean_single_well_sorted[['strain', 'mean_fraction_strains', 'mean_ratio_to_add', 'n_barcodes']]
        .drop_duplicates(subset='strain')
    )

    # Per-strain library metadata. `shortname` and `bloom_lab_plasmid_log_id` are
    # barcode-level (they carry a trailing _bc1/_bc2), so strip the suffix from
    # shortname to get the strain-level name the subpool rules below key off, and
    # drop the plasmid id, which has no strain-level equivalent.
    lib_id_df = pd.read_csv(viral_library_csv).drop(
        columns=[
            'barcode',
            'bloom_lab_plasmid_log_id',
            'nt_sequence_HA_ectodomain',
            'protein_sequence_HA_ectodomain',
        ]
    )
    lib_id_df['shortname'] = lib_id_df['shortname'].str.replace(r'_bc\d+$', '', regex=True)

    lib_per_strain = lib_id_df.drop_duplicates()
    # One row per strain must now be lossless; if any field still varies within a
    # strain this fails rather than silently keeping an arbitrary barcode's value.
    _varying = lib_per_strain.loc[lib_per_strain['strain'].duplicated(keep=False)]
    assert _varying.empty, (
        f"library metadata varies within a strain:\n{_varying.sort_values('strain')}"
    )

    assert set(strain_fracs['strain']) == set(lib_per_strain['strain']), (
        "strain mismatch between counts data and viral library"
    )
    repool_df = (
        strain_fracs
        .merge(lib_per_strain, on='strain', how='left', validate='one_to_one')
        .reset_index(drop=True)
    )
    assert len(repool_df) == strain_fracs['strain'].nunique()
    return (repool_df,)


@app.cell
def _(np, over_rep_strains, repool_df):
    # Add in subpool information
    conditions = [
        repool_df['shortname'].str.startswith('flu-seqneut-2026_H1N1'),
        repool_df['shortname'].str.startswith('flu-seqneut-2026_H3N2'),
        (~repool_df['shortname'].str.startswith('flu-seqneut-2026')) & repool_df['shortname'].str.contains('H1N1'),
        (~repool_df['shortname'].str.startswith('flu-seqneut-2026')) & repool_df['shortname'].str.contains('H3N2'),
    ]

    choices = [
        'flu-seqneut-2026_h1',
        'flu-seqneut-2026_h3',
        'old_h1_vax',
        'old_h3_vax',
    ]
    # Fail loudly on an unclassifiable shortname rather than silently assigning it
    # to a nameless subpool, where it would drop out of the per-subpool summary.
    repool_df['subpool'] = np.select(conditions, choices, default="")
    _unclassified = repool_df.loc[repool_df['subpool'] == "", 'shortname']
    assert _unclassified.empty, f"shortnames match no subpool rule: {list(_unclassified)}"

    repool_df_filtered = (
        repool_df[repool_df['strain'].isin(over_rep_strains['strain'])]
        [['strain', 'subtype', 'shortname', 'subpool', 'mean_fraction_strains']]
        .sort_values('mean_fraction_strains', ascending=False)
        .reset_index(drop=True)
    )
    repool_df_filtered
    return (repool_df_filtered,)


@app.cell
def _(pd, repool_df, repool_df_filtered):
    # Rate, not raw count: the subpools differ several-fold in size, so the number
    # of over-represented strains alone makes large subpools look worse than they
    # are. Normalize by subpool size to compare like with like.
    _n_over = repool_df_filtered['subpool'].value_counts().rename('n_over')
    _n_total = repool_df['subpool'].value_counts().rename('n_total')
    subpool_summary = (
        pd.concat([_n_over, _n_total], axis=1)
        .fillna(0)
        .astype(int)
        .assign(pct_over=lambda x: (x['n_over'] / x['n_total'] * 100).round(1))
        .sort_values('pct_over', ascending=False)
    )
    subpool_summary
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Here it appears that most of the over-represented strains were those with the highest titer in the equal volume pool (relatively `low x_volume_to_add` values from the repooling math performed previously).

    Over-representation is concentrated in the two *old vaccine* subpools: a large minority of the strains in each of `old_h3_vax` and `old_h1_vax` are over-represented, versus only a small percentage of the strains in each of the much larger `flu-seqneut-2026_h1` and `flu-seqneut-2026_h3` subpools (see the `pct_over` column above). This is consistent with the old vaccine strains having grown to the highest titers in the initial equal-volume pool.

    Both old-vaccine subpools therefore warrant remaking, not `old_h3_vax` alone. The remaining over-represented strains are scattered thinly through the two 2026 subpools, where they are a small fraction of each; whether to also adjust those individually is a separate call.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Cross-reference the over-represented strains against the volumes used in the previous repool, to confirm they line up with the low `x_volume_to_add` (i.e. highest-titer) end.
    """)
    return


@app.cell
def _(pd, repool_df_filtered, repooling_math):
    initial_repool_df=pd.read_csv(repooling_math)
    initial_repool_df
    repool_df_filtered.merge(initial_repool_df[['strain', 'x_volume_to_add']], on='strain', how='left')
    return


if __name__ == "__main__":
    app.run()
