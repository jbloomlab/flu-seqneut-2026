# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.23.14"
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
        repooling_math = context["output"]["repooling_math"]
        dropped_strains = context["output"]["dropped_strains"]
    else:
        # Interactive stub: fill in parameters.
        viral_library_csv = '../../data/viral_libraries/flu-seqneut-2026-barcode-to-strain-designed.csv'
        neut_standard_set_csv = '../../data/neut_standard_sets/loes2023_neut_standards.csv'
        samplesfile = '../../data/miscellaneous_plates/2026-07-15_equal_vol_pool.csv'
        platedir = '../../results/miscellaneous_plates/20260715_equal_vol_pool/'
        repooling_math = '../results/pooling_math/2026-07-15_repooling_math.csv'
        dropped_strains = '../results/pooling_math/dropped_strains.csv'

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
        dropped_strains,
        neut_standard_set_csv,
        platedir,
        repooling_math,
        samplesfile,
        viral_library_csv,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Process plate counts to get ratios of variants in initial pool

    An initial pooled library was made by adding equal volumes of all variants and then infecting this pool on MDCK-SIAT1 cells. These infections were done by serially diluting a starting 50 uL volume of pool. Barcodes were then isolated and sequenced so that we can determine representation of each variant in the pool, as well as the appropriate library pool dilution (i.e., MOI) to use in neutralization assays.

    The plots generated by this notebook are interactive, so you can mouseover points for details, use the mouse-scroll to zoom and pan, and use interactive dropdowns at the bottom of the plots.

    ## Setup
    Import Python modules:
    """)
    return


@app.cell
def _():
    import altair as alt

    import numpy as np
    import string
    import pandas as pd
    from os.path import join

    _ = alt.data_transformers.disable_max_rows()

    # Basic color palette
    color_palette = [
        '#345995', #blue
        '#03cea4', #teal
        '#ca1551', #red
        '#eac435', #yellow
                   ]
    return alt, color_palette, np, pd, string


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

    counts
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
    # drop wells failing QC
    avg_barcode_counts_per_well_drops = list(avg_barcode_counts.query("fails_qc")["well"])
    return


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
    # drop wells failing QC
    min_neut_standard_frac_per_well_drops = list(
        neut_standard_fracs.query("fails_qc")["well"]
    )
    return


@app.cell
def _(alt, neut_standard_fracs):
    # Scatterplot of the same data as above, plotted by dilution factor
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
    ## Rebalancing strains contained in the library
    Viruses were rescued and blind passaged individually. To make the initial pool, we added equal volumes of all strains together and re-infected MDCK-SIAT1 cells. Now we can assess the contribution of each strain to the pool, and determine how much should be added of each virus to achieve more equal balancing.

    Each of the 2-4 viral barcodes associated with each strain were pooled prior to rescue, so they cannot be balanced.
    """)
    return


@app.cell
def _(counts):
    # Get summed barcode counts for all strains across all wells
    straincounts_allbarcodes = (
        counts.groupby(['sample', 'sample_well', 'strain', 'dilution_factor', 'serum', 'well'])['count']
        .sum()
        .reset_index()
    )

    # Get sum of all virus/barcode counts per well
    sumperwell = (
        straincounts_allbarcodes.groupby(['sample', 'dilution_factor', 'serum', 'well'])['count']
        .sum()
        .reset_index()
        .rename(columns={'count': 'counts_perwell'})
    )

    # Merge dataframes and calculate fraction of each well devoted to each strain
    merged_df = straincounts_allbarcodes.merge(sumperwell, on=['sample','dilution_factor','serum','well'])
    # Count how many barcodes were used for each strain (this can vary, e.g. 1-4)
    barcodes_per_strain = (
        counts[['strain', 'barcode']]
        .drop_duplicates()
        .groupby('strain')
        .size()
        .rename('n_barcodes')
        .reset_index()
    )
    print(barcodes_per_strain['n_barcodes'].value_counts())
    merged_df = merged_df.merge(barcodes_per_strain, on='strain', validate='many_to_one')
    merged_df['fraction_strain'] = (
        merged_df['count'] / merged_df['counts_perwell'] / merged_df['n_barcodes']
    )
    #### OLD CODE ASSUMED 2 BARCODES PER WELL
    # merged_df['fraction_strain'] = merged_df['count'] /merged_df['counts_perwell'] / 2
    ####
    merged_df
    return (merged_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We now have this fraction of reads devoted to all strains calculated for all wells. However, ideally we should just focus on those wells containing dilutions that we would use for actual neutralization assays. We should choose a set of replicate wells where the fraction of neutralization standard reads begins to increase linearly with the increasing reciprocal dilution factor. See plots above for choosing these wells.
    """)
    return


@app.cell
def _(merged_df):
    # Choose wells
    well1 = 'A8'
    well2 = 'B8'

    # Choose a pair of replicate wells near the beginning of the linear range
    single_well = merged_df.loc[merged_df['sample'].str.contains(f'{well1}-|{well2}-')]
    return single_well, well1, well2


@app.cell
def _(single_well):
    # Calculate mean fraction strain across both wells
    mean_df = single_well.groupby(['strain'])['fraction_strain'].mean().to_frame().rename(columns = {'fraction_strain': 'mean_fraction_strains'}).reset_index()
    mean_single_well = single_well.merge(mean_df, on = 'strain', how = 'left')

    # calcualte ratios to add for equal pool
    num_strains = len(mean_single_well.strain.unique())
    mean_single_well['ratio_to_add'] = (1/num_strains)/mean_single_well['fraction_strain']
    mean_single_well['mean_ratio_to_add'] = (1/num_strains)/mean_single_well['mean_fraction_strains']

    mean_single_well['est_tcid50'] = (mean_single_well['mean_fraction_strains']*25000)*76

    print(f'this library has {num_strains} total strains')
    print('stats where there isnt 0 of a virus...')
    print(mean_single_well.query('mean_ratio_to_add != inf')[['mean_ratio_to_add']].describe())

    print('\nviruses with 0 titer...')
    print(mean_single_well.query('mean_ratio_to_add == inf').strain.unique())

    RATIO_CUTOFF = 250
    print(f'\nviruses with >0 titer but ratio >={RATIO_CUTOFF} to increase...')
    print(mean_single_well.query('mean_ratio_to_add != inf').query(f'mean_ratio_to_add >= {RATIO_CUTOFF}').strain.unique())
    return mean_single_well, num_strains


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Re-pooling calculations
    """)
    return


@app.cell
def _(pd, viral_library_csv):
    # Factor to multiply ratios by
    repool_factor = 15

    # Get library IDs
    lib_id_df=pd.read_csv(viral_library_csv)
    return lib_id_df, repool_factor


@app.cell
def _(lib_id_df, mean_single_well, pd, repool_factor, well2):
    # Make repool dataframe using chosen well 
    repool_df = (mean_single_well
                 .query('mean_ratio_to_add != inf')
                 .query(f'well == "{well2}"')
                 [['strain','fraction_strain','mean_ratio_to_add']]
                 .assign(x_volume_to_add = lambda x: x['mean_ratio_to_add'] * repool_factor)
                 .merge(lib_id_df, how='outer')
                 .assign(
                     subtype = lambda x: x['strain'].str.rsplit('_', n=1).str[-1],
                     number = lambda x: pd.to_numeric(
                         x['shortname'].str.replace('flu-seqneut-2026_','').str.split('_').str[1],  # removed inner lambda
                         errors='coerce').fillna(1e6).astype(int),  # use a big number to push NaNs to bottom
                     strain_id = lambda x: x['subtype'] + '_' + x['number'].astype(str)
                 )
                 .sort_values(by=['subtype', 'number'], ascending=True)
                 .drop(columns=['number', 'barcode', 'nt_sequence_HA_ectodomain','protein_sequence_HA_ectodomain'])  # number is now temporary
                 .drop_duplicates()
                 .dropna(subset=['fraction_strain'])
                 .reset_index(drop=True)
    )

    # drop strains
    strains_to_drop = [
        'A/Bangkok/P2323/2025_H3N2', # low titer
        'A/California/LACPHL-INF02113/2025_H3N2', # low titer
        'A/England/1845724/2025_H3N2', # low titer
        'A/Galicia/GA-CHUAC-449/2025_H1N1', # low titer
    ]
    trimmed_repool_df = repool_df[~repool_df['strain'].isin(strains_to_drop)]
    trimmed_repool_df
    return repool_df, strains_to_drop, trimmed_repool_df


@app.cell
def _(
    dropped_strains,
    os,
    repool_df,
    repooling_math,
    strains_to_drop,
    trimmed_repool_df,
):
    # Save pooling math
    pooling_df = trimmed_repool_df.sort_values(by='x_volume_to_add', ascending=False).reset_index(drop=True)[['strain_id','x_volume_to_add','strain']]
    pooling_df['number'] = pooling_df['strain_id'].str.split('_').str[1].fillna(1e6).astype(int)
    pooling_df['subtype'] = pooling_df['strain_id'].str.split('_').str[0]
    pooling_df = (pooling_df
        .sort_values(by=['subtype','number'], ascending=True)
        .drop(columns=['number', 'subtype'])
                 )

    repooling_dir = os.path.dirname(repooling_math)
    os.makedirs(repooling_dir, exist_ok=True)
    pooling_df.to_csv(repooling_math, index=False)

    # Save info on dropped strains
    dropped_repool_df = repool_df[repool_df['strain'].isin(strains_to_drop)]
    dropped_repool_df.to_csv(dropped_strains,index=False)

    # Summarize top 10 lowest titer strains
    print('Displaying the top 10 lowest titer strains that made the cut...')
    trimmed_repool_df.sort_values(by='x_volume_to_add', ascending=False).head(10).reset_index(drop=True)
    return


@app.cell
def _(repool_factor, trimmed_repool_df):
    print(f'Adding {repool_factor}x of each strain ratio...')
    _sum = (
        sum(trimmed_repool_df.mean_ratio_to_add) * repool_factor
         )
    print(_sum, 'uL total pool')
    print('This means adding strains in volumes ranging from...')
    print(trimmed_repool_df.x_volume_to_add.min(), 'uL to ', trimmed_repool_df.x_volume_to_add.max(), 'uL')
    print('Assuming worse case scenario of 1:16 on 150k cells...')
    volume_per_plate = (50/16)*110
    print(volume_per_plate, 'uL per plate')
    print('About how many plates can I run?')
    print((_sum-10)/volume_per_plate)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Visualize barcode- and strain-level balancing in the current pool
    """)
    return


@app.cell
def _(alt, mean_single_well, num_strains, pd):
    # Plot the current fraction of each strain in the pool
    strains_chart = (
        alt.Chart(mean_single_well)
        .encode(
            alt.X(
                "fraction_strain",
                scale=alt.Scale(nice=False, padding=3),
            ),
            alt.Y("strain"),

            tooltip = ['strain', 'fraction_strain', 'est_tcid50'],
        )
    ).mark_bar(height={"band": 0.85}).properties(
            height=alt.Step(10),
            width=250,
            title="",
        ).properties(
            height = alt.Step(10),
            width = 200,
            title = "Strain representation, initial pool")

    # add veritcal line where we would expect equal representation of all barcodes in pool
    expected_line = alt.Chart(
        pd.DataFrame({'x': [1/num_strains]})
    ).mark_rule(strokeDash = [2,2], strokeWidth = 2).encode(x = 'x')

    # plot both barcode counts and expected line
    strains_chart + expected_line
    return


@app.cell
def _(alt, color_palette, counts, string, well1, well2):
    # Each barcode fraction across strains
    all_barcode_counts = counts[['strain', 'barcode', 'count', 'well']].dropna()
    single_well_all_barcode_counts = all_barcode_counts[all_barcode_counts['well'].isin([f'{well1}',f'{well2}'])]

    # Get tidy single well means
    tidy_single_well = single_well_all_barcode_counts[['strain','barcode','count']].groupby(['strain', 'barcode']).mean().reset_index()
    # Get sums for each strain
    strain_sums_df = tidy_single_well.groupby('strain').sum().rename(columns = {'count': 'strain_count_sum'}).reset_index()
    # Merge and calculate per strain the fraction represented by each barcode
    tidy_single_well = tidy_single_well.merge(strain_sums_df[['strain', 'strain_count_sum']], 
                           on = ['strain'],
                           validate="many_to_one",
                          )
    tidy_single_well['per_strain_fraction_barcode'] = tidy_single_well['count'] / tidy_single_well['strain_count_sum']
    tidy_single_well['barcode_letter'] = tidy_single_well.groupby('strain').cumcount().apply(lambda x: string.ascii_uppercase[x])

    # Plot as colored bar chart
    bar_chart = alt.Chart(tidy_single_well).mark_bar(height={"band": 0.85}).encode(
        x = 'per_strain_fraction_barcode',
        y = 'strain',
        color=alt.Color('barcode_letter', legend=None).scale(range=color_palette),
        tooltip = ['strain', 'per_strain_fraction_barcode', 'barcode'],
    ).configure_axis(grid=False).properties(
            height = alt.Step(10),
            width = 200,
            title = "Barcode fraction for each strain, initial pool")

    bar_chart
    return


if __name__ == "__main__":
    app.run()
