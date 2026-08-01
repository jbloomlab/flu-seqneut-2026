import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell
def _(mo):
    mo.md(r"""
    # Analysis notebook for single virus per well infections
    Author: Caroline Kikawa.

    See README and inline annotation for details.

    Barcodes are classified by `barcode_status`: whether the barcode's expected
    strain lives in the *same* well, an *adjacent* well, a *distant* well, or is
    *not in plate* at all. `barcode_status` is computed independently per plate
    (adjacency only makes sense within a single plate's layout). Plots are
    saved to `results/` alongside this notebook.

    This notebook supports multiple plates — add or edit entries in
    `plate_configs` below as new plates come in.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import os
    from pathlib import Path

    import pandas as pd
    import altair as alt

    alt.data_transformers.disable_max_rows()
    return Path, alt, mo, os, pd


@app.cell
def _(Path, mo, os):
    # Marimo path to notebook
    notebook_directory: Path = mo.notebook_dir()

    # Neut standard set is shared across all plates below (loes2023 per config.yml)
    neut_standard_file = 'data/neut_standard_sets/loes2023_neut_standards.csv'

    # One entry per plate. viral_library_file is per-plate since it can differ
    # across runs (per config.yml).
    plate_configs = [
        {
            "plate": "20260730_single_well_H3_strains",
            "countsdir": './results/miscellaneous_plates/20260730_single_well_H3_strains/',
            "samplesfile": 'data/miscellaneous_plates/20260730_single_well_H3_strains.csv',
            "viral_library_file": 'data/viral_libraries/flu-seqneut-2026-barcode-to-strain-designed.csv',
        },
        {
            "plate": "20260730_single_well_H1_vax_old_strains",
            "countsdir": './results/miscellaneous_plates/20260730_single_well_H1_vax_old_strains/',
            "samplesfile": 'data/miscellaneous_plates/20260730_single_well_H1_vax_old_strains.csv',
            "viral_library_file": 'data/viral_libraries/flu-seqneut-2026-barcode-to-strain-designed.csv',
        },
    ]

    # Output path for saved plots (created next to the notebook)
    output_dir = notebook_directory / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)

    def find_count_and_fate_csvs(countsdir):
        """Return (count_csvs, fate_csvs) full paths for one plate's countsdir."""
        count_csvs, fate_csvs = [], []
        for f in os.listdir(countsdir):
            location = countsdir + f
            if "_counts" in f:
                count_csvs.append(location)
            elif "_fates" in f:
                fate_csvs.append(location)
        return count_csvs, fate_csvs
    return (
        find_count_and_fate_csvs,
        neut_standard_file,
        output_dir,
        plate_configs,
    )


@app.cell
def _(pd):
    def load_samples(samplesfile):
        """Load one plate's samples file into a samples_df with a 'sample' column."""
        samples_df = pd.read_csv(samplesfile)
        samples_df.drop(columns=['fastq'], inplace=True)
        samples_df['sample'] = samples_df.apply(
            lambda x: '-'.join(x.astype(str)), axis=1
        )
        return samples_df
    return (load_samples,)


@app.cell
def _(
    alt,
    find_count_and_fate_csvs,
    load_samples,
    output_dir,
    pd,
    plate_configs,
):
    # Aggregate and plot barcode fates across wells, for each plate
    all_fates = []
    for _cfg in plate_configs:
        _samples_df = load_samples(_cfg["samplesfile"])
        _samples = _samples_df["sample"].unique().tolist()
        print(f"{_cfg['plate']}: {len(_samples)} barcode runs.")

        _, _fate_csvs = find_count_and_fate_csvs(_cfg["countsdir"])

        _plate_fates = (
            pd.concat([pd.read_csv(f).assign(well=f.replace(_cfg["countsdir"], '').replace('_fates.csv', '')) for f in _fate_csvs])
            .merge(_samples_df, validate="many_to_one", on="well")
            .assign(
                fate_counts=lambda x: x.groupby("fate")["count"].transform("sum"),
            )
            .query("fate_counts > 0")[  # only keep fates with at least one count
                ["fate", "count", "well", "dilution_factor"]
            ]
            .assign(plate=_cfg["plate"])
        )
        all_fates.append(_plate_fates)

    fates = pd.concat(all_fates, ignore_index=True)
    assert len(fates) == len(fates.drop_duplicates())

    fates_chart = (
        alt.Chart(fates)
        .encode(
            alt.X("count", scale=alt.Scale(nice=False, padding=3)),
            alt.Y("well", title=None),
            alt.Color("fate", sort=sorted(fates["fate"].unique(), reverse=True)),
            alt.Order("fate", sort="descending"),
            alt.Row("plate"),
            tooltip=fates.columns.tolist(),
        )
        .mark_bar(height={"band": 0.85})
        .properties(
            height=alt.Step(10),
            width=200,
            title="Barcode parsing for initial titering plate",
        )
        .resolve_scale(y="independent")
        .configure_axis(grid=False)
    )

    fates_chart.save(output_dir / "fates_chart.html")

    fates_chart
    return


@app.cell
def _(
    find_count_and_fate_csvs,
    load_samples,
    neut_standard_file,
    pd,
    plate_configs,
):
    # Aggregate barcode counts across wells for each plate, classifying each
    # barcode as expected in the well, from a virus in an adjacent well, a
    # distant well, or not in the plate. barcode_status is computed
    # independently within each plate.

    def well_distance(
        w1: int,
        w2: int,
        n_cols: int = 12,
        n_rows: int = 8,
        row_first: bool = True,
    ) -> str:
        """
        Classify the relationship between two 1-indexed well numbers in a plate layout.

        If row_first=True (default): number left-to-right across columns, then top-to-bottom across rows.
          i.e., row-major: (row, col) -> idx = row*n_cols + col + 1

        If row_first=False: number top-to-bottom down rows, then left-to-right across columns.
          i.e., column-major: (row, col) -> idx = col*n_rows + row + 1

        Returns: "same", "adjacent" (including diagonals), or "distant".
        """
        total = n_cols * n_rows
        if not (1 <= w1 <= total) or not (1 <= w2 <= total):
            raise ValueError(f"Well numbers must be in [1, {total}]: {w1=}, {w2=}")

        if w1 == w2:
            return "same"

        def to_rc(w: int) -> tuple[int, int]:
            w0 = w - 1  # 0-index
            if row_first:
                r, c = divmod(w0, n_cols)   # row-major
            else:
                c, r = divmod(w0, n_rows)   # column-major
            return r, c

        r1, c1 = to_rc(w1)
        r2, c2 = to_rc(w2)

        return "adjacent" if (abs(r1 - r2) <= 1 and abs(c1 - c2) <= 1) else "distant"

    def classify_plate_counts(cfg, neut_standard_file):
        """Build the barcode_status-classified counts dataframe for one plate."""
        samples_df = load_samples(cfg["samplesfile"])
        samples = samples_df["sample"].unique().tolist()

        count_csvs, _ = find_count_and_fate_csvs(cfg["countsdir"])

        plate_counts = (
            pd.concat([pd.read_csv(c).assign(well=c.replace(cfg["countsdir"], '').replace('_counts.csv', '')) for c, s in zip(count_csvs, samples)])
            .merge(samples_df, validate="many_to_one", on="well")
            .drop(columns=["replicate"])
        )

        # classify barcodes as viral or neut standard
        viral_library = pd.read_csv(cfg["viral_library_file"])
        barcode_class = pd.concat(
            [
                viral_library[["barcode", "strain"]],
                pd.read_csv(neut_standard_file)[["barcode"]].assign(strain="neut_standard"),
            ],
            ignore_index=True,
        )

        assert set(plate_counts["barcode"]) == set(barcode_class["barcode"])
        plate_counts = (
            plate_counts
            .merge(barcode_class, on="barcode", validate="m:1")
            [["barcode", "count", "well", "strain"]]
            .rename(columns={"strain": "barcode_strain"})
        )

        # Map each well name to its designated/expected strain. Two well-naming
        # schemes appear in practice:
        #   - numbered wells like 'flu-seqneut-2026-H3N2_7', derived from the
        #     library's `shortname` column (e.g. 'flu-seqneut-2026_H3N2_7_bc1'
        #     with the trailing '_bcN' stripped)
        #   - wells named directly after the strain, e.g.
        #     'A-Wisconsin-67-2022_H1N1', derived from the library's `strain`
        #     column with '/' replaced by '-' (e.g. 'A/Wisconsin/67/2022_H1N1')
        # Both schemes use dashes and underscores inconsistently (e.g. '-H3N2'
        # vs '_H3N2'), so wells are matched on a normalized key with all
        # underscores converted to dashes.
        def normalize_well(w):
            return w.replace('_', '-')

        # Only derive a well name from `shortname` for rows that actually
        # follow the '..._bcN' convention -- some legacy library rows use
        # other shortname suffixes (e.g. '_VS5', '_32'), and blindly slicing
        # off the last 4 characters would corrupt those and create bogus
        # well-name collisions.
        has_bc_suffix = viral_library["shortname"].str.contains(r'_bc\d+$', regex=True)

        # A few wells don't follow either standard naming convention (e.g. a
        # missing subtype suffix); map them manually to their designated strain.
        well_key_overrides = {
            "A-Netherlands-1739-2023": "A/Netherlands/1739/2023_H1N1",
        }

        well_to_strain = pd.concat(
            [
                viral_library[has_bc_suffix].assign(
                    well_key=lambda x: x["shortname"].str.replace(r'_bc\d+$', '', regex=True).map(normalize_well)
                )[["strain", "well_key"]],
                viral_library.assign(
                    well_key=lambda x: x["strain"].str.replace('/', '-', regex=False).map(normalize_well)
                )[["strain", "well_key"]],
                pd.DataFrame(
                    {"well_key": list(well_key_overrides.keys()), "strain": list(well_key_overrides.values())}
                ),
            ],
            ignore_index=True,
        ).drop_duplicates()
        assert well_to_strain["well_key"].is_unique, "well_key maps to more than one strain"

        # Assign well_number from each well's position in samples_df BEFORE
        # excluding any unmatched wells, so remaining wells keep their correct
        # physical plate position (excluding first would shift everything
        # after it).
        well_strain_number_all = (
            samples_df[["well"]].reset_index(names="well_number")
            .assign(well_number=lambda x: x["well_number"] + 1)
            .assign(well_key=lambda x: x["well"].map(normalize_well))
        )

        unmatched_well_keys = sorted(set(well_strain_number_all["well_key"]) - set(well_to_strain["well_key"]))
        if unmatched_well_keys:
            unmatched_names = well_strain_number_all.loc[
                well_strain_number_all["well_key"].isin(unmatched_well_keys), "well"
            ].tolist()
            print(f"WARNING [{cfg['plate']}]: excluding well(s) with no matching strain in the library: {unmatched_names}")

        well_strain_number = well_strain_number_all.merge(well_to_strain, on="well_key", how="inner")

        plate_counts["well_key"] = plate_counts["well"].map(normalize_well)
        plate_counts = plate_counts[~plate_counts["well_key"].isin(unmatched_well_keys)]

        plate_counts = (
            plate_counts
            .merge(well_to_strain, on="well_key", validate="m:1")
            .rename(columns={"strain": "well_strain"})
            .drop(columns=["well_key"])
            .merge(well_strain_number[["well", "well_number"]], on="well", validate="m:1")
            .merge(
                well_strain_number[["strain", "well_number"]].rename(
                    columns={"strain": "barcode_strain", "well_number": "barcode_well_number"}
                ),
                on="barcode_strain",
                validate="many_to_one",
                how="left",
            )
        )

        plate_counts["barcode_status"] = plate_counts.apply(
            lambda row: (
                "not in plate"
                if pd.isnull(row["barcode_well_number"])
                else well_distance(row["well_number"], row["barcode_well_number"], row_first=False)
            ),
            axis=1,
        )

        plate_counts["plate"] = cfg["plate"]
        return plate_counts

    counts = pd.concat(
        [classify_plate_counts(cfg, neut_standard_file) for cfg in plate_configs],
        ignore_index=True,
    )

    counts
    return (counts,)


@app.cell
def _(alt, counts, output_dir, pd):
    # Plot how many counts are of each type (same, adjacent, distant, not in
    # plate), faceted by plate

    count_by_barcode_status = (
        counts
        .query("barcode_strain != 'neut_standard'")  # this chart is about virus, not neut standard, counts
        .sort_values("count", ascending=False)
        .groupby(["plate", "well", "well_strain", "well_number", "barcode_status"], as_index=False)
        .aggregate(
            count=pd.NamedAgg("count", "sum"),
            top_strain_w_status=pd.NamedAgg("barcode_strain", "first"),
            top_strain_w_status_counts=pd.NamedAgg("count", "first"),
        )
        .assign(
            well_counts=lambda x: x.groupby(["plate", "well"])["count"].transform("sum"),
            frac=lambda x: x["count"] / x["well_counts"],
        )
        .melt(
            id_vars=["plate", "well", "well_number", "well_strain", "barcode_status", "top_strain_w_status", "top_strain_w_status_counts", "well_counts"],
            value_vars=["count", "frac"],
            var_name="stat_type",
            value_name="stat",
        )
    )

    count_by_barcode_status_chart = (
        alt.Chart(count_by_barcode_status)
        .encode(
            alt.X("stat"),
            alt.Y("well", sort=alt.SortField("well_number")),
            alt.Color("barcode_status"),
            alt.Column("stat_type", title=None),
            alt.Row("plate"),
            tooltip=["plate", "well", "well_number",
                     "well_strain", "barcode_status", "top_strain_w_status", "top_strain_w_status_counts", "well_counts"],
        )
        .mark_bar()
        .resolve_scale(x="independent", y="independent")
        .properties(
            height=alt.Step(10),
            width=200,
            title="Barcode counts from virus in same well, adjacent well, distant well, or not in plate",
        )
    )

    count_by_barcode_status_chart.save(output_dir / "count_by_barcode_status.html")

    count_by_barcode_status_chart
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
