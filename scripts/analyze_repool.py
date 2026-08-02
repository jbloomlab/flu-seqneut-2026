"""Analyze how well a balanced re-pool achieved equal strain representation.

Writes a self-contained HTML report, a CSV of every strain's measured representation and
the volume to add for a corrective re-pool, and a CSV of the strains dropped from that
re-pooling. Everything specific to a re-pool comes from the `analyze_repools` section of
`config.yml` and is passed in via `snakemake.params`, so this script is not specific to
any one re-pool.

This is the companion to `analyze_pool.py`, which measures the *initial, equal-volume*
pool. Two differences between the two are deliberate and easy to "fix" into bugs:

  - **`fraction_strain` is not divided by the number of barcodes here.** `count` is already
    summed over a strain's barcodes, so the per-well fractions sum to 1 across strains, and
    each strain is compared against an absolute `1 / n_strains` expectation. Dividing by
    `n_barcodes` rescales strains unequally, breaks that sum (it comes to about 0.5 for
    this library), and makes the comparison meaningless. `analyze_pool.py` does divide,
    because it balances *per barcode*: there, each barcode is a replicate in the
    neutralization assays. Same input, different question. The sum-to-1 assertion below is
    what guards this.

  - **The corrective volume is `V_previous / fraction_strain`, not `1 / fraction_strain`.**
    The pool being measured was itself mixed at unequal per-strain volumes, so a strain's
    share of the reads is proportional to `V_previous x titer` rather than to titer alone.
    Recovering the titer therefore needs the previous volumes, and leaving them out gives
    volumes that reproduce the current imbalance almost exactly rather than correcting it.
    This is only valid because the corrective pool is mixed **from the individual strain
    stocks**, as the previous re-pool was. If it were instead made by topping up the
    existing pool with aliquots of itself, the right rule would be `1 / fraction_strain`.
    Nothing in the data reveals which was done, so the premise is stated here and in the
    report rather than inferred.

For the same reason, `analyze_pool.py`'s relative-titer and re-pool-dilution calculation is
deliberately *not* carried over: it assumes the pool it measured was mixed at equal volumes,
which is false here, and it would report a confidently wrong dilution. The corrective pool
has to be titrated.

Ordering constraint: the report is built by appending to `report` in order, and one slot is
reserved by index for the warnings, so nothing may be inserted ahead of that slot.

"""

import math
import re
import sys
import textwrap
from pathlib import Path

import altair as alt
import markdown
import pandas as pd

sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

alt.data_transformers.disable_max_rows()

# pixels allowed for y-axis labels, enough for the longest strain name
Y_LABEL_LIMIT = 300

# Diluting the pool should divide the viral counts by the dilution factor while leaving the
# neutralization standard alone, so the ratio of viral to neutralization-standard counts
# falls along a line of slope IDEAL_SLOPE against the dilution factor on log-log axes for as
# long as that holds. These are the same values `analyze_pool.py` uses, and the same
# judgement calls; see that script for why the runs are short and local.
#
# This block is duplicated from `analyze_pool.py` on purpose rather than shared. The two
# scripts have not diverged yet, so it is not yet clear which parts are genuinely common,
# and a shared statistic that quietly drifts is worse here than a little duplication. If a
# change ever has to be made to both copies, that is the signal to extract them.
IDEAL_SLOPE = -1.0
LINEARITY_TOL = 0.3
LINEARITY_Z = 3.0
ONSET_WINDOW_DILUTIONS = 3

# factor by which the odds axis is padded beyond the wells plotted on it
ODDS_AXIS_PAD = 1.5

# Colors for the barcodes of a strain, carried over from the notebook this analysis
# replaces. There are four because no strain in the library has more than four barcodes;
# Altair recycles them if one ever does, which would only make two of a strain's barcodes
# share a color rather than mislabel anything.
BARCODE_PALETTE = ["#345995", "#03cea4", "#ca1551", "#eac435"]

# How a strain's representation is labelled and colored: too scarce, about right, too
# abundant. Ordered so the legend reads from under to over rather than alphabetically.
REPRESENTATION_ORDER = ["under", "ok", "over"]
REPRESENTATION_COLORS = ["#345995", "#bbbbbb", "#ca1551"]

# reason recorded for strains that cannot be re-pooled as they had no counts
NO_COUNTS_REASON = "no counts in the analyzed wells"

# columns required of the previous pool's re-pooling CSV. Nothing in `snakemake` declares
# this cross-rule contract, so it is checked explicitly: a rename in `analyze_pool.py` should
# name itself here rather than surfacing as a `KeyError` inside a merge.
PREVIOUS_POOL_COLS = ["strain", "volume_to_add_uL"]

# See `analyze_pool.py`; this is its page template, duplicated for the same reason as the
# constants above.
PAGE_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: sans-serif; max-width: 55em; margin: 2em auto; padding: 0 1em; }}
  table {{ border-collapse: collapse; margin: 1em 0; font-size: 90%; }}
  th, td {{ border: 1px solid #ccc; padding: 0.2em 0.6em; text-align: left; }}
  h1, h2 {{ border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }}
</style>
</head>
<body>
{body}
</body>
</html>
"""

# --- read and validate the configuration -----------------------------------------------

repool = snakemake.wildcards.repool
date = snakemake.params.date
repool_config = snakemake.params.repool_config

required_config_keys = {
    "miscellaneous_plate",
    "previous_pool",
    "description",
    "linear_range_wells",
    "min_avg_barcode_count_per_well",
    "min_neut_standard_frac_per_well",
    "over_representation_factor",
    "under_representation_factor",
    "total_pool_volume",
    "min_pipettable_volume_uL",
    "dilution_steps",
    "subpools",
    "strains_to_drop",
}
if set(repool_config) != required_config_keys:
    raise ValueError(
        f"configuration for re-pool {repool} must have exactly the keys "
        f"{sorted(required_config_keys)}, but has {sorted(repool_config)}"
    )

# wells fixed in the configuration, or None to detect them from the dilution series below
configured_wells = repool_config["linear_range_wells"]
if configured_wells == "calculate":
    configured_wells = None
elif isinstance(configured_wells, str) or not configured_wells:
    raise ValueError(
        f"'linear_range_wells' for re-pool {repool} must be \"calculate\" or a list of "
        f"wells, but is {configured_wells}"
    )
elif len(configured_wells) != len(set(configured_wells)):
    raise ValueError(
        f"'linear_range_wells' for re-pool {repool} has duplicates: {configured_wells}"
    )

# Both are factors away from equal representation, applied in opposite directions: a strain
# is over-represented above `equal_share * over`, and under-represented below
# `equal_share / under`. A factor below one would put the threshold on the wrong side of
# equal representation and silently invert the meaning of the flag.
over_rep_factor = repool_config["over_representation_factor"]
under_rep_factor = repool_config["under_representation_factor"]
for _name, _factor in [
    ("over_representation_factor", over_rep_factor),
    ("under_representation_factor", under_rep_factor),
]:
    if not _factor >= 1:
        raise ValueError(
            f"'{_name}' for re-pool {repool} must be at least 1, but is {_factor}"
        )

min_pipettable_volume = repool_config["min_pipettable_volume_uL"]
if not min_pipettable_volume > 0:
    raise ValueError(
        f"'min_pipettable_volume_uL' for re-pool {repool} must be positive, but is "
        f"{min_pipettable_volume}"
    )
dilution_steps = sorted(repool_config["dilution_steps"])
if any(step <= 1 for step in dilution_steps):
    raise ValueError(
        f"'dilution_steps' for re-pool {repool} must all be greater than one, but are "
        f"{repool_config['dilution_steps']}"
    )

# Each subpool gives the regex that assigns strains to it and whether it is being made
# again from the individual strain stocks. Order matters and is preserved from the
# configuration, as a strain takes the first subpool whose regex it matches.
subpool_regexes = {}
subpools_to_remake = []
for _subpool, _spec in repool_config["subpools"].items():
    if set(_spec) != {"regex", "remake"}:
        raise ValueError(
            f"subpool '{_subpool}' of re-pool {repool} must have exactly the keys "
            f"['regex', 'remake'], but has {sorted(_spec)}"
        )
    if not isinstance(_spec["remake"], bool):
        raise ValueError(
            f"'remake' for subpool '{_subpool}' of re-pool {repool} must be true or "
            f"false, but is {_spec['remake']!r}"
        )
    subpool_regexes[_subpool] = _spec["regex"]
    if _spec["remake"]:
        subpools_to_remake.append(_subpool)

print(
    f"Analyzing re-pool {repool} from {date}, "
    f"linear_range_wells={repool_config['linear_range_wells']}"
)
print(f"Subpools to remake from strain stocks: {subpools_to_remake or 'none'}")

# --- build the report ------------------------------------------------------------------

report = []

# Problems that make the results questionable but not impossible to compute. They are
# collected here, logged as they are found, and rendered at the top of the report, as a
# report that fails to build cannot be used to diagnose the data that broke it.
warnings = []


def add_warning(text):
    """Record a warning to log and to show at the top of the report."""
    print(f"WARNING: {text}")
    warnings.append(text)


def add_markdown(text):
    """Render markdown ``text`` and add it to the report.

    ``textwrap.dedent`` strips the whitespace *common* to every line, so keep all the lines
    of a literal passed here at one indentation, and give a multi-line value its own call
    rather than interpolating it into an indented literal.

    """
    report.append(
        markdown.markdown(textwrap.dedent(text).strip(), extensions=["tables"])
    )


def add_chart(chart):
    """Add an Altair ``chart`` to the report."""
    report.append(chart.to_html(fullhtml=False, output_div=f"chart{len(report)}"))


def add_table(df):
    """Add the data frame ``df`` to the report as a table."""
    report.append(df.to_html(index=False, float_format=lambda v: f"{v:.4g}"))


add_markdown(f"""
    # Balance of library re-pool `{repool}`

    Analysis of the barcode sequencing of the miscellaneous plate
    `{repool_config["miscellaneous_plate"]}` from {date}.
    The plots are interactive: mouse over points for details.

    ## Experimental description
    """)
add_markdown(repool_config["description"])  # added separately as it is multi-line

# Reserved by index and filled in at the very end, once every warning has been collected.
# Nothing may be inserted into `report` ahead of it, only appended after it.
warnings_index = len(report)
report.append("")

# --- read the input data ---------------------------------------------------------------

samples = pd.read_csv(snakemake.input.samples_csv).drop(columns="fastq")
if "dilution_factor" not in samples.columns:
    raise ValueError(f"{snakemake.input.samples_csv} lacks a 'dilution_factor' column")
samples["sample"] = samples.astype(str).agg("-".join, axis=1)
assert samples["sample"].nunique() == len(samples), "samples are not uniquely named"

if configured_wells:
    missing_wells = set(configured_wells) - set(samples["well"])
    if missing_wells:
        raise ValueError(
            f"'linear_range_wells' for re-pool {repool} not on the plate: "
            f"{sorted(missing_wells)}"
        )

# order wells by the amount of pool they got, for the per-well plots
sample_order = samples.sort_values(["dilution_factor", "well"])["sample"].tolist()


def well_of(csv_file, suffix):
    """Get the well that a per-well ``suffix`` CSV of the plate is for."""
    return Path(csv_file).stem.removesuffix(suffix)


def read_per_well_csvs(csv_files, suffix):
    """Read per-well CSVs into one data frame merged with the sample information."""
    return pd.concat(
        [pd.read_csv(f).assign(well=well_of(f, suffix)) for f in csv_files],
        ignore_index=True,
    ).merge(samples, on="well", validate="many_to_one", how="left")


# --- fates of the sequencing reads -----------------------------------------------------

add_markdown("""
    ## Fates of the sequencing reads

    The "fates" of the reads parsed for each well. If most reads are not "valid barcode",
    that could indicate a problem with the sequencing or with the barcodes being parsed.
    """)

fates = read_per_well_csvs(snakemake.input.fates, "_fates")
fates = fates[fates.groupby("fate")["count"].transform("sum") > 0]  # drop unused fates

add_chart(
    alt.Chart(fates)
    .encode(
        alt.X("count", scale=alt.Scale(nice=False, padding=3)),
        alt.Y("sample", title=None, sort=sample_order),
        alt.Color("fate", sort=sorted(fates["fate"].unique(), reverse=True)),
        alt.Order("fate", sort="descending"),
        tooltip=fates.columns.tolist(),
    )
    .mark_bar(height={"band": 0.85})
    .properties(height=alt.Step(10), width=200)
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

# --- barcode counts and their QC -------------------------------------------------------

counts = read_per_well_csvs(snakemake.input.counts, "_counts")[
    ["barcode", "count", "well", "sample", "dilution_factor"]
]

viral_library = pd.read_csv(snakemake.input.viral_library)
barcode_class = pd.concat(
    [
        viral_library[["barcode", "strain"]].assign(neut_standard=False),
        pd.read_csv(snakemake.input.neut_standard_set)[["barcode"]].assign(
            neut_standard=True,
            strain=pd.NA,
        ),
    ],
    ignore_index=True,
)
if set(counts["barcode"]) != set(barcode_class["barcode"]):
    raise ValueError(
        "barcodes in the counts do not match those in the viral library and neutralization "
        "standard set"
    )
counts = counts.merge(barcode_class, on="barcode", validate="many_to_one")

add_markdown(f"""
    ## Barcode counts per well

    Wells need enough counts per barcode for the composition of the pool to be estimated
    accurately, and enough counts from the neutralization standard for the amount of pool
    in the well to be measured. Wells are flagged when they have an average of fewer than
    {repool_config["min_avg_barcode_count_per_well"]} counts per barcode or when fewer than
    {repool_config["min_neut_standard_frac_per_well"]} of their counts are from the
    neutralization standard.
    """)

well_qc = (
    counts.assign(
        neut_standard_count=lambda x: x["count"] * x["neut_standard"].astype(int)
    )
    .groupby(["well", "sample", "dilution_factor"], as_index=False)
    .aggregate(
        avg_count=pd.NamedAgg("count", "mean"),
        total_count=pd.NamedAgg("count", "sum"),
        neut_standard_count=pd.NamedAgg("neut_standard_count", "sum"),
    )
    .assign(
        # a well with no counts at all gets a null fraction rather than a zero division
        neut_standard_frac=lambda x: (
            x["neut_standard_count"] / x["total_count"].where(x["total_count"] > 0)
        ),
        fails_count_qc=lambda x: (
            x["avg_count"] < repool_config["min_avg_barcode_count_per_well"]
        ),
        # a null fraction cannot be compared, so treat it as failing
        fails_neut_standard_qc=lambda x: ~(
            x["neut_standard_frac"] >= repool_config["min_neut_standard_frac_per_well"]
        ),
        fails_qc=lambda x: x["fails_count_qc"] | x["fails_neut_standard_qc"],
        viral_count=lambda x: x["total_count"] - x["neut_standard_count"],
        neut_standard_odds=lambda x: (
            x["neut_standard_count"] / x["viral_count"].where(x["viral_count"] > 0)
        ),
        viral_to_neut_standard=lambda x: (
            x["viral_count"]
            / x["neut_standard_count"].where(x["neut_standard_count"] > 0)
        ),
        fittable=lambda x: ~x["fails_qc"]
        & (x["viral_count"] > 0)
        & (x["neut_standard_count"] > 0),
    )
)

empty_wells = list(well_qc.query("total_count <= 0")["well"])
if empty_wells:
    add_warning(f"these wells have no barcode counts at all: {empty_wells}")

# Hovering a bar outlines that well in both panels. The param has to be declared on the
# concatenated chart rather than a panel to be in scope for both, and needs `empty=False`
# or every bar matches when nothing is hovered.
well_hover = alt.selection_point(
    name="well_hover",
    fields=["well"],
    on="pointerover",
    clear="pointerout",
    empty=False,
)

well_qc_panels = []
for qc_col, x_col, x_title in [
    ("fails_count_qc", "avg_count", "average counts per barcode in well"),
    (
        "fails_neut_standard_qc",
        "neut_standard_frac",
        "fraction of counts from neutralization standard",
    ),
]:
    first_panel = not well_qc_panels  # only the first labels the shared y-axis
    well_qc_panels.append(
        alt.Chart(well_qc)
        .encode(
            alt.X(
                x_col,
                title=x_title,
                scale=alt.Scale(nice=False, padding=3),
                axis=alt.Axis(tickCount=4),
            ),
            alt.Y(
                "sample",
                title=None,
                sort=sample_order,
                scale=alt.Scale(domain=sample_order),
                axis=alt.Axis(labelLimit=Y_LABEL_LIMIT) if first_panel else None,
            ),
            alt.Fill(qc_col, legend=alt.Legend(titleLimit=500)),
            stroke=alt.value("black"),
            strokeWidth=alt.condition(well_hover, alt.value(2), alt.value(0)),
            tooltip=[
                (
                    alt.Tooltip(c, format=".3g")
                    if pd.api.types.is_float_dtype(well_qc[c])
                    else c
                )
                for c in well_qc.columns
            ],
        )
        .mark_bar(height={"band": 0.85})
        .properties(height=alt.Step(10), width=250)
    )

add_chart(
    alt.hconcat(*well_qc_panels, spacing=45)
    .add_params(well_hover)
    .resolve_scale(y="shared")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

# --- the linear range of the dilution series -------------------------------------------

add_markdown(f"""
    ## The linear range of the dilution series

    Diluting the pool should divide its viral counts by the dilution factor while leaving
    the neutralization standard alone, so the neutralization standard's share of the
    counts, in odds form, rises along a straight line against the dilution factor on
    log-log axes for as long as that holds. It is over the dilutions where that is true
    that the composition of the pool can be read.

    Runs of {ONSET_WINDOW_DILUTIONS} consecutive dilutions are read from concentrated to
    dilute, and the first whose fitted slope is within {LINEARITY_TOL} of ideal marks where
    the pool starts responding. Wells failing either QC threshold above are drawn below but
    left out of the fits, as is any well with no viral counts, whose odds are undefined.
    """)


def log_log_slope(x, y):
    """Least-squares slope of log10(``y``) on log10(``x``), and its standard error."""
    x, y = x.map(math.log10), y.map(math.log10)
    if len(x) < 3:
        raise ValueError("need at least three points to fit a slope and its error")
    x_centered, y_centered = x - x.mean(), y - y.mean()
    sxx = (x_centered**2).sum()
    if sxx == 0:
        raise ValueError("cannot fit a slope to a single distinct dilution")
    slope = (x_centered * y_centered).sum() / sxx
    sse = ((y_centered - slope * x_centered) ** 2).sum()
    return slope, math.sqrt(sse / (len(x) - 2) / sxx)


def fit_dilution_windows(fittable_wells):
    """Fit the slope of each run of ``ONSET_WINDOW_DILUTIONS`` consecutive dilutions."""
    dilutions = sorted(fittable_wells["dilution_factor"].unique())
    if len(dilutions) < ONSET_WINDOW_DILUTIONS:
        add_warning(
            f"only {len(dilutions)} dilutions are usable, fewer than the "
            f"{ONSET_WINDOW_DILUTIONS} that a run of dilutions needs, so the linear range "
            "cannot be located"
        )
        return []

    windows = []
    for i in range(len(dilutions) - ONSET_WINDOW_DILUTIONS + 1):
        window = dilutions[i : i + ONSET_WINDOW_DILUTIONS]
        in_window = fittable_wells[fittable_wells["dilution_factor"].isin(window)]
        slope, slope_se = log_log_slope(
            in_window["dilution_factor"], in_window["viral_to_neut_standard"]
        )
        deviation = abs(slope - IDEAL_SLOPE)
        windows.append(
            {
                "window": window,
                "slope": slope,
                "slope_se": slope_se,
                "deviation": deviation,
                "within_tol": deviation <= LINEARITY_TOL,
                "within_tol_widened": deviation
                <= LINEARITY_TOL + LINEARITY_Z * slope_se,
            }
        )
    print(f"slope of each run of dilutions, against an ideal of {IDEAL_SLOPE}:")
    for w in windows:
        print(
            f"    {w['window'][0]:>6}-{w['window'][-1]:<6} slope {w['slope']:6.3f} "
            f"+/- {w['slope_se']:5.3f}, off by {w['deviation']:5.3f}"
            f"{'  <-- within tolerance' if w['within_tol'] else ''}"
        )
    return windows


def linear_range_onset(windows):
    """Take the first run of dilutions whose slope is within tolerance of ideal."""
    for within in ["within_tol", "within_tol_widened"]:
        qualifying = [w for w in windows if w[within]]
        if qualifying:
            if within == "within_tol_widened":
                add_warning(
                    "no run of dilutions has a slope within "
                    f"{LINEARITY_TOL} of the ideal {IDEAL_SLOPE}, so the linear range is "
                    f"placed using a run admitted only once the tolerance is widened by "
                    f"{LINEARITY_Z} standard errors of its slope"
                )
            return qualifying[0]
    add_warning(
        f"no run of {ONSET_WINDOW_DILUTIONS} dilutions has a slope near the ideal "
        f"{IDEAL_SLOPE}, so the pool never clearly responds to being diluted"
    )
    return None


fittable_wells = well_qc.query("fittable")
if not len(fittable_wells):
    raise ValueError(
        f"no well of re-pool {repool} passes QC, so there is nothing to analyze"
    )
windows = fit_dilution_windows(fittable_wells)
onset = linear_range_onset(windows) if windows else None

if windows:
    add_markdown("The slope fit to each run of dilutions:")
    add_table(
        pd.DataFrame(
            [
                {
                    "dilutions": f"{w['window'][0]}-{w['window'][-1]}",
                    "slope": w["slope"],
                    "standard error": w["slope_se"],
                    "off ideal by": w["deviation"],
                    "within tolerance": w["within_tol"],
                    "chosen": w is onset,
                }
                for w in windows
            ]
        )
    )

# The composition is read from the middle dilution of the run that located the linear
# range, so that there is a dilution of margin on either side of the wells used.
if onset:
    calculated_dilution = onset["window"][len(onset["window"]) // 2]
    calculated_wells = sorted(
        well_qc[well_qc["dilution_factor"] == calculated_dilution]["well"]
    )
else:
    calculated_dilution, calculated_wells = None, None

if configured_wells:
    analysis_wells = list(configured_wells)
elif calculated_wells:
    analysis_wells = calculated_wells
else:
    raise ValueError(
        f"cannot find the linear range for re-pool {repool}, so set 'linear_range_wells' "
        'to the wells to read the composition from rather than to "calculate"'
    )
print(f"Reading the composition from wells {analysis_wells}")

odds_domain = [
    float(fittable_wells["neut_standard_odds"].min() / ODDS_AXIS_PAD),
    float(fittable_wells["neut_standard_odds"].max() * ODDS_AXIS_PAD),
]
dilution_axis = alt.X(
    "dilution_factor",
    scale=alt.Scale(type="log"),
    title="pool reciprocal dilution factor",
)
odds_axis = alt.Y(
    "neut_standard_odds",
    scale=alt.Scale(type="log", domain=odds_domain, nice=False),
    title="neutralization standard odds, frac / (1 - frac)",
)

odds_layers = [
    alt.Chart(well_qc)
    .mark_circle(size=60, clip=True)
    .encode(
        dilution_axis,
        odds_axis,
        alt.Color("fails_qc"),
        tooltip=[
            "well",
            "dilution_factor",
            alt.Tooltip("neut_standard_frac", format=".3g"),
            alt.Tooltip("neut_standard_odds", format=".3g"),
            "viral_count",
            "neut_standard_count",
            "fails_count_qc",
            "fails_neut_standard_qc",
        ],
    )
]

if onset:
    in_window = fittable_wells[fittable_wells["dilution_factor"].isin(onset["window"])]
    ideal_intercept = (
        in_window["neut_standard_odds"].map(math.log10)
        + IDEAL_SLOPE * in_window["dilution_factor"].map(math.log10)
    ).mean()

    def ideal_line(dilutions):
        """The ideal line over ``dilutions``, clipped to the limits set by the wells."""
        line = pd.DataFrame({"dilution_factor": sorted(dilutions)}).assign(
            neut_standard_odds=lambda x: x["dilution_factor"].map(
                lambda d: 10 ** (ideal_intercept - IDEAL_SLOPE * math.log10(d))
            )
        )
        return line[line["neut_standard_odds"].between(*odds_domain)]

    odds_layers += [
        alt.Chart(ideal_line(fittable_wells["dilution_factor"].unique()))
        .mark_line(color="black", opacity=0.25, strokeWidth=2)
        .encode(dilution_axis, odds_axis),
        alt.Chart(ideal_line(onset["window"]))
        .mark_line(color="black", strokeWidth=2)
        .encode(dilution_axis, odds_axis),
        alt.Chart(pd.DataFrame({"dilution_factor": [calculated_dilution]}))
        .mark_rule(color="black", strokeDash=[2, 2])
        .encode(dilution_axis),
    ]

add_chart(alt.layer(*odds_layers).configure_axisX(labelOverlap="greedy"))

if onset:
    add_markdown(f"""
        The pool starts responding to dilution at **{onset["window"][0]}-fold**: over
        {onset["window"][0]} to {onset["window"][-1]} the slope is
        {onset["slope"]:.3f} +/- {onset["slope_se"]:.3f}, off the ideal {IDEAL_SLOPE} by
        {onset["deviation"]:.3f}. The composition would be calculated from the **middle**
        dilution of that run, **{calculated_dilution}-fold** (wells
        {", ".join(f"`{w}`" for w in calculated_wells)}, marked by the dashed rule).
        """)

if configured_wells:
    add_markdown(f"""
        That calculation is **overridden** by `linear_range_wells` in `config.yml`, which
        sets the wells to {", ".join(f"`{w}`" for w in analysis_wells)}. The composition
        below is read from those wells.
        """)
    if calculated_wells and set(analysis_wells) != set(calculated_wells):
        add_warning(
            f"the configured `linear_range_wells` {analysis_wells} are not the wells "
            f"{calculated_wells} at the {calculated_dilution}-fold dilution that would "
            "have been calculated from the linear range"
        )
else:
    add_markdown("The composition below is read from those wells.")

failing_analysis_wells = list(
    well_qc[
        well_qc["well"].isin(analysis_wells)
        & (well_qc["fails_count_qc"] | well_qc["fails_neut_standard_qc"])
    ]["well"]
)
if failing_analysis_wells:
    add_warning(
        f"the wells used for the re-pooling math fail QC: {failing_analysis_wells}. "
        "The composition measured below, and everything computed from it, may not be "
        "reliable; choose different wells in `config.yml`."
    )

# a well with no counts gives no information on the composition, so it cannot contribute
analysis_wells_with_counts = [w for w in analysis_wells if w not in empty_wells]
if not analysis_wells_with_counts:
    raise ValueError(
        f"none of the wells used for the re-pooling math for re-pool {repool} have any "
        f"barcode counts: {analysis_wells}"
    )
if analysis_wells_with_counts != analysis_wells:
    add_warning(
        "the composition below is measured from only "
        f"{analysis_wells_with_counts}, as the other wells chosen for it have no counts."
    )

# --- representation of each strain in the re-pool --------------------------------------

viral_counts = counts[counts["strain"].notna()]
assert (
    counts.loc[counts["neut_standard"], "strain"].isna().all()
), "neutralization standard barcodes unexpectedly carry a strain"

n_barcodes = (
    viral_counts[["strain", "barcode"]]
    .drop_duplicates()
    .groupby("strain", as_index=False)
    .size()
    .rename(columns={"size": "n_barcodes"})
)
print(f"Barcodes per strain:\n{n_barcodes['n_barcodes'].value_counts()}")

strain_fracs = (
    viral_counts.groupby(["well", "sample", "strain"], as_index=False)
    .aggregate(count=pd.NamedAgg("count", "sum"))
    .assign(well_count=lambda x: x.groupby("well")["count"].transform("sum"))
    .merge(n_barcodes, on="strain", validate="many_to_one")
    # The strain's share of the viral counts. `count` is already summed over the strain's
    # barcodes, so this is NOT divided by `n_barcodes`; see the module docstring.
    .assign(fraction_strain=lambda x: x["count"] / x["well_count"])
)

# The guard for the decision above: shares that are summed over barcodes must sum to one
# across strains in every well. Dividing by `n_barcodes` breaks this immediately.
_frac_sums = strain_fracs.groupby("well")["fraction_strain"].sum()
assert (
    _frac_sums.sub(1).abs().lt(1e-9).all()
), f"'fraction_strain' does not sum to 1 per well:\n{_frac_sums}"

# strain-level short names, dropping the per-barcode suffix
strain_names = (
    viral_library.assign(
        shortname=lambda x: x["shortname"].str.replace(r"_bc\d+$", "", regex=True)
    )[["strain", "shortname"]]
    .drop_duplicates()
    .reset_index(drop=True)
)
if strain_names["strain"].nunique() != len(strain_names):
    raise ValueError(
        "viral library metadata varies within a strain; 'shortname' must be constant "
        "across a strain's barcodes"
    )


def assign_subpool(shortname):
    """Get the subpool for ``shortname``, the first one whose regex it matches."""
    for subpool, regex in subpool_regexes.items():
        if re.search(regex, shortname):
            return subpool
    raise ValueError(f"{shortname} does not match any regex in 'subpools'")


def sort_key(names):
    """Sort key that orders the numbers in ``names`` numerically."""
    return names.str.replace(r"\d+", lambda m: m.group().zfill(10), regex=True)


strain_means = (
    strain_fracs[strain_fracs["well"].isin(analysis_wells_with_counts)]
    .groupby("strain", as_index=False)
    .aggregate(
        mean_fraction_strain=pd.NamedAgg("fraction_strain", "mean"),
        n_wells=pd.NamedAgg("well", "nunique"),
    )
    .merge(strain_names, on="strain", validate="one_to_one")
    .assign(subpool=lambda x: x["shortname"].map(assign_subpool))
    .sort_values(["subpool", "shortname"], key=sort_key)
    .reset_index(drop=True)
)
assert (
    strain_means["n_wells"] == len(analysis_wells_with_counts)
).all(), "wells missing a strain"
n_strains = len(strain_means)
equal_share = 1 / n_strains

# The charts order strains by how much of the pool they make up, most abundant at the top,
# so the imbalance reads off the shape of the bars and the over-represented strains gather
# together. The data frames keep their subpool ordering, which is what the CSVs are used
# for at the bench.
strain_order = (
    strain_means.sort_values("mean_fraction_strain", ascending=False)["strain"]
    .astype(str)
    .tolist()
)

# --- the previous pool's volumes -------------------------------------------------------

previous_repooling = pd.read_csv(snakemake.input.previous_repooling_math)
missing_cols = [c for c in PREVIOUS_POOL_COLS if c not in previous_repooling.columns]
if missing_cols:
    raise ValueError(
        f"{snakemake.input.previous_repooling_math}, the re-pooling math of pool "
        f"{repool_config['previous_pool']} produced by the `analyze_pool` rule, lacks the "
        f"columns {missing_cols}; a change to that rule's output needs this script updated"
    )

# Every strain went into the pool at some positive volume. A zero or negative one would
# otherwise reach the division and the logarithm below as a bare "math domain error" that
# names neither the strain nor the file it came from.
_bad_volumes = previous_repooling.loc[
    ~(previous_repooling["volume_to_add_uL"] > 0), "strain"
]
if len(_bad_volumes):
    raise ValueError(
        f"{snakemake.input.previous_repooling_math} gives a 'volume_to_add_uL' that is not "
        f"positive for {sorted(_bad_volumes)}; every strain in the pool was added at some "
        "positive volume"
    )

# The strains measured here and the strains that were re-pooled must be the same set. They
# are maintained independently -- `drop_strains.csv` builds the viral library, while
# `strains_to_drop` in `analyze_pools` decides what got re-pooled -- so this is asserted
# rather than assumed, and names the divergence if the two lists ever drift apart.
measured_strains = set(strain_means["strain"])
previous_strains = set(previous_repooling["strain"])
if measured_strains != previous_strains:
    raise ValueError(
        "the strains in this re-pool do not match the strains that were re-pooled into "
        f"it. Measured but not re-pooled: {sorted(measured_strains - previous_strains)}. "
        f"Re-pooled but not measured: {sorted(previous_strains - measured_strains)}"
    )

strain_means = strain_means.merge(
    previous_repooling[PREVIOUS_POOL_COLS].rename(
        columns={"volume_to_add_uL": "previous_volume_uL"}
    ),
    on="strain",
    validate="one_to_one",
)

add_markdown(f"""
    ## Representation of each strain in the re-pool

    The left panel shows each strain's share of the viral counts (the neutralization
    standard is excluded) in wells
    {", ".join(f"`{w}`" for w in analysis_wells_with_counts)}. A strain's counts are summed
    over its barcodes, so the shares sum to one across strains and the dashed line at
    {equal_share:.2%} marks where every strain would be equally represented. The black
    points are the individual wells, so replicate disagreement stays visible rather than
    being averaged away. The right panel shows how a strain's counts are split among its
    barcodes, which should be roughly even unless a barcode is poorly represented in the
    rescued virus.

    A strain is called **over-represented** when its share exceeds that equal share by more
    than the `over_representation_factor` of {over_rep_factor} set in `config.yml` (above
    {over_rep_factor * equal_share:.2%}), and **under-represented** when its share falls
    below that share divided by the `under_representation_factor` of {under_rep_factor}
    (below {equal_share / under_rep_factor:.2%}). The two thresholds are reciprocal, so
    being this many fold too abundant and this many fold too scarce count as equally far
    from equal representation.
    """)

strain_means["over_represented"] = strain_means["mean_fraction_strain"] > (
    over_rep_factor * equal_share
)
strain_means["under_represented"] = strain_means["mean_fraction_strain"] < (
    equal_share / under_rep_factor
)
# A strain cannot be both, as `equal_share / under <= equal_share <= equal_share * over`
# whenever both factors are at least one, which is enforced above.
assert not (
    strain_means["over_represented"] & strain_means["under_represented"]
).any(), "a strain is flagged as both over- and under-represented"

strain_means["representation"] = (
    strain_means["over_represented"]
    .map({True: "over", False: "ok"})
    .where(~strain_means["under_represented"], "under")
)

# hovering a bar in either panel outlines that strain in both
strain_hover = alt.selection_point(
    name="strain_hover",
    fields=["strain"],
    on="pointerover",
    clear="pointerout",
    empty=False,
)

# Individual wells overlaid on the mean, so replicate disagreement stays visible rather
# than being averaged away.
replicate_fracs = strain_fracs[strain_fracs["well"].isin(analysis_wells_with_counts)][
    ["strain", "well", "fraction_strain"]
]

# How each strain's counts split among its barcodes, for the right panel. Averaged over the
# analysis wells first so a well with more reads does not dominate the split.
barcode_fracs = (
    viral_counts[viral_counts["well"].isin(analysis_wells_with_counts)]
    .groupby(["strain", "barcode"], as_index=False)
    .aggregate(count=pd.NamedAgg("count", "mean"))
    .assign(
        strain_count=lambda x: x.groupby("strain")["count"].transform("sum"),
        fraction_barcode=lambda x: x["count"] / x["strain_count"],
        barcode_number=lambda x: x.groupby("strain").cumcount() + 1,
    )
)

strain_panel = (
    alt.Chart(strain_means)
    .mark_bar(height={"band": 0.85})
    .encode(
        alt.X(
            "mean_fraction_strain",
            title="fraction of pool from each strain",
            scale=alt.Scale(nice=False, padding=3),
        ),
        alt.Y(
            "strain",
            title=None,
            sort=strain_order,
            scale=alt.Scale(domain=strain_order),
            axis=alt.Axis(labelLimit=Y_LABEL_LIMIT),
        ),
        alt.Fill(
            "representation",
            title="representation",
            scale=alt.Scale(domain=REPRESENTATION_ORDER, range=REPRESENTATION_COLORS),
        ),
        stroke=alt.value("black"),
        strokeWidth=alt.condition(strain_hover, alt.value(2), alt.value(0)),
        tooltip=[
            "strain",
            "shortname",
            "subpool",
            alt.Tooltip("mean_fraction_strain", format=".3g"),
            "representation",
        ],
    )
    + alt.Chart(replicate_fracs)
    .mark_point(size=15, filled=True, opacity=0.9, color="black")
    .encode(
        alt.X("fraction_strain"),
        alt.Y("strain", sort=strain_order, scale=alt.Scale(domain=strain_order)),
        tooltip=["strain", "well", alt.Tooltip("fraction_strain", format=".3g")],
    )
    # the solid line marks where every strain would be equally represented, and the dashed
    # lines the two thresholds either side of it
    + alt.Chart(pd.DataFrame({"equal representation": [equal_share]}))
    .mark_rule(strokeWidth=2, color="red")
    .encode(alt.X("equal representation"))
    + alt.Chart(
        pd.DataFrame(
            {
                "threshold": [
                    equal_share / under_rep_factor,
                    equal_share * over_rep_factor,
                ]
            }
        )
    )
    .mark_rule(strokeDash=[2, 2], strokeWidth=1.5, color="red", opacity=0.6)
    .encode(alt.X("threshold"))
)

barcode_panel = (
    alt.Chart(barcode_fracs)
    .mark_bar(height={"band": 0.85})
    .encode(
        # A strain's barcodes stack into one full-width bar, so what is read here is where
        # the color boundaries fall rather than the length of any one segment. `stack` and
        # `order` are set explicitly: the default stacking is the same, but the segment
        # order would otherwise follow the data rather than the barcode number, which makes
        # the colors jump between strains and the panel hard to scan.
        alt.X(
            "fraction_barcode",
            title="fraction of strain from each barcode",
            stack="zero",
            scale=alt.Scale(domain=[0, 1], nice=False),
        ),
        alt.Order("barcode_number:N"),
        # y-axis labels come from the strain panel that shares this axis
        alt.Y(
            "strain",
            title=None,
            sort=strain_order,
            scale=alt.Scale(domain=strain_order),
            axis=None,
        ),
        alt.Fill(
            "barcode_number:N",
            legend=None,
            scale=alt.Scale(range=BARCODE_PALETTE),
        ),
        # White rather than transparent when not hovered, so the boundary between a
        # strain's barcodes stays visible at this row height even when two segments are
        # similar colors. The encoding channel overrides any stroke set on the mark, so
        # the hairline has to be the unhovered branch of this condition.
        stroke=alt.condition(strain_hover, alt.value("black"), alt.value("white")),
        strokeWidth=alt.condition(strain_hover, alt.value(2), alt.value(0.3)),
        tooltip=[
            "strain",
            "barcode",
            alt.Tooltip("fraction_barcode", format=".3g"),
        ],
    )
)

add_chart(
    alt.hconcat(
        *(
            p.properties(height=alt.Step(10), width=240)
            for p in [strain_panel, barcode_panel]
        )
    )
    .add_params(strain_hover)
    # The y scale is shared so the two panels line up strain for strain. The fill scale
    # must NOT be: the panels fill by different fields (`representation` on the left,
    # `barcode_number` on the right), and `hconcat` shares fill by default, which merges
    # the two domains onto one scale. The barcode values then land where that scale has no
    # color and the whole right panel is drawn white on white, which looks like the panel
    # failing to render rather than a color bug.
    .resolve_scale(y="shared", fill="independent")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

n_over = int(strain_means["over_represented"].sum())
n_under = int(strain_means["under_represented"].sum())
most = strain_means.nlargest(1, "mean_fraction_strain").iloc[0]
least = strain_means.nsmallest(1, "mean_fraction_strain").iloc[0]
add_markdown(f"""
    **{n_over}** of the {n_strains} strains are over-represented and **{n_under}** are
    under-represented. The most abundant, `{most["strain"]}`, is at
    {most["mean_fraction_strain"]:.2%} of the pool, which is
    **{most["mean_fraction_strain"] / equal_share:.1f}x** its equal share; the least
    abundant, `{least["strain"]}`, is at {least["mean_fraction_strain"]:.2%}, or
    **{least["mean_fraction_strain"] / equal_share:.2f}x**.
    """)


def representation_table(df):
    """Table of strains with their share of the pool relative to an equal share."""
    return df.assign(x_expected=lambda x: x["mean_fraction_strain"] / equal_share)[
        ["shortname", "strain", "subpool", "mean_fraction_strain", "x_expected"]
    ].reset_index(drop=True)


add_markdown("The 10 most over-represented strains:")
add_table(representation_table(strain_means.nlargest(10, "mean_fraction_strain")))
add_markdown("The 10 most under-represented strains:")
add_table(representation_table(strain_means.nsmallest(10, "mean_fraction_strain")))

# --- representation by subpool ---------------------------------------------------------

add_markdown("""
    ### Representation by subpool

    As *rates* rather than counts: the subpools differ several-fold in size, so the number
    of flagged strains alone would make the large subpools look worse than they are.
    """)

subpool_summary = (
    strain_means.groupby("subpool", as_index=False)
    .aggregate(
        n_strains=pd.NamedAgg("strain", "nunique"),
        n_over_represented=pd.NamedAgg("over_represented", "sum"),
        n_under_represented=pd.NamedAgg("under_represented", "sum"),
        mean_fraction_strain=pd.NamedAgg("mean_fraction_strain", "mean"),
    )
    .assign(
        percent_over_represented=lambda x: (
            100 * x["n_over_represented"] / x["n_strains"]
        ),
        percent_under_represented=lambda x: (
            100 * x["n_under_represented"] / x["n_strains"]
        ),
    )
    .sort_values("percent_over_represented", ascending=False)
    .reset_index(drop=True)
)
add_table(subpool_summary)

# --- cross-reference against the previous pool's volumes -------------------------------

add_markdown(f"""
    ### Against the volumes that made this pool

    Each strain's representation here against the volume of it that went into the pool,
    taken from the `analyze_pool` report for `{repool_config["previous_pool"]}`. A strain
    that needed little volume was growing well, so if the imbalance tracks stock titer the
    over-represented strains gather at the low-volume end.
    """)

add_chart(
    alt.Chart(strain_means)
    .mark_circle(size=60)
    .encode(
        alt.X(
            "previous_volume_uL",
            scale=alt.Scale(type="log"),
            title="volume added to make this pool (uL)",
        ),
        alt.Y(
            "mean_fraction_strain",
            scale=alt.Scale(type="log"),
            title="fraction of pool",
        ),
        alt.Color(
            "representation",
            title="representation",
            scale=alt.Scale(domain=REPRESENTATION_ORDER, range=REPRESENTATION_COLORS),
        ),
        tooltip=[
            "strain",
            "shortname",
            "subpool",
            "representation",
            alt.Tooltip("mean_fraction_strain", format=".3g"),
            alt.Tooltip("previous_volume_uL", format=".4g"),
        ],
    )
    .properties(height=300, width=400)
    .configure_axis(grid=False)
)

# --- corrective re-pooling volumes -----------------------------------------------------

add_markdown(f"""
    ## Corrective re-pooling

    The strains are held in subpools, which are combined to make the final pool, so the
    imbalance has two separate causes and two separate remedies:

    1. **Within a subpool**, strains can be out of balance with each other. Fixing that
       means making the subpool again from the individual strain stocks, which is
       pipetting work proportional to the number of strains in it. A subpool is remade
       only where `remake` is set for it in `config.yml`.
    2. **Between subpools**, a subpool can contribute too much or too little to the final
       pool. Fixing that is just a matter of how much of each subpool is added when they
       are combined, so it is done for every subpool whether or not it is being remade.

    The target is an equal share for every strain: each subpool contributes in proportion
    to the number of strains it holds, so that all {n_strains} strains aim for the same
    {equal_share:.2%} of the final pool.

    The arithmetic rests on the same premise as before, stated rather than inferred
    because nothing in the data reveals it: a strain's share of the reads is proportional
    to the volume of it that went in *times* the titer of its stock, so the stock titer is
    proportional to `fraction_strain / previous_volume`. Within a remade subpool the volume
    that equalizes its strains is therefore proportional to
    `previous_volume / fraction_strain`.

    Note this predicts the *composition* only. It says nothing about the titer of the
    corrective pool, so unlike the initial pool no dilution is recommended here: the
    corrective pool has to be titrated.
    """)

# Strains excluded from the corrective re-pooling: those named in the config, plus any with
# no counts, which would need an infinite volume. Determined before the ratios and volumes,
# which are computed over just the strains that are kept.
drop_reasons = dict(repool_config["strains_to_drop"])
unknown_drops = set(drop_reasons) - set(strain_means["strain"])
if unknown_drops:
    raise ValueError(f"'strains_to_drop' not in the library: {sorted(unknown_drops)}")

no_counts_strains = [
    strain
    for strain in strain_means.query("mean_fraction_strain <= 0")["strain"]
    if strain not in drop_reasons
]
for strain in no_counts_strains:
    drop_reasons[strain] = NO_COUNTS_REASON
if no_counts_strains:
    add_warning(
        f"{len(no_counts_strains)} strains have no counts and so cannot be re-pooled, "
        f"but are not in `strains_to_drop` for re-pool `{repool}` in `config.yml`. They "
        "have been dropped anyway, and should be added there with the reason each one "
        f"failed: {', '.join(f'`{s}`' for s in no_counts_strains)}"
    )

kept = ~strain_means["strain"].isin(drop_reasons)

unknown_remakes = set(subpools_to_remake) - set(strain_means["subpool"])
if unknown_remakes:
    raise ValueError(
        f"'remake' is set for subpools {sorted(unknown_remakes)} of re-pool {repool}, "
        "which no strain belongs to"
    )

strain_means["remake_subpool"] = strain_means["subpool"].isin(subpools_to_remake)

# Stage 1: what each subpool will be made of.
#
# A remade subpool is built again from the strain stocks, at volumes that equalize its
# strains, so its composition comes out flat. A subpool left alone keeps whatever
# composition it already has, which is the strains' measured shares renormalized within it.
# Either way `subpool_composition` is a strain's share *of its own subpool*, summing to one
# within each.
_share = (
    strain_means["mean_fraction_strain"]
    .where(~strain_means["remake_subpool"], 1.0)
    .where(kept)
)
strain_means["subpool_composition"] = _share / _share.groupby(
    strain_means["subpool"]
).transform("sum")

# Stage 2: how much of each subpool goes into the final pool. Every strain is aiming for
# the same share, so a subpool contributes in proportion to how many strains it holds.
subpool_strains = strain_means[kept].groupby("subpool")["strain"].nunique()
subpool_fraction_of_pool = subpool_strains / subpool_strains.sum()

# A strain's share of the *final* pool is then its share of its subpool times that
# subpool's share of the pool.
strain_means["predicted_fraction_strain"] = strain_means[
    "subpool_composition"
] * strain_means["subpool"].map(subpool_fraction_of_pool)

# Volumes to pipette, for the subpools being remade. `ratio_to_add` is the volume of a
# strain relative to a typical strain *of its own subpool*, normalized to geometric mean
# one there, so the numbers stay comparable within the subpool actually being made.
#
# Float NaN rather than `pd.NA`: a column initialized with `pd.NA` is `object` dtype, which
# silently survives the arithmetic below and then breaks `nlargest` further down.
strain_means["ratio_to_add"] = float("nan")
strain_means["neat_volume_uL"] = float("nan")
strain_means["dilution_factor"] = float("nan")
strain_means["volume_to_add_uL"] = float("nan")


def dilution_for(neat_volume):
    """Smallest configured dilution bringing ``neat_volume`` up to the minimum, or 1."""
    if neat_volume >= min_pipettable_volume:
        return 1
    for step in dilution_steps:
        if neat_volume * step >= min_pipettable_volume:
            return step
    return dilution_steps[-1]  # nothing is enough; warned about below


# how much more liquid each remade subpool ends up holding because of those dilutions
subpool_dilution_bulking = {}

for _subpool in subpools_to_remake:
    _in = kept & (strain_means["subpool"] == _subpool)
    _vol = (
        strain_means.loc[_in, "previous_volume_uL"]
        / strain_means.loc[_in, "mean_fraction_strain"]
    )
    _ratio = _vol / math.exp(_vol.map(math.log).mean())
    strain_means.loc[_in, "ratio_to_add"] = _ratio
    # scaled so this subpool makes the volume its share of the final pool calls for
    _subpool_volume = (
        repool_config["total_pool_volume"] * subpool_fraction_of_pool[_subpool]
    )
    _neat = _ratio * (_subpool_volume / _ratio.sum())
    strain_means.loc[_in, "neat_volume_uL"] = _neat

    # A strain needing too little of its stock to pipette is diluted instead, and that
    # much more of the dilution added. The virus delivered is identical, so the subpool's
    # composition does not change; only its volume grows, by the extra carrier liquid.
    _factor = _neat.map(dilution_for)
    strain_means.loc[_in, "dilution_factor"] = _factor
    strain_means.loc[_in, "volume_to_add_uL"] = _neat * _factor
    subpool_dilution_bulking[_subpool] = (_neat * _factor).sum() / _neat.sum()

    _still_short = _neat * _factor < min_pipettable_volume
    if _still_short.any():
        add_warning(
            f"{int(_still_short.sum())} strains of subpool `{_subpool}` still need less "
            f"than {min_pipettable_volume} uL after the largest configured dilution of "
            f"1:{dilution_steps[-1]}: "
            f"{', '.join(f'`{s}`' for s in strain_means.loc[_in][_still_short]['strain'])}. "
            "Add a larger step to `dilution_steps` in `config.yml`."
        )

# Combining the subpools: a remade subpool is more dilute than planned by exactly the
# extra liquid its dilutions added, so proportionally more of it has to go in for it to
# deliver the amount of virus its share calls for. Subpools that are not remade are
# unaffected, and the shares are renormalized so they still sum to one.
subpool_combining_fraction = subpool_fraction_of_pool * pd.Series(
    {s: subpool_dilution_bulking.get(s, 1.0) for s in subpool_fraction_of_pool.index}
)
subpool_combining_fraction /= subpool_combining_fraction.sum()

dropped = (
    strain_means[~kept]
    .assign(reason=lambda x: x["strain"].map(drop_reasons))[
        ["shortname", "strain", "subpool", "mean_fraction_strain", "reason"]
    ]
    .reset_index(drop=True)
)
repooling_math = strain_means[kept][
    [
        "shortname",
        "strain",
        "subpool",
        "remake_subpool",
        "mean_fraction_strain",
        "representation",
        "previous_volume_uL",
        "subpool_composition",
        "predicted_fraction_strain",
        "ratio_to_add",
        "neat_volume_uL",
        "dilution_factor",
        "volume_to_add_uL",
    ]
].reset_index(drop=True)
# Only the strains of a remade subpool get a volume; the rest are not being repipetted.
_remade = repooling_math["remake_subpool"]
assert repooling_math.loc[_remade, "volume_to_add_uL"].notna().all()
assert (repooling_math.loc[_remade, "volume_to_add_uL"] < float("inf")).all()
assert repooling_math.loc[~_remade, "volume_to_add_uL"].isna().all()
# Diluting changes how much liquid is added, never how much virus, so the volume of stock
# a strain contributes must be what the undiluted calculation asked for.
_stock = (
    repooling_math.loc[_remade, "volume_to_add_uL"]
    / repooling_math.loc[_remade, "dilution_factor"]
)
assert (
    (_stock - repooling_math.loc[_remade, "neat_volume_uL"]).abs() < 1e-9
).all(), "diluting a strain changed the amount of its stock going into the pool"

# Check the predicted composition is internally consistent: a strain's share of its subpool
# times its subpool's share of the pool, summed over everything, must be one.
_pred = repooling_math["predicted_fraction_strain"].astype(float)
assert (
    abs(_pred.sum() - 1) < 1e-9
), f"the predicted composition sums to {_pred.sum():.6f} rather than one"

# For a subpool being remade, forward-simulate the volumes actually being pipetted and
# check they flatten it. Mixing volume `V` of a stock of titer `t` gives shares
# proportional to `V x t`, and `t` is `fraction_strain / previous_volume`.
#
# Be clear about what this does and does not prove. The volumes are built from the same two
# columns the implied titer is built from, so the arithmetic reduces to a constant and the
# check passes for any input data. It is a self-consistency check on this script's
# arithmetic, not evidence that the premise behind it holds. It does catch the mistake worth
# catching -- using `1 / fraction_strain`, the rule that is right for the equal-volume pool
# `analyze_pool` measures -- which blows the spread up by orders of magnitude while the
# volumes still sum correctly with geometric mean one.
#
# What it cannot detect, because nothing in the data reveals any of them: stock titers that
# drifted between the two poolings, `previous_volume_uL` differing from what was actually
# pipetted, or a subpool being mixed from something other than the strain stocks.
for _subpool in subpools_to_remake:
    _r = repooling_math[repooling_math["subpool"] == _subpool]
    _titer = _r["mean_fraction_strain"] / _r["previous_volume_uL"]
    # the neat volume, since that is the amount of stock delivered whether or not the
    # strain was diluted first
    _p = _r["neat_volume_uL"].astype(float) * _titer
    _p = _p / _p.sum()
    _spread = _p.max() / _p.min()
    print(f"Remade subpool {_subpool}: predicted spread {_spread:.4f}x")
    assert _spread < 1.001, (
        f"the volumes for subpool {_subpool} do not flatten it: predicted representation "
        f"spans {_spread:.3f}x, expected 1.000x"
    )

_predicted_spread = _pred.max() / _pred.min()
_current_spread = (
    strain_means.loc[kept, "mean_fraction_strain"].max()
    / strain_means.loc[kept, "mean_fraction_strain"].min()
)
print(
    f"Whole pool: {_current_spread:.1f}x now, {_predicted_spread:.2f}x predicted after "
    "the corrective re-pool"
)

_remade_names = ", ".join(f"`{s}`" for s in subpools_to_remake) or "none"
add_markdown(f"""
    ### What this gives

    Subpools being remade from strain stocks: {_remade_names}. Every subpool is then
    combined in the proportion below.

    Across the {len(repooling_math)} strains that are kept (see below for the
    {len(dropped)} that are dropped), the most and least abundant strain currently differ
    **{_current_spread:.1f}-fold**. After the corrective re-pool they are predicted to
    differ **{_predicted_spread:.2f}-fold**.

    The volumes to pipette are in
    `{Path(snakemake.output.subpool_repooling_math).name}/`: one CSV per subpool being
    remade, most volume first, where `volume_to_add_uL` is what to pipette and
    `dilution_factor` how far the stock is diluted first, so a strain with a factor of 1
    is added neat; plus `combine_subpools.csv`, holding the fractions the finished
    subpools are combined in. The measurement those came from is in
    `{Path(snakemake.output.repooling_math).name}`, which also carries
    `predicted_fraction_strain`, what each strain is expected to come out at. The strains
    dropped from the re-pooling are in
    `{Path(snakemake.output.dropped_strains).name}`.

    **What this assumes.** The predictions hold only if each remade subpool is mixed from
    the individual strain stocks, those stocks have not changed titer since the previous
    pool was made, and the volumes recorded for that pool are what was actually pipetted.
    None of those can be checked from the sequencing data, so they are premises rather than
    findings. If any does not hold the numbers will be wrong in a way that looks perfectly
    reasonable, and the corrective pool should be sequenced again to find out. Nothing here
    predicts the *titer* of the corrective pool either, so it has to be titrated rather
    than diluted by calculation.
    """)

add_markdown(f"""
    ### Combining the subpools

    **`fraction_of_pool` is the number to work from**: multiply it by whatever volume of
    pool is being made. The fractions sum to one, so a 100 uL test pool is
    {", ".join(
        f"{100 * subpool_combining_fraction[_s]:.2f} uL of `{_s}`"
        for _s in sorted(subpool_combining_fraction.index)
    )}. `volume_uL` is the same thing for the full
    {repool_config["total_pool_volume"]:.0f} uL.

    The fraction already accounts for the dilutions below. A subpool holding diluted
    strains carries proportionally more liquid for the same amount of virus, so it goes in
    at proportionally greater volume; `equal_share_fraction` is what the fraction would be
    without that adjustment, and is shown only so the difference is visible.
    """)
subpool_plan = (
    repooling_math.groupby("subpool", as_index=False)
    .aggregate(
        n_strains=pd.NamedAgg("strain", "nunique"),
        remake=pd.NamedAgg("remake_subpool", "first"),
        current_fraction_of_pool=pd.NamedAgg("mean_fraction_strain", "sum"),
    )
    .assign(
        equal_share_fraction=lambda x: x["subpool"].map(subpool_fraction_of_pool),
        fraction_of_pool=lambda x: x["subpool"].map(subpool_combining_fraction),
        volume_uL=lambda x: (
            x["fraction_of_pool"] * repool_config["total_pool_volume"]
        ),
    )
    .sort_values("subpool")
    .reset_index(drop=True)
)
add_table(subpool_plan)

if subpools_to_remake:
    for _subpool in subpools_to_remake:
        _r = repooling_math[repooling_math["subpool"] == _subpool]
        _diluted = _r[_r["dilution_factor"] > 1]
        add_markdown(f"""
            #### Remaking `{_subpool}`

            Made from {len(_r)} strains, giving
            {_r["volume_to_add_uL"].sum():.0f} uL of subpool.
            """)
        if len(_diluted):
            add_markdown(f"""
                {len(_diluted)} of them need less than {min_pipettable_volume} uL of stock,
                so they are diluted first and that much more of the dilution added. The
                virus delivered is the same either way; the subpool just ends up
                {subpool_dilution_bulking[_subpool]:.3f}x its undiluted volume, which is
                why proportionally more of it goes into the final pool above.
                """)
            add_table(
                _diluted.assign(
                    dilution=lambda x: "1:"
                    + x["dilution_factor"].astype(int).astype(str)
                )[
                    [
                        "shortname",
                        "strain",
                        "neat_volume_uL",
                        "dilution",
                        "volume_to_add_uL",
                    ]
                ].reset_index(
                    drop=True
                )
            )
        else:
            add_markdown(f"""
                Every strain needs at least {min_pipettable_volume} uL of stock, so none
                has to be diluted first.
                """)
        add_markdown("The 5 strains needing the most volume:")
        add_table(
            _r.nlargest(5, "volume_to_add_uL")[
                ["shortname", "strain", "mean_fraction_strain", "volume_to_add_uL"]
            ].reset_index(drop=True)
        )
else:
    add_markdown("""
        No subpool is being remade, so the correction is entirely in the proportions
        above. Set `remake` for a subpool in `config.yml` to also rebalance the strains
        within it.
        """)

if len(dropped):
    add_markdown("""
        ### Dropped strains

        These strains are excluded from the corrective re-pooling for the reasons given.
        """)
    add_table(dropped)
else:
    add_markdown("No strains are dropped from the corrective re-pooling.")

# --- write the outputs -----------------------------------------------------------------

if warnings:
    report[warnings_index] = markdown.markdown(
        "## Warnings\n\n"
        + "\n".join(f" - {warning}" for warning in warnings)
        + "\n\nThese do not stop the analysis, but should be looked into."
    )
    print(f"Finished with {len(warnings)} warnings")

# written to four significant digits, which is far finer than anything can be pipetted
# and keeps a negligible change in the numbers from rewriting every line of the CSVs
repooling_math.to_csv(snakemake.output.repooling_math, index=False, float_format="%.4g")
dropped.to_csv(snakemake.output.dropped_strains, index=False, float_format="%.4g")

# One CSV per subpool being remade, holding just what is needed to pipette it: the strain,
# the volume, and the ratio that volume came from. The whole-pool CSV above keeps the
# measurement the volumes were derived from.
_subpool_dir = Path(snakemake.output.subpool_repooling_math)
_subpool_dir.mkdir(parents=True, exist_ok=True)

# how the finished subpools are combined, as fractions so any volume of pool can be made
subpool_plan.to_csv(
    _subpool_dir / "combine_subpools.csv", index=False, float_format="%.4g"
)
print(f"Wrote the subpool combining fractions to {_subpool_dir}/combine_subpools.csv")
for _subpool in subpools_to_remake:
    _out = _subpool_dir / f"{_subpool}_repooling_math.csv"
    (
        repooling_math[repooling_math["subpool"] == _subpool]
        .assign(
            dilution_factor=lambda x: x["dilution_factor"].astype(int),
        )[
            [
                "shortname",
                "strain",
                "subpool",
                "ratio_to_add",
                "neat_volume_uL",
                "dilution_factor",
                "volume_to_add_uL",
            ]
        ]
        .sort_values("volume_to_add_uL", ascending=False)
        .reset_index(drop=True)
        .to_csv(_out, index=False, float_format="%.4g")
    )
    print(f"Wrote the re-pooling volumes for {_subpool} to {_out}")

Path(snakemake.output.html).write_text(
    PAGE_TEMPLATE.format(
        title=f"Balance of library re-pool {repool}",
        body="\n".join(report),
    )
)

print(f"Wrote the report to {snakemake.output.html}")
