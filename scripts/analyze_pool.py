"""Analyze the composition of a library pool and compute the re-pooling volumes.

Writes a self-contained HTML report, a CSV with the volume of each strain to add to
re-pool the library so that all strains are equally represented, and a CSV of the
strains dropped from that re-pooling. Everything specific to a pool comes from the
`analyze_pools` section of `config.yml` and is passed in via `snakemake.params`, so this
script is not specific to any one pool.

The dilution predicted for the re-pooled library assumes that the pool being analyzed was
made by mixing **equal volumes** of each strain's stock. That is what makes a strain's
share of the reads proportional to the titer of its stock, and so what lets the re-pool's
titer be computed as the volume-weighted mean of those shares. A pool mixed any other way
breaks that step, and nothing in the data gives it away: the script would report a
confidently wrong dilution. Check how a pool was made before adding it to `analyze_pools`.

Three ordering constraints, none of them obvious from any one section on its own:

  - The wells the composition is read from are not known until the linear range has been
    located, so that section has to come before the invalid barcodes, the composition of
    the pool, and everything downstream of them.
  - The strains to drop are determined before the ratios and volumes, because the volumes
    are scaled so that the strains that are *kept* sum to `total_pool_volume`.
  - The report is built by appending to `report` in order, and one slot in it is reserved
    by index for the warnings, so nothing may be inserted ahead of that slot.

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
# falls along a line of slope IDEAL_SLOPE against the dilution factor on log-log axes for
# as long as that holds. IDEAL_SLOPE is what the assay should do rather than a choice, but
# the rest are judgement calls, inherited from the `flu-seqneut-fhCART` analysis of the same
# kind of plate rather than tuned here: where the pool starts responding is judged over runs
# of ONSET_WINDOW_DILUTIONS consecutive dilutions, and a run's fitted slope counts as ideal
# when it is within LINEARITY_TOL, widened by LINEARITY_Z standard errors of the slope so
# that a noisy fit is not called nonlinear on the strength of its noise. The runs are short
# and local on purpose, as judging by a window running to the most dilute end of the series
# averages the wells that behave with those that do not.
#
# These values decide which dilution is chosen, so the report gives the slope of every run
# and not just the chosen one; set `linear_range_wells` to fix the wells where the choice
# they make looks wrong.
IDEAL_SLOPE = -1.0
LINEARITY_TOL = 0.3
LINEARITY_Z = 3.0
ONSET_WINDOW_DILUTIONS = 3

# factor by which the odds axis is padded beyond the wells plotted on it
ODDS_AXIS_PAD = 1.5

# An invalid barcode within this Hamming distance of a known barcode is read as sequencing
# error off that barcode, and one further away as a genuinely different barcode, meaning
# material that does not belong in the pool rather than noise.
MAX_SEQUENCING_ERROR_HAMMING = 1

# The column of the invalid-barcode CSVs holding that distance. `bacode` is a typo in
# `seqneut-pipeline`, not here; if a later version of the pipeline fixes it, this is the
# one place to change.
HAMMING_COL = "closest_valid_bacode_hamming_distance"

# reason recorded for strains that cannot be re-pooled as they had no counts
NO_COUNTS_REASON = "no counts in the analyzed wells"

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

pool = snakemake.wildcards.pool
date = snakemake.params.date
pool_config = snakemake.params.pool_config

required_config_keys = {
    "miscellaneous_plate",
    "description",
    "linear_range_wells",
    "min_avg_barcode_count_per_well",
    "min_neut_standard_frac_per_well",
    "total_pool_volume",
    "subpools",
    "strains_to_drop",
}
if set(pool_config) != required_config_keys:
    raise ValueError(
        f"configuration for pool {pool} must have exactly the keys "
        f"{sorted(required_config_keys)}, but has {sorted(pool_config)}"
    )

# wells fixed in the configuration, or None to detect them from the dilution series below
configured_wells = pool_config["linear_range_wells"]
if configured_wells == "calculate":
    configured_wells = None
elif isinstance(configured_wells, str) or not configured_wells:
    raise ValueError(
        f"'linear_range_wells' for pool {pool} must be \"calculate\" or a list of wells, "
        f"but is {configured_wells}"
    )
elif len(configured_wells) != len(set(configured_wells)):
    raise ValueError(
        f"'linear_range_wells' for pool {pool} has duplicates: {configured_wells}"
    )

print(
    f"Analyzing pool {pool} from {date}, linear_range_wells={pool_config['linear_range_wells']}"
)

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
    rather than interpolating it into an indented literal: its unindented continuation
    lines make the common prefix empty, nothing is stripped, and markdown then renders the
    whole block as a code block. Indent a formula by four spaces beyond the rest so that
    four survive the dedent and it does render as one.

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
    # Composition of library pool `{pool}`

    Analysis of the barcode sequencing of the miscellaneous plate
    `{pool_config["miscellaneous_plate"]}` from {date}.
    The plots are interactive: mouse over points for details.

    ## Experimental description
    """)
add_markdown(pool_config["description"])  # added separately as it is multi-line

# Reserved by index and filled in at the very end, once every warning has been collected.
# Nothing may be inserted into `report` ahead of it, only appended after it.
warnings_index = len(report)
report.append("")

# --- read the input data -----------------------------------------------------------

samples = pd.read_csv(snakemake.input.samples_csv).drop(columns="fastq")
if "dilution_factor" not in samples.columns:
    raise ValueError(f"{snakemake.input.samples_csv} lacks a 'dilution_factor' column")
samples["sample"] = samples.astype(str).agg("-".join, axis=1)
assert samples["sample"].nunique() == len(samples), "samples are not uniquely named"

if configured_wells:
    missing_wells = set(configured_wells) - set(samples["well"])
    if missing_wells:
        raise ValueError(
            f"'linear_range_wells' for pool {pool} not on the plate: "
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


# --- fates of the sequencing reads ---------------------------------------------------

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

# --- barcode counts and their QC -----------------------------------------------------

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
    {pool_config["min_avg_barcode_count_per_well"]} counts per barcode or when fewer than
    {pool_config["min_neut_standard_frac_per_well"]} of their counts are from the
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
            x["avg_count"] < pool_config["min_avg_barcode_count_per_well"]
        ),
        # a null fraction cannot be compared, so treat it as failing
        fails_neut_standard_qc=lambda x: ~(
            x["neut_standard_frac"] >= pool_config["min_neut_standard_frac_per_well"]
        ),
        fails_qc=lambda x: x["fails_count_qc"] | x["fails_neut_standard_qc"],
        viral_count=lambda x: x["total_count"] - x["neut_standard_count"],
        # the neutralization standard's share of the counts in odds form, which is exactly
        # (neutralization standard) / (viral), and its reciprocal. Diluting the pool leaves
        # the neutralization standard alone while dividing the viral counts, so the
        # reciprocal falls along a line of slope IDEAL_SLOPE on log-log axes and the odds
        # rise along one of the opposite slope, for as long as the pool responds to being
        # diluted. The fraction itself just flattens towards one, which is why it is the
        # odds that are plotted and fit below.
        neut_standard_odds=lambda x: (
            x["neut_standard_count"] / x["viral_count"].where(x["viral_count"] > 0)
        ),
        viral_to_neut_standard=lambda x: (
            x["viral_count"]
            / x["neut_standard_count"].where(x["neut_standard_count"] > 0)
        ),
        # Only wells that pass the QC above are fit below. The counts also have to be
        # positive for the odds and their log to be defined at all, which is a matter of
        # arithmetic rather than a further threshold: a well with no viral counts has a
        # neutralization standard fraction of one, which passes any minimum.
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
                # few enough ticks that their labels cannot collide at this panel width
                axis=alt.Axis(tickCount=4),
            ),
            alt.Y(
                "sample",
                title=None,
                sort=sample_order,
                scale=alt.Scale(domain=sample_order),
                axis=alt.Axis(labelLimit=Y_LABEL_LIMIT) if first_panel else None,
            ),
            # `fill` rather than `color`, which for a bar can also set the stroke that
            # the hover highlight below uses
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
    # spaced apart so that the last tick label of one panel and the first of the next
    # cannot collide, which `labelOverlap` cannot prevent as it works within one axis
    alt.hconcat(*well_qc_panels, spacing=45)
    .add_params(well_hover)
    .resolve_scale(y="shared")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

# --- the linear range of the dilution series -----------------------------------------

add_markdown(f"""
    ## The linear range of the dilution series

    Diluting the pool should divide its viral counts by the dilution factor while leaving
    the neutralization standard alone, so the neutralization standard's share of the
    counts, in odds form, rises along a straight line of slope 1 against the dilution
    factor on log-log axes for as long as that holds. It is over the dilutions where that
    is true that the composition of the pool can be read, since a well too concentrated
    for its counts to fall with dilution measures the assay rather than the pool.

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
    """Fit the slope of each run of ``ONSET_WINDOW_DILUTIONS`` consecutive dilutions.

    The runs slide one dilution at a time and are returned from concentrated to dilute, so
    that each is judged on its own wells rather than on an average with the rest of the
    series.

    """
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
    """Take the first run of dilutions whose slope is within tolerance of ideal.

    Strict tolerance is preferred over the widening by the error on the slope, as a run
    admitted only by its error bars is not evidence that the pool responds to dilution.
    Returns ``None`` if no run qualifies either way.

    """
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
        f"no well of pool {pool} passes QC, so there is nothing to analyze"
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
# range, rather than from its most concentrated edge, so that there is a dilution of margin
# on either side of the wells used: at the edge the pool has only just started responding,
# and a slightly misplaced onset would put the wells outside the linear range altogether.
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
        f"cannot find the linear range for pool {pool}, so set 'linear_range_wells' to "
        'the wells to read the composition from rather than to "calculate"'
    )
print(f"Reading the composition from wells {analysis_wells}")

# Limits from the wells rather than from everything drawn, as the ideal line runs the
# length of the series and reaches odds the wells never take.
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
    # every well is drawn, but the limits come from the wells that are fit, so a failing
    # well far off the scale is clipped rather than compressing all the others
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
    # A line of the ideal slope rather than a fit to the wells: what the odds would do if
    # each dilution step divided the viral counts by that step. Only its height is fit, over
    # the run of dilutions that located the linear range, so it is drawn solid there, where
    # the wells are claimed to follow it, and faint elsewhere, where they are not.
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
        # the dilution the composition is read from
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
        {onset["deviation"]:.3f}. The solid line is the ideal slope placed over that run,
        drawn faint across the rest of the series.

        The composition is calculated from the **middle** dilution of that run,
        **{calculated_dilution}-fold** (wells
        {", ".join(f"`{w}`" for w in calculated_wells)}, marked by the dashed rule), rather
        than from its most concentrated edge. That leaves a dilution of margin on either
        side of the wells used: at the edge of the run the pool has only just started
        responding, so a slightly misplaced onset would put the wells outside the linear
        range altogether.
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
        f"none of the wells used for the re-pooling math for pool {pool} have any "
        f"barcode counts: {analysis_wells}"
    )
if analysis_wells_with_counts != analysis_wells:
    add_warning(
        "the composition below is measured from only "
        f"{analysis_wells_with_counts}, as the other wells chosen for it have no counts."
    )

# --- invalid barcodes ------------------------------------------------------------------

add_markdown(f"""
    ## Invalid barcodes

    Reads that parsed as a barcode but match neither the viral library nor the
    neutralization standard, over the wells the composition is read from. One within
    {MAX_SEQUENCING_ERROR_HAMMING} nucleotide of a known barcode is sequencing error off
    that barcode, while one further away is a genuinely different barcode, and so material
    that does not belong in the pool rather than noise. The distance to the closest known
    barcode is the one computed by `seqneut-pipeline`, not recomputed here; that barcode
    itself is named only where it is close enough to be what the reads came from.
    """)

invalid = read_per_well_csvs(snakemake.input.invalid, "_invalid")
missing_cols = {"barcode", "count", "closest_valid_barcode", HAMMING_COL} - set(
    invalid.columns
)
if missing_cols:
    raise ValueError(
        f"the invalid barcode CSVs lack the columns {sorted(missing_cols)}; a version of "
        "`seqneut-pipeline` that renames them needs this script updated"
    )
invalid = invalid[invalid["well"].isin(analysis_wells_with_counts)]

# the two outputs of the pipeline must agree on how many reads were invalid
invalid_per_well = invalid.groupby("well")["count"].sum()
fates_per_well = (
    fates[fates["fate"] == "invalid barcode"].set_index("well")["count"].to_dict()
)
for well in analysis_wells_with_counts:
    if invalid_per_well.get(well, 0) != fates_per_well.get(well, 0):
        raise ValueError(
            f"well {well} has {invalid_per_well.get(well, 0)} counts in its invalid "
            f"barcode CSV but {fates_per_well.get(well, 0)} reads with the 'invalid "
            "barcode' fate"
        )

valid_count = counts[counts["well"].isin(analysis_wells_with_counts)]["count"].sum()
invalid_count = invalid["count"].sum()

if not invalid_count:
    add_markdown("None of the reads in these wells had an invalid barcode.")
else:
    # sorted so that taking the first of each barcode takes its closest known barcode,
    # which can be a different one in each well when two are equally close
    invalid_barcodes = (
        invalid.sort_values(HAMMING_COL)
        .groupby("barcode", as_index=False)
        .aggregate(
            count=pd.NamedAgg("count", "sum"),
            hamming_distance=pd.NamedAgg(HAMMING_COL, "first"),
            closest_valid_barcode=pd.NamedAgg("closest_valid_barcode", "first"),
        )
        .assign(
            sequencing_error=lambda x: (
                x["hamming_distance"] <= MAX_SEQUENCING_ERROR_HAMMING
            )
        )
    )
    error_count = invalid_barcodes.query("sequencing_error")["count"].sum()

    add_markdown(f"""
        **{invalid_count / (valid_count + invalid_count):.2%}** of the
        {valid_count + invalid_count:,} parsed reads in wells
        {", ".join(f"`{w}`" for w in analysis_wells_with_counts)} had an invalid barcode,
        spread over {len(invalid_barcodes):,} distinct barcodes. Of those reads,
        {error_count / invalid_count:.1%} are within
        {MAX_SEQUENCING_ERROR_HAMMING} nucleotide of a known barcode and the remaining
        {1 - error_count / invalid_count:.1%} are novel.

        The most abundant of them:
        """)
    add_table(
        invalid_barcodes.nlargest(5, "count")
        .merge(
            barcode_class.rename(columns={"barcode": "closest_valid_barcode"}),
            on="closest_valid_barcode",
            how="left",
            validate="many_to_one",
        )
        .assign(
            fraction_of_reads=lambda x: x["count"] / (valid_count + invalid_count),
            barcode_is=lambda x: x["sequencing_error"].map(
                {True: "sequencing error", False: "novel"}
            ),
            # The closest known barcode is only reported for a barcode close enough to be
            # sequencing error off it. For a novel barcode the closest match is a
            # coincidence rather than what the reads came from, and naming it would imply
            # an association with that strain that the distance rules out.
            closest_is=lambda x: x["strain"]
            .where(~x["neut_standard"], "neutralization standard")
            .where(x["sequencing_error"], ""),
            closest_valid_barcode=lambda x: x["closest_valid_barcode"].where(
                x["sequencing_error"], ""
            ),
        )[
            [
                "barcode",
                "count",
                "fraction_of_reads",
                "barcode_is",
                "hamming_distance",
                "closest_valid_barcode",
                "closest_is",
            ]
        ]
    )

# --- representation of each strain in the pool ---------------------------------------

viral_counts = counts[counts["strain"].notna()]
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
    .assign(
        # the strain's share of the viral counts, and that share per barcode, which is
        # what is balanced as each barcode is a replicate in the neutralization assays
        fraction_strain_counts=lambda x: x["count"] / x["well_count"],
        fraction_strain=lambda x: x["fraction_strain_counts"] / x["n_barcodes"],
    )
)

# strain-level short names
strain_names = (
    viral_library.rename(columns={"shortname_strain": "shortname"})[
        ["strain", "shortname"]
    ]
    .drop_duplicates()
    .reset_index(drop=True)
)
if strain_names["strain"].nunique() != len(strain_names):
    raise ValueError("strains in the viral library lack a unique 'shortname'")


def assign_subpool(shortname):
    """Get the subpool for ``shortname``, the first one whose regex it matches."""
    for subpool, regex in pool_config["subpools"].items():
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
        mean_fraction_strain_counts=pd.NamedAgg("fraction_strain_counts", "mean"),
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
strain_order = strain_means["strain"].tolist()

# Strains excluded from the re-pooling: those named in the config, plus any that had no
# counts and so would need an infinite volume. Determined before the ratios and volumes
# below, which are computed over just the strains that are kept.
drop_reasons = dict(pool_config["strains_to_drop"])
unknown_drops = set(drop_reasons) - set(strain_means["strain"])
if unknown_drops:
    raise ValueError(f"'strains_to_drop' not in the library: {sorted(unknown_drops)}")

# A strain with no counts cannot be re-pooled, so it is dropped whether or not the
# configuration says to. That is not an error, but it does need to be acted on, as
# `strains_to_drop` should record why each strain failed.
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
        f"but are not in `strains_to_drop` for pool `{pool}` in `config.yml`. They have "
        "been dropped anyway, and should be added there with the reason each one "
        f"failed: {', '.join(f'`{s}`' for s in no_counts_strains)}"
    )

kept = ~strain_means["strain"].isin(drop_reasons)

# `ratio_to_add` is normalized so that its geometric mean over the kept strains is one,
# making it the volume of a strain relative to that of a typical strain.
#
# Its scale carries no information: the volumes below are scaled so the kept strains sum to
# `total_pool_volume`, so any constant factor in the ratio cancels there. Changing how it
# is normalized therefore changes the reported `ratio_to_add` column and nothing else, and
# the volumes are unaffected. This is also why it does not matter whether an equal share is
# taken over every measured strain or only the ones being re-pooled, which is a question
# that would otherwise be worth arguing about.
geometric_mean_fraction = math.exp(
    strain_means.loc[kept, "mean_fraction_strain"].map(math.log).mean()
)
strain_means["ratio_to_add"] = (
    geometric_mean_fraction / strain_means["mean_fraction_strain"]
)
strain_means["volume_to_add_uL"] = strain_means["ratio_to_add"] * (
    pool_config["total_pool_volume"] / strain_means.loc[kept, "ratio_to_add"].sum()
)

add_markdown(f"""
    ## Representation of each strain in the pool

    The left panel shows each strain's fraction of the viral counts (the neutralization
    standard is excluded) in wells
    {", ".join(f"`{w}`" for w in analysis_wells_with_counts)}, divided by the number of
    barcodes for that strain: the barcodes of a strain are pooled before it is rescued, so
    they cannot be balanced against each other. The right panel shows how a strain's
    counts are split among its barcodes, which should be roughly even unless a barcode is
    poorly represented in the rescued virus.
    """)

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

# hovering a bar in either panel outlines that strain in both; see `well_hover` above
strain_hover = alt.selection_point(
    name="strain_hover",
    fields=["strain"],
    on="pointerover",
    clear="pointerout",
    empty=False,
)

strain_panel = (
    alt.Chart(strain_means)
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
        stroke=alt.value("black"),
        strokeWidth=alt.condition(strain_hover, alt.value(2), alt.value(0)),
        tooltip=[
            "strain",
            "shortname",
            "subpool",
            alt.Tooltip("mean_fraction_strain", format=".3g"),
        ],
    )
    .mark_bar(height={"band": 0.85}, fill="gray")
    # the dashed line marks where every strain would be equally represented
    + alt.Chart(pd.DataFrame({"equal representation": [1 / n_strains]}))
    .mark_rule(strokeDash=[2, 2], strokeWidth=2, color="red")
    .encode(alt.X("equal representation"))
)

barcode_panel = (
    alt.Chart(barcode_fracs)
    .mark_bar(height={"band": 0.85})
    .encode(
        alt.X("fraction_barcode", title="fraction of strain from each barcode"),
        # y-axis labels come from the strain panel that shares this axis
        alt.Y(
            "strain",
            title=None,
            sort=strain_order,
            scale=alt.Scale(domain=strain_order),
            axis=None,
        ),
        alt.Fill("barcode_number:N", legend=None),
        stroke=alt.value("black"),
        strokeWidth=alt.condition(strain_hover, alt.value(2), alt.value(0)),
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
    .resolve_scale(y="shared")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

# --- re-pooling calculations ---------------------------------------------------------

add_markdown(f"""
    ## Re-pooling calculations

    Each strain is added at a volume proportional to the reciprocal of its current
    representation, so that all of the strains that are kept end up equally represented.
    That volume relative to a typical strain is `ratio_to_add`, normalized so its
    geometric mean over the kept strains is one: a strain with a `ratio_to_add` of 2 needs
    twice the volume of a typical strain, and one with 0.5 needs half. The volumes are
    then scaled to sum to the `total_pool_volume` of
    {pool_config["total_pool_volume"]} uL set in `config.yml`.
    """)

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
        "mean_fraction_strain",
        "ratio_to_add",
        "volume_to_add_uL",
    ]
].reset_index(drop=True)
assert repooling_math["volume_to_add_uL"].notna().all()
assert (repooling_math["volume_to_add_uL"] < float("inf")).all()

total_volume = repooling_math["volume_to_add_uL"].sum()

add_markdown(f"""
    Re-pooling the {len(repooling_math)} strains that are kept (see below for the
    {len(dropped)} that are dropped) gives **{total_volume:.0f} uL** of pool, adding
    between {repooling_math["volume_to_add_uL"].min():.2g} and
    {repooling_math["volume_to_add_uL"].max():.0f} uL of each strain.

    The volume of each strain to add is in
    `{Path(snakemake.output.repooling_math).name}`, and the strains dropped from the
    re-pooling are in `{Path(snakemake.output.dropped_strains).name}`.

    ### Dilution to use for the re-pooled library
    """)

# The chosen wells define the dilution of the current pool that gives the desired
# infection, so they must all be at the same dilution.
analysis_dilutions = set(
    samples[samples["well"].isin(analysis_wells_with_counts)]["dilution_factor"]
)
if len(analysis_dilutions) != 1:
    raise ValueError(
        f"the 'wells' for pool {pool} must share a 'dilution_factor', but have "
        f"{sorted(analysis_dilutions)}"
    )
current_dilution = analysis_dilutions.pop()

# The current pool was made by mixing equal volumes of each strain, so a strain's share
# of the counts is proportional to the titer of its stock, and the pool's titer is the
# mean of those stock titers. The re-pool mixes volume `V` of each strain instead, so its
# titer is the `V`-weighted mean of the same stock titers.
kept_means = strain_means[kept]
relative_titer = (
    n_strains
    * (kept_means["volume_to_add_uL"] * kept_means["mean_fraction_strain_counts"]).sum()
    / kept_means["volume_to_add_uL"].sum()
)
repool_dilution = current_dilution * relative_titer
print(f"Re-pool is {relative_titer:.4f} of the current pool's titer")

add_markdown(f"""
    Wells {", ".join(f"`{w}`" for w in analysis_wells_with_counts)} were chosen as being
    in the linear range, so a **{current_dilution:g}-fold** dilution of the current pool
    gives the amount of infection that the neutralization assays should use. The re-pooled
    library will not have the same titer, as balancing the strains means adding a lot of
    volume of the poorly growing ones.

    The current pool was made by mixing equal volumes of all {n_strains} strains, so each
    strain's share `g` of the counts is proportional to the titer of its stock, and the
    current pool's titer is the mean of those {n_strains} stock titers. The re-pool
    instead mixes volume `V` of each strain, so its titer is the `V`-weighted mean of the
    same stock titers, and the ratio of the two is:

        re-pool titer / current titer = {n_strains} x sum(V x g) / sum(V) = {relative_titer:.3f}

    So the re-pool is expected to be **{1 / relative_titer:.2f}-fold weaker**, and
    therefore needs proportionally less diluting to give the same infection:

        {current_dilution:g} x {relative_titer:.3f} = {repool_dilution:.0f}

    **Use the re-pooled library at about a {repool_dilution:.0f}-fold dilution.**

    This assumes that the counts are proportional to infectivity, which is why wells in
    the linear range are used, and that the strain stocks have not changed titer since the
    current pool was made. It is a starting point rather than a substitute for titrating
    the re-pooled library.
    """)

add_markdown("""
    ### Volume of each subpool

    The strains are combined in subpools that are then combined into the final pool, so
    that a problem with one subpool does not require remaking all of them.
    """)

subpool_summary = (
    repooling_math.groupby("subpool", as_index=False)
    .aggregate(
        n_strains=pd.NamedAgg("strain", "nunique"),
        volume_uL=pd.NamedAgg("volume_to_add_uL", "sum"),
    )
    .assign(fraction_of_pool=lambda x: x["volume_uL"] / x["volume_uL"].sum())
)
add_table(subpool_summary)

add_markdown("""
    ### Strains needing the most and least volume

    The strains needing the most volume are the poorest growing ones that are kept, and
    they dominate the volume of the pool. Those needing the least are the ones already
    best represented in the current pool.
    """)
add_markdown("The 5 strains needing the most volume:")
add_table(repooling_math.nlargest(5, "volume_to_add_uL"))
add_markdown("The 5 strains needing the least volume:")
add_table(repooling_math.nsmallest(5, "volume_to_add_uL"))

add_markdown("""
    ### Dropped strains

    These strains are excluded from the re-pooling for the reasons given.
    """)
add_table(dropped)

# --- write the outputs ---------------------------------------------------------------

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
Path(snakemake.output.html).write_text(
    PAGE_TEMPLATE.format(
        title=f"Composition of library pool {pool}",
        body="\n".join(report),
    )
)

print(f"Wrote the report to {snakemake.output.html}")
