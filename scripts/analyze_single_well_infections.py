"""Analyze the plates that infect each well with a single virus of the library.

Writes a self-contained HTML report and a CSV of what every well held. These plates
measure two things at once: how well each strain grows on its own, and whether virus from
one well ends up in another.

A well should hold the strain named for it and the neutralization standard spiked into it,
and nothing else. Every parsed read of a well is therefore put into one of five categories
that partition it:

    own strain          a barcode of the strain the samples CSV assigns to the well
    neut standard       a barcode of the neutralization standard
    other strains       a barcode of a *different* library strain, which does not belong
    unmatched, near     matched neither whitelist, but is within
                        `MAX_OWN_STRAIN_ERROR_HAMMING` of a barcode of the well's own
                        strain, so it is that strain's own material read imperfectly
    unmatched, far      matched neither whitelist and is further than that from the own
                        strain, so it cannot be attributed to the strain that belongs there

The titer of a strain is its own reads against the neutralization standard in its well,
`own / neut standard`. Only those two are in that ratio, so a contaminated well's titer is
not inflated by whatever else was in the well; that is reported separately instead.

Because `well` is the physical position on the plate, contamination can be traced rather
than only counted: a read of another library strain is attributed to the well that strain
was grown in and to how far away that well is, counting a step onto a diagonal neighbor as
one. Material from a neighboring well is spillover; material from far away, or from a
strain grown on another plate of the group, arrived some other way.

Two things to hold on to when reading the report:

  - A well whose own strain barely grew shows a large foreign *fraction* from a small
    absolute spill, so every fraction is reported beside the counts behind it.
  - A strain's reads are summed over its barcodes, of which strains have differing numbers,
    so a strain with more barcodes contributes more reads at the same titer.

One group of plates is analyzed together and reported in one HTML. The barcodes of all of
them are counted against one whitelist, the group's `viral_library`; the separate
`final_viral_library` is read only to mark which of the strains assayed here were carried
forward, and is never used in a calculation. Everything specific to a group comes from the
`analyze_single_well_infections` section of `config.yml` via `snakemake.params`, so this
script is not specific to any one group.

"""

import sys
import textwrap
from html import escape
from pathlib import Path

import altair as alt
import markdown
import numpy as np
import pandas as pd

sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

alt.data_transformers.disable_max_rows()

# pixels allowed for axis labels, enough for the longest strain name
LABEL_LIMIT = 300

# The plate layout that a well name is read as, so that `B6` becomes a row and a column and
# two wells can be told how far apart they are. A 96-well plate, which is what these
# infections are done in; a well outside it raises rather than being placed somewhere
# arbitrary.
PLATE_ROWS = "ABCDEFGH"
PLATE_COLUMNS = 12

# An unmatched barcode this close to one of the barcodes of the strain in its own well is
# that strain's own material misread rather than something that does not belong. A distance
# to the barcodes of *one strain*, which is a different question from the distance to the
# whole whitelist below, so the two are separate thresholds even at the same value.
MAX_OWN_STRAIN_ERROR_HAMMING = 1

# An unmatched barcode this close to a barcode of either whitelist is read as sequencing
# error off that barcode, and one further away as a genuinely different barcode, meaning
# material that does not belong on the plate at all rather than noise.
MAX_SEQUENCING_ERROR_HAMMING = 1

# The column of the invalid-barcode CSVs holding that distance. `bacode` is a typo in
# `seqneut-pipeline`, not here; if a later version of the pipeline fixes it, this is the
# one place to change.
HAMMING_COL = "closest_valid_bacode_hamming_distance"

# Distances further than this from a well are reported as one group rather than each on its
# own, there being nothing to tell them apart: what matters is whether material came from a
# neighboring well, not whether it came from five wells away or from nine.
MAX_NAMED_DISTANCE = 2

# Other strains and unmatched barcodes listed for each well the two reporting thresholds
# pick out. A display limit, affecting only how much of a well is broken out and never what
# is computed, so it lives here rather than in the configuration.
N_TOP_BARCODES = 3

# Parsed reads a well needs for the fractions in this report to rest on anything. Only warns,
# and so changes nothing that is computed or shown, which is why it lives here rather than in
# the configuration. No average over barcodes, unlike the pool analyses: a well is expected
# to hold one strain's barcodes rather than the whole library.
MIN_PARSED_READS_PER_WELL = 1000

# Neutralization standard reads a well needs for its titer to be measured with any precision,
# the standard's reads being that ratio's denominator. Used only to flag titers in the last
# section, which is the one place the amount of standard matters quantitatively, and only to
# warn, so it is a constant rather than configuration. The standard's *share* of a well is not
# checked anywhere: with one strain per well that share is `1 / (1 + titer)`, so it carries no
# information the titer does not.
MIN_NEUT_STANDARD_READS_FOR_TITER = 100

# How a strain assayed here is marked wherever it appears: by the color of the stroke
# around its bars, and by the color of its text in the tables. The dropped color has to read
# as body text as well as a thin outline, and to stay clear of the three fill palettes the
# panels use, so it is neither pale nor a member of `tableau10`, `dark2`, or `set1`.
CARRIED_FORWARD_LABELS = {True: "carried forward", False: "dropped"}
RETAINED_COLOR = "black"
DROPPED_COLOR = "#cc00cc"

# Width of that stroke, ordinarily and when the strain is hovered. Nonzero either way, the
# stroke being what carries the mark. Thin around a bar, where it would otherwise crowd the
# divisions of a stacked one, and heavier around a point, where a hairline would vanish.
STROKE_WIDTH = 0.5
STROKE_WIDTH_HOVERED = 2.5
POINT_STROKE_WIDTH = 1.5
POINT_STROKE_WIDTH_HOVERED = 3

# What each fate of a read means. Only the fates that have counts are explained in the
# report, and a fate the pipeline writes that is missing here raises: a fate nobody can
# interpret is worse than one nobody sees.
FATE_EXPLANATIONS = {
    "valid barcode": (
        "the barcode was read and is in the viral library or the neutralization standard "
        "set, so the read is counted below"
    ),
    "invalid barcode": (
        "the barcode was read but is in neither, so the read counts towards neither "
        "strain nor standard; these are the reads the unmatched barcodes section is about"
    ),
    "low quality barcode": (
        "the barcode was found but some of its nucleotides were sequenced too poorly to "
        "be trusted"
    ),
    "unparseable barcode": (
        "the sequence flanking the barcode was found but the barcode itself could not be "
        "read from the position it should be in"
    ),
    "invalid outer flank": (
        "the outer flanking sequence did not match. Where an outer barcode is used to "
        "separate samples that share a sequencing run, reads belonging to the other "
        "samples are what fails this way, so some of these reads are expected rather than "
        "a problem"
    ),
    "read too short": "the read was too short to contain the barcode and its flanks",
    "failed chastity filter": "the read failed the Illumina chastity filter",
}

# The four categories that partition the virus in a well: the column each is computed in,
# mapped to the label it carries in the report, in the order they stack. Ordered here rather
# than left to sort by name, which would separate the two unmatched categories.
#
# The neutralization standard is not among them, and every fraction in this report is a
# fraction of the *non-neut-standard* parsed reads, named
# `NON_NEUT_STANDARD_READS` below. The standard is a fixed spike rather than part of what a
# well grew, so including it would make each of these fractions depend on how much virus the
# well happened to hold. What the standard is for is the titer, in the last section.
COMPOSITION_CATEGORIES = {
    "own_strain_reads": "own strain",
    "other_strain_reads": "other library strains",
    "unmatched_near_expected_reads": (
        f"within {MAX_OWN_STRAIN_ERROR_HAMMING} nt of strain's or neut standard barcode"
    ),
    "unmatched_far_expected_reads": (
        f"over {MAX_OWN_STRAIN_ERROR_HAMMING} nt from strain's or neut standard barcode"
    ),
}

# the non-neut-standard parsed reads the categories above are taken as fractions of, spelled
# with hyphens as it is a column of the output CSV as well as of the frame here
NON_NEUT_STANDARD_READS = "non-neut-standard_parsed_reads"

# Of those categories, the one that is the well's own strain, which is to say the only one
# that belongs there. Everything that does not belong is the rest, named by exclusion so that
# a category added above joins that sum rather than being silently left out of it.
BELONGS_CATEGORIES = ["own_strain_reads"]

# and, of what does not belong, the part that is not the own strain's material at all: the
# same thing without the unmatched reads close enough to be its barcodes misread. This is
# what the wells are ordered by below.
EXPECTED_MISREAD_CATEGORY = "unmatched_near_expected_reads"

PAGE_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: sans-serif; max-width: 55em; margin: 2em auto; padding: 0 1em; }}
  table {{ border-collapse: collapse; margin: 1em 0; font-size: 90%; }}
  th, td {{ border: 1px solid #ccc; padding: 0.2em 0.6em; text-align: left; }}
  /* the last column of these tables packs several findings into one cell */
  table.small-last-column td:last-child {{ font-size: 80%; }}
  /* rows for strains not carried forward, marked as their bars are */
  span.dropped {{ color: {dropped_color}; }}
  h1, h2 {{ border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }}
</style>
</head>
<body>
{body}
</body>
</html>
"""

# --- read and validate the configuration -----------------------------------------------

single_well = snakemake.wildcards.single_well
plate_dates = snakemake.params.plate_dates
single_well_config = snakemake.params.single_well_config

required_config_keys = {
    "miscellaneous_plates",
    "description",
    "final_viral_library",
    "max_other_strain_frac",
    "max_unmatched_frac",
}
if set(single_well_config) != required_config_keys:
    raise ValueError(
        f"configuration for {single_well} must have exactly the keys "
        f"{sorted(required_config_keys)}, but has {sorted(single_well_config)}"
    )

plates = list(single_well_config["miscellaneous_plates"])
if len(plates) != len(set(plates)):
    raise ValueError(
        f"'miscellaneous_plates' for {single_well} repeats a plate: {plates}"
    )
if set(plates) != set(plate_dates):
    raise ValueError(
        f"the plates {plates} of {single_well} do not match the plates the rule passed "
        f"dates for, {sorted(plate_dates)}"
    )

print(f"Analyzing single-well infections {single_well} over plates {plates}")

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
    whole block as a code block.

    """
    report.append(
        markdown.markdown(textwrap.dedent(text).strip(), extensions=["tables"])
    )


def add_chart(chart):
    """Add an Altair ``chart`` to the report."""
    report.append(chart.to_html(fullhtml=False, output_div=f"chart{len(report)}"))


def add_table(df, small_last_column=False, dropped=None):
    """Add the data frame ``df`` to the report as a table.

    ``small_last_column`` shrinks the type in the last column, for the tables whose last
    column packs several findings into one cell and so runs much longer than the rest.

    ``dropped`` is a boolean mask of rows whose strain was not carried forward, whose text is
    written in `DROPPED_COLOR` to mark them as their bars are marked. Every cell is rendered
    to text and escaped here rather than by ``to_html``, which cannot be left to escape the
    markup that carries that color.

    """
    text = df.map(
        lambda v: (
            "" if pd.isna(v) else escape(f"{v:.4g}" if isinstance(v, float) else str(v))
        )
    )
    if dropped is not None and dropped.any():
        text.loc[dropped] = text.loc[dropped].map(
            lambda v: f'<span class="dropped">{v}</span>'
        )
    report.append(
        text.to_html(
            index=False,
            escape=False,
            classes="small-last-column" if small_last_column else None,
        )
    )


add_markdown(f"""
    # Single virus per well infections `{single_well}`

    Analysis of the barcode sequencing of
    {", ".join(f"`{plate}` ({plate_dates[plate]})" for plate in plates)}, whose wells were
    each infected with a single virus of the library rather than with the pool. The plots
    are interactive: mouse over a bar for details.

    ## Experimental description
    """)
add_markdown(single_well_config["description"])  # added separately as it is multi-line

# Reserved by index and filled in at the very end, once every warning has been collected.
# Nothing may be inserted into `report` ahead of it, only appended after it.
warnings_index = len(report)
report.append("")

# --- read the input data ---------------------------------------------------------------


def well_position(well):
    """Get the (row index, column) of ``well``, a well name such as ``B6``."""
    row, column = well[:1], well[1:]
    if row not in PLATE_ROWS or not column.isdigit():
        raise ValueError(
            f"well '{well}' is not a row letter followed by a column number, which a "
            f"{len(PLATE_ROWS)} x {PLATE_COLUMNS} plate needs it to be"
        )
    if not 1 <= int(column) <= PLATE_COLUMNS:
        raise ValueError(f"well '{well}' is outside a {PLATE_COLUMNS}-column plate")
    return PLATE_ROWS.index(row), int(column)


def wells_apart(well_1, well_2):
    """Get how many wells apart ``well_1`` and ``well_2`` are on the plate.

    A step onto any of the eight surrounding wells counts as one, diagonals included, so
    two wells one apart are the wells that touch. That is the right measure for liquid
    carried between wells, which has no reason to prefer a row to a column or either to a
    diagonal.

    """
    (row_1, column_1), (row_2, column_2) = well_position(well_1), well_position(well_2)
    return max(abs(row_1 - row_2), abs(column_1 - column_2))


def plate_and_well(csv_file, suffix):
    """Get the (plate, well) that a per-well ``suffix`` CSV is for, from its path."""
    path = Path(csv_file)
    return path.parent.name, path.stem.removesuffix(suffix)


def read_per_well_csvs(csv_files, suffix, required_columns):
    """Read per-well CSVs of every plate into one data frame labeled by plate and well.

    The wells are checked from the file *paths* rather than from the rows read, since a CSV
    of this kind can legitimately be empty: a well with no unmatched barcodes has nothing
    to list, and that is not the same as a well nobody counted.

    ``required_columns`` is checked on each file as it is read rather than on the frame they
    are concatenated into. A column missing from one file and present in the rest survives
    the concatenation as nulls, and would then be read as a measurement of zero rather than
    as the absence of one.

    """
    read_for = [plate_and_well(csv_file, suffix) for csv_file in csv_files]
    if len(set(read_for)) != len(read_for):
        raise ValueError(f"two of the '{suffix}' CSVs are for the same plate and well")
    expected = set(map(tuple, samples[["plate", "well"]].values))
    if set(read_for) != expected:
        raise ValueError(
            f"the '{suffix}' CSVs are for wells the samples CSVs do not list, or the "
            f"other way round: {sorted(set(read_for) ^ expected)}"
        )
    frames = []
    for csv_file, (plate, well) in zip(csv_files, read_for, strict=True):
        frame = pd.read_csv(csv_file)
        missing = set(required_columns) - set(frame.columns)
        if missing:
            raise ValueError(
                f"{csv_file} lacks the columns {sorted(missing)}; a version of "
                "`seqneut-pipeline` that renames them needs this script updated"
            )
        frames.append(frame.assign(plate=plate, well=well))
    return pd.concat(frames, ignore_index=True)


samples = pd.concat(
    [
        pd.read_csv(csv_file).assign(plate=plate)
        for (plate, csv_file) in zip(plates, snakemake.input.samples_csvs, strict=True)
    ],
    ignore_index=True,
)
missing_columns = {"well", "strain", "dilution_factor"} - set(samples.columns)
if missing_columns:
    raise ValueError(
        f"the samples CSVs lack the columns {sorted(missing_columns)}, which name the "
        "plate position of each sample and the virus infected into it"
    )
samples = samples[["plate", "well", "strain", "dilution_factor"]]
samples[["row", "column"]] = [well_position(well) for well in samples["well"]]

if samples.duplicated(["plate", "well"]).any():
    raise ValueError("a plate has two samples in the same well")
if samples.duplicated(["plate", "strain"]).any():
    raise ValueError(
        "a plate has the same strain in more than one well, so a read of that strain "
        "cannot be traced back to the well it came from"
    )

# Where each well falls in the order the plate was filled, down each column and then across.
# The charts below all order their rows by what they are about rather than by position, so this
# is here to put the output CSV in plate order, which is the order a plate is read in.
samples = samples.sort_values(["plate", "column", "row"]).reset_index(drop=True)
samples["position"] = samples.groupby("plate").cumcount()

viral_library = pd.read_csv(snakemake.input.viral_library)
neut_standard_barcodes = pd.read_csv(snakemake.input.neut_standard_set)[
    "barcode"
].tolist()
barcode_class = pd.concat(
    [
        viral_library[["barcode", "strain"]].assign(neut_standard=False),
        pd.DataFrame({"barcode": neut_standard_barcodes}).assign(
            neut_standard=True,
            strain=pd.NA,
        ),
    ],
    ignore_index=True,
)

not_in_library = sorted(set(samples["strain"]) - set(viral_library["strain"]))
if not_in_library:
    raise ValueError(
        f"the samples CSVs assign strains to wells that are not in "
        f"{snakemake.input.viral_library}, so they have no barcodes to count and cannot be "
        f"analyzed: {not_in_library}"
    )

# barcodes of each strain, needed to measure how far an unmatched barcode is from the
# barcodes of the strain that is supposed to be in its well
strain_barcodes = viral_library.groupby("strain")["barcode"].apply(list).to_dict()
n_barcodes = viral_library.groupby("strain")["barcode"].nunique()

counts = read_per_well_csvs(snakemake.input.counts, "_counts", ["barcode", "count"])[
    ["plate", "well", "barcode", "count"]
]
if set(counts["barcode"]) != set(barcode_class["barcode"]):
    raise ValueError(
        "barcodes in the counts do not match those in the viral library and "
        "neutralization standard set"
    )
counts = counts.merge(barcode_class, on="barcode", validate="many_to_one").merge(
    samples[["plate", "well", "strain"]].rename(columns={"strain": "own_strain"}),
    on=["plate", "well"],
    validate="many_to_one",
)

fates = read_per_well_csvs(snakemake.input.fates, "_fates", ["fate", "count"])
unexplained_fates = sorted(
    set(fates[fates["count"] > 0]["fate"]) - set(FATE_EXPLANATIONS)
)
if unexplained_fates:
    raise ValueError(
        f"`seqneut-pipeline` reported the fates {unexplained_fates}, which this script has "
        "no explanation for; add them to `FATE_EXPLANATIONS`"
    )

invalid = read_per_well_csvs(
    snakemake.input.invalid,
    "_invalid",
    ["barcode", "count", "closest_valid_barcode", HAMMING_COL],
)

# --- which strains are assayed here ------------------------------------------------------

# Every titer below is a ratio to the same fixed spike, so two wells are comparable only if
# they were at the same dilution, and they are all ranked against each other in one chart.
_dilutions = sorted(samples["dilution_factor"].unique())
if len(_dilutions) != 1:
    raise ValueError(
        f"the wells of this group are at more than one dilution factor {_dilutions}, so "
        "their titers are not comparable as computed and would each have to be scaled by "
        "their own dilution factor before being ranked against one another"
    )

# How each strain of the counted library fared in the library carried forward, and whether it
# is assayed here at all. Computed from the strains and the two library CSVs alone, so it is
# available before anything is counted and can mark every bar and table row below. Over the
# whole counted library rather than the assayed strains alone, so that the ones not assayed can
# be reported with their status too.
final_viral_library = pd.read_csv(snakemake.input.final_viral_library)
final_n_barcodes = final_viral_library.groupby("strain")["barcode"].nunique()
strain_status = (
    viral_library[["strain"]]
    .drop_duplicates()
    .assign(
        assayed=lambda x: x["strain"].isin(set(samples["strain"])),
        n_barcodes=lambda x: x["strain"].map(n_barcodes),
        n_barcodes_final=lambda x: x["strain"]
        .map(final_n_barcodes)
        .fillna(0)
        .astype(int),
    )
    .assign(
        final_library=lambda x: np.where(
            x["n_barcodes_final"] == 0,
            "dropped from the library",
            np.where(
                x["n_barcodes_final"] < x["n_barcodes"],
                "kept, some barcodes dropped",
                "kept",
            ),
        ),
        carried_forward=lambda x: (x["n_barcodes_final"] > 0).map(
            CARRIED_FORWARD_LABELS
        ),
    )
    .reset_index(drop=True)
)

_counted_name = Path(snakemake.input.viral_library).name
_final_name = Path(snakemake.input.final_viral_library).name
_not_assayed = strain_status.query("not assayed").sort_values("strain")
_assayed = strain_status.query("assayed")
_dropped = sorted(_assayed.query("n_barcodes_final == 0")["strain"])
_partial = _assayed.query("final_library == 'kept, some barcodes dropped'").sort_values(
    "strain"
)

add_markdown(f"""
    ## Strains assayed here

    {len(samples)} wells across {len(plates)} plates hold one strain each, all at a dilution
    factor of {_dilutions[0]}. The barcodes are counted against `{_counted_name}`, which holds
    {len(strain_status)} strains.
    """)
add_markdown(f"""
    **{len(_not_assayed)} of those strains are not assayed here**, so there is no data on
    them. Whether each was carried through to `{_final_name}`, the library used for the titer
    measurements, is given alongside:
    """)
add_markdown(
    "\n".join(
        f" - `{row.strain}` \u2014 "
        + ("carried through" if row.n_barcodes_final else "not carried through either")
        for row in _not_assayed.itertuples()
    )
    or " - none"
)
add_markdown(f"""
    **{len(_dropped)} of the {len(_assayed)} strains that are assayed were not carried through
    to `{_final_name}`**. Their bars are outlined and their table rows written in
    <span class="dropped">this color</span> throughout, so that what was measured here can be
    read against the decision to drop them:
    """)
add_markdown("\n".join(f" - `{strain}`" for strain in _dropped) or " - none")
if len(_partial):
    add_markdown(f"""
        A further {len(_partial)} strains were carried through with fewer barcodes than were
        counted here, which matters because a strain's reads are summed over its barcodes.
        These count as carried forward and are not marked:
        """)
    add_markdown(
        "\n".join(
            f" - `{row.strain}` \u2014 {row.n_barcodes} barcodes counted here, "
            f"{row.n_barcodes_final} carried through"
            for row in _partial.itertuples()
        )
    )

# --- tally the reads of each well --------------------------------------------------------


def as_barcode_array(barcodes):
    """Get ``barcodes`` as a 2D array of one nucleotide per cell, for comparing them."""
    lengths = {len(barcode) for barcode in barcodes}
    if len(lengths) != 1:
        raise ValueError(f"the barcodes are not all one length: {sorted(lengths)}")
    return np.frombuffer("".join(barcodes).encode(), dtype="S1").reshape(
        len(barcodes), -1
    )


def min_hamming(barcodes, targets):
    """Get the smallest Hamming distance from each of ``barcodes`` to any of ``targets``."""
    barcodes, targets = list(barcodes), list(targets)
    if not barcodes:
        return np.array([], dtype=int)
    if not targets:
        raise ValueError("no target barcodes to measure a distance to")
    if len(barcodes[0]) != len(targets[0]):
        raise ValueError("the barcodes and the targets are of different lengths")
    return (
        (
            as_barcode_array(barcodes)[:, None, :]
            != as_barcode_array(targets)[None, :, :]
        )
        .sum(axis=2)
        .min(axis=1)
    )


# How far each unmatched barcode is from the barcodes that belong in the well it turned up
# in: those of the strain infected into it, and those of the neutralization standard spiked
# into it. Within `MAX_OWN_STRAIN_ERROR_HAMMING` of either, the read is a misread of
# something that belongs there rather than anything that does not; further away it cannot be
# put down to either, whether it is another strain's material misread or something foreign to
# the library outright.
#
# The distance to the own strain is computed a well at a time and written back by index, so
# that each barcode is measured against the strain of its own well. The standard is the same
# for every well of the group, so that distance is one call over the whole frame.
invalid = invalid.merge(
    samples[["plate", "well", "strain"]].rename(columns={"strain": "own_strain"}),
    on=["plate", "well"],
    validate="many_to_one",
)
invalid["own_strain_distance"] = pd.Series(np.nan, index=invalid.index)
for _own_strain, _group in invalid.groupby("own_strain", sort=False):
    invalid.loc[_group.index, "own_strain_distance"] = min_hamming(
        _group["barcode"], strain_barcodes[_own_strain]
    )
if invalid["own_strain_distance"].isnull().any():
    raise ValueError("an unmatched barcode was not measured against its own strain")
invalid["neut_standard_distance"] = min_hamming(
    invalid["barcode"], neut_standard_barcodes
)
invalid["expected_distance"] = invalid[
    ["own_strain_distance", "neut_standard_distance"]
].min(axis=1)

# `strain_status` covers the whole counted library; only the columns describing a strain are
# wanted here, `assayed` being true of every well by construction.
well_reads = samples.merge(
    strain_status[
        ["strain", "n_barcodes", "n_barcodes_final", "final_library", "carried_forward"]
    ],
    on="strain",
    validate="many_to_one",
)


def reads_where(frame, mask, name):
    """Sum ``frame``'s counts per well where ``mask`` holds, as a column of every well.

    Reindexed onto every well rather than merged, so that a well with no reads of this kind
    gets 0 rather than a null that would change the column's type.

    """
    return (
        frame[mask]
        .groupby(["plate", "well"])["count"]
        .sum()
        .reindex(pd.MultiIndex.from_frame(well_reads[["plate", "well"]]), fill_value=0)
        .rename(name)
        .reset_index(drop=True)
    )


well_reads["neut_standard_reads"] = reads_where(
    counts, counts["neut_standard"], "neut_standard_reads"
)
well_reads["own_strain_reads"] = reads_where(
    counts, counts["strain"] == counts["own_strain"], "own_strain_reads"
)
well_reads["other_strain_reads"] = reads_where(
    counts,
    ~counts["neut_standard"] & (counts["strain"] != counts["own_strain"]),
    "other_strain_reads",
)
_misread = invalid["expected_distance"] <= MAX_OWN_STRAIN_ERROR_HAMMING
well_reads["unmatched_near_expected_reads"] = reads_where(
    invalid, _misread, "unmatched_near_expected_reads"
)
well_reads["unmatched_far_expected_reads"] = reads_where(
    invalid, ~_misread, "unmatched_far_expected_reads"
)
well_reads[NON_NEUT_STANDARD_READS] = well_reads[list(COMPOSITION_CATEGORIES)].sum(
    axis=1
)
well_reads["parsed_reads"] = (
    well_reads[NON_NEUT_STANDARD_READS] + well_reads["neut_standard_reads"]
)

# The categories are built from two of the three files the pipeline writes per well, so
# they have to add back up to what the third one says became of the reads. A category that
# double-counted or dropped a read would otherwise go unnoticed.
for fate, columns in [
    (
        "valid barcode",
        ["neut_standard_reads", "own_strain_reads", "other_strain_reads"],
    ),
    (
        "invalid barcode",
        ["unmatched_near_expected_reads", "unmatched_far_expected_reads"],
    ),
]:
    from_fates = reads_where(fates, fates["fate"] == fate, "from_fates")
    disagree = well_reads.loc[
        well_reads[columns].sum(axis=1) != from_fates, ["plate", "well"]
    ]
    if len(disagree):
        raise ValueError(
            f"the reads counted in {columns} do not match the '{fate}' fate in "
            f"{disagree.to_dict(orient='records')}"
        )

# How much data each well produced, below which what is measured in it should not be trusted.
# Only warns, and drops the well from nothing: every fraction in this report is shown beside
# the counts behind it, so a shallow well is more informative shown with a caveat than hidden.
# Subsumes a well with no reads at all, the degenerate case of the same condition.
well_reads["too_shallow"] = well_reads["parsed_reads"] < MIN_PARSED_READS_PER_WELL
if well_reads["too_shallow"].any():
    add_warning(
        f"these wells have fewer than {MIN_PARSED_READS_PER_WELL} parsed reads, which is "
        "the depth every fraction in this report rests on, so what is reported for them "
        "rests on very little data: "
        f"{well_reads.loc[well_reads['too_shallow'], ['plate', 'well', 'strain', 'parsed_reads']].to_dict(orient='records')}"
    )

# every fraction is of the non-neut-standard parsed reads, so the categories sum to one
for column in COMPOSITION_CATEGORIES:
    well_reads[f"frac_{column}"] = well_reads[column] / well_reads[
        NON_NEUT_STANDARD_READS
    ].where(well_reads[NON_NEUT_STANDARD_READS] > 0)
category_sums = well_reads[[f"frac_{c}" for c in COMPOSITION_CATEGORIES]].sum(axis=1)
if not np.allclose(category_sums[well_reads[NON_NEUT_STANDARD_READS] > 0], 1):
    raise ValueError(
        "the categories of a well do not account for all of its non-neut-standard parsed "
        f"reads: {category_sums[~np.isclose(category_sums, 1)].to_dict()}"
    )

# What does not belong in the well at all, which is what the wells are ordered by: everything
# other than the own strain, less the unmatched reads close enough to be the own strain's own
# barcodes misread. Only an ordering, not a threshold, as a well can reach it through either of
# two unrelated problems, which the two sections below report on the quantity each is about.
_foreign = [
    column
    for column in COMPOSITION_CATEGORIES
    if column not in BELONGS_CATEGORIES and column != EXPECTED_MISREAD_CATEGORY
]
well_reads["frac_foreign"] = well_reads[[f"frac_{c}" for c in _foreign]].sum(axis=1)
well_reads["frac_unmatched"] = well_reads[
    [f"frac_{c}" for c in COMPOSITION_CATEGORIES if c.startswith("unmatched")]
].sum(axis=1)
well_reads["flagged_other_strain"] = (
    well_reads["frac_other_strain_reads"] > single_well_config["max_other_strain_frac"]
)
well_reads["flagged_unmatched"] = (
    well_reads["frac_unmatched"] > single_well_config["max_unmatched_frac"]
)

add_markdown("""
    ## Reads per well and what each well held

    One row per well, across every plate of this group, labeled by the strain the well was
    infected with. Hovering a bar outlines that strain in all three panels, and its tooltip
    names the plate and well the row is from.

    The **left panel** is what became of every read sequenced for the well. Where the reads
    went:
    """)
add_markdown(
    "\n".join(
        f" - **{fate}**: {FATE_EXPLANATIONS[fate]}"
        for fate in FATE_EXPLANATIONS
        if fate in set(fates[fates["count"] > 0]["fate"])
    )
)
add_markdown(f"""
    The **middle panel** is what the well's virus was, as fractions of the reads that could be
    parsed as having a barcode, regardless of whether or not that barcode is supposed to be in
    the viral library, excluding neutralization standard barcodes. The four categories
    partition those reads and so sum to one:

     - **own strain**: a barcode of the strain this well was infected with.
     - **other library strains**: a barcode of a *different* strain of the library, which does
       not belong in this well at all.
     - **within {MAX_OWN_STRAIN_ERROR_HAMMING} nt of strain's or neut standard barcode**:
       matched no known barcode, but is that close to one belonging in this well, so it is a
       misread of something that belongs rather than contamination.
     - **over {MAX_OWN_STRAIN_ERROR_HAMMING} nt from strain's or neut standard barcode**:
       matched no known barcode and is further away than that, so it cannot be put down to
       either. This holds other strains' material misread as well as barcodes foreign to the
       library outright.

    The **right panel** takes the other library strains of the middle panel and says where they
    came from, on the same fraction scale, so its bars are small by construction. Every read of
    another strain is attributed to the well that strain was grown in, and to how far that well
    is from the one the read turned up in; a step onto any of the eight surrounding wells counts
    as one, diagonals included, so wells one apart are the wells that touch. Material one well
    away is liquid carried between neighboring wells. Material from further off, from a well on
    another plate of this group, or from a strain grown nowhere here has to have arrived some
    other way, and distances above {MAX_NAMED_DISTANCE} wells are grouped together since there
    is nothing to tell them apart.

    All three panels are ordered by how much of a well is not the own strain's material at all,
    largest first: the other library strains plus the further unmatched reads. Reads misread off
    a barcode that belongs in the well are deliberately left out of that sum. Read the left
    panel's counts alongside the fractions, as a well whose own strain barely grew shows a large
    fraction of foreign material from a very small amount of it.
    """)
if well_reads["too_shallow"].any():
    add_markdown(f"""
        Wells below {MIN_PARSED_READS_PER_WELL} parsed reads are named in the warnings at the
        top of this report. That is the depth every fraction in this report rests on, so their
        fraction bars should not be read closely; they are drawn along with the rest because
        the counts behind them are visible in the left panel. There is no threshold on counts
        per barcode, a well being expected to hold one strain's barcodes and not the library.
        """)

fates_plotted = fates[fates.groupby("fate")["count"].transform("sum") > 0].merge(
    well_reads[["plate", "well", "strain", "frac_foreign", "carried_forward"]],
    on=["plate", "well"],
    validate="many_to_one",
)
fate_order = sorted(fates_plotted["fate"].unique(), reverse=True)

composition_long = (
    well_reads.melt(
        id_vars=[
            "plate",
            "well",
            "strain",
            "frac_foreign",
            "carried_forward",
            NON_NEUT_STANDARD_READS,
        ],
        value_vars=list(COMPOSITION_CATEGORIES),
        var_name="category",
        value_name="reads",
    )
    .assign(fraction=lambda x: x["reads"] / x[NON_NEUT_STANDARD_READS])
    .replace({"category": COMPOSITION_CATEGORIES})
)
# the bars stack in the order the categories are declared in, which the label alone cannot
# convey to Vega, so the position in that order rides along as a field to sort on
composition_long["category_order"] = composition_long["category"].map(
    {label: i for (i, label) in enumerate(COMPOSITION_CATEGORIES.values())}
)

# --- trace where the other strains in a well came from ------------------------------------

# the wells each strain was grown in, at most one per plate, which is what lets a read of
# another strain be traced back to where that strain actually was
_strain_homes = {}
for _row in samples.itertuples():
    _strain_homes.setdefault(_row.strain, []).append((_row.plate, _row.well))

_other_strain_reads = counts[
    ~counts["neut_standard"]
    & (counts["strain"] != counts["own_strain"])
    & (counts["count"] > 0)
].copy()

_origins = []
for row in _other_strain_reads.itertuples():
    homes = _strain_homes.get(row.strain, [])
    same_plate = [well for (plate, well) in homes if plate == row.plate]
    if same_plate:
        distance = wells_apart(row.well, same_plate[0])
        origin = (
            f"{distance} well{'s' if distance > 1 else ''} away"
            if distance <= MAX_NAMED_DISTANCE
            else f"over {MAX_NAMED_DISTANCE} wells away"
        )
        _origins.append((origin, same_plate[0], distance, f"{same_plate[0]}, {origin}"))
    elif homes:
        _origins.append(
            (
                "another plate here",
                homes[0][1],
                pd.NA,
                f"{homes[0][1]} of {homes[0][0]}",
            )
        )
    else:
        _origins.append(("not assayed here", "", pd.NA, "not assayed here"))

_other_strain_reads[["origin", "source_well", "wells_away", "came_from"]] = (
    pd.DataFrame(_origins, index=_other_strain_reads.index)
)

# The reads of each other strain in each well, with where that strain was grown and that
# strain's share of the well. Everything downstream is built from this rather than from the
# per-barcode rows again: the table below takes the largest few of each well, and the chart's
# third panel sums them within an origin.
#
# Grouped by strain rather than by where it came from, since `wells_away` is null for a strain
# grown on another plate and pandas drops rows with a null grouping key: the origin columns are
# constant within a (well, strain) and so are carried along instead of grouped on.
other_strain_totals = (
    _other_strain_reads.groupby(["plate", "well", "strain"], as_index=False)
    .aggregate(
        reads=pd.NamedAgg("count", "sum"),
        origin=pd.NamedAgg("origin", "first"),
        source_well=pd.NamedAgg("source_well", "first"),
        wells_away=pd.NamedAgg("wells_away", "first"),
        came_from=pd.NamedAgg("came_from", "first"),
    )
    .rename(columns={"strain": "other_strain"})
    .merge(
        well_reads[["plate", "well", NON_NEUT_STANDARD_READS]],
        on=["plate", "well"],
        validate="many_to_one",
    )
    .assign(fraction_of_reads=lambda x: x["reads"] / x[NON_NEUT_STANDARD_READS])
    .sort_values("reads", ascending=False)
)


def top_other_strains(wells, n_top):
    """Get the ``n_top`` other strains contributing most of each of ``wells``' reads."""
    return (
        other_strain_totals.merge(
            wells[["plate", "well"]], on=["plate", "well"], validate="many_to_one"
        )
        .groupby(["plate", "well"])
        .head(n_top)
    )


# ordered nearest first, with the origins that are not a distance last, so the stacking and
# the legend both read from "came from next door" towards "came from somewhere else"
origin_order = [
    origin
    for origin in (
        [
            f"{d} well{'s' if d > 1 else ''} away"
            for d in range(1, MAX_NAMED_DISTANCE + 1)
        ]
        + [
            f"over {MAX_NAMED_DISTANCE} wells away",
            "another plate here",
            "not assayed here",
        ]
    )
    if origin in set(other_strain_totals["origin"])
]

# What the third panel draws: the other strains of a well gathered by where they came from.
# The fraction is taken once from the summed reads rather than by summing the per-strain
# fractions above, which would accumulate a rounding error for no reason.
origins_per_well = (
    other_strain_totals.groupby(["plate", "well", "origin"], as_index=False)
    .aggregate(reads=pd.NamedAgg("reads", "sum"))
    .merge(
        well_reads[
            [
                "plate",
                "well",
                "strain",
                "carried_forward",
                "frac_foreign",
                NON_NEUT_STANDARD_READS,
            ]
        ],
        on=["plate", "well"],
        validate="many_to_one",
    )
    .assign(
        fraction_of_reads=lambda x: x["reads"] / x[NON_NEUT_STANDARD_READS],
        # as with the composition above, the stacking order has to reach Vega as a field
        origin_order=lambda x: x["origin"].map(
            {origin: i for (i, origin) in enumerate(origin_order)}
        ),
    )
)


# --- the reads-per-well chart -------------------------------------------------------------

# Hovering a bar outlines that strain in every panel. The param has to be declared on the
# concatenated chart rather than on a panel to be in scope for all of them, and needs
# `empty=False` or every bar matches when nothing is hovered.
_strain_hover = alt.selection_point(
    name="strain_hover",
    fields=["strain"],
    on="pointerover",
    clear="pointerout",
    empty=False,
)
# Every panel marks a strain the same way, by the color of the stroke around its bars, so the
# encoding is shared and draws one legend for the whole chart rather than one per panel.
_stroke = alt.Stroke(
    "carried_forward",
    title="in the library carried forward",
    scale=alt.Scale(
        domain=list(CARRIED_FORWARD_LABELS.values()),
        range=[RETAINED_COLOR, DROPPED_COLOR],
    ),
    legend=alt.Legend(orient="top", columns=1),
)
# every panel shares this axis, so all of them take the same ordering of strains
_y = alt.Y(
    "strain",
    title=None,
    sort=alt.EncodingSortField("frac_foreign", op="min", order="descending"),
)
# the panels are given palettes that cannot be mistaken for each other, their categories
# being unrelated and their scales resolved independently
_fates_panel = (
    alt.Chart(fates_plotted)
    .encode(
        alt.X("count", title="reads", scale=alt.Scale(nice=False, padding=3)),
        _y.axis(alt.Axis(labelLimit=LABEL_LIMIT)),
        alt.Fill(
            "fate",
            sort=fate_order,
            scale=alt.Scale(scheme="tableau10"),
            legend=alt.Legend(orient="top", columns=1),
        ),
        alt.Order("fate", sort="descending"),
        _stroke,
        strokeWidth=alt.condition(
            _strain_hover, alt.value(STROKE_WIDTH_HOVERED), alt.value(STROKE_WIDTH)
        ),
        tooltip=["plate", "well", "strain", "fate", "count"],
    )
    .mark_bar(height={"band": 0.85})
)
_composition_panel = (
    alt.Chart(composition_long)
    .encode(
        alt.X(
            "fraction",
            title="fraction of the non-neut-standard reads",
            stack=True,
            # to one, so that a bar's length is read against the whole of a well rather than
            # against whatever the widest bar happens to reach
            scale=alt.Scale(domain=[0, 1], nice=False),
            axis=alt.Axis(tickCount=4),
        ),
        # the strain labels come from the panel that shares this axis
        _y.axis(None),
        alt.Fill(
            "category",
            title="barcode category",
            sort=list(COMPOSITION_CATEGORIES.values()),
            scale=alt.Scale(scheme="dark2"),
            legend=alt.Legend(labelLimit=LABEL_LIMIT, orient="top", columns=1),
        ),
        alt.Order("category_order", sort="ascending"),
        _stroke,
        strokeWidth=alt.condition(
            _strain_hover, alt.value(STROKE_WIDTH_HOVERED), alt.value(STROKE_WIDTH)
        ),
        tooltip=[
            "plate",
            "well",
            "strain",
            "category",
            "reads",
            alt.Tooltip("fraction", format=".3g"),
            NON_NEUT_STANDARD_READS,
        ],
    )
    .mark_bar(height={"band": 0.85})
)
_origins_panel = (
    alt.Chart(origins_per_well)
    .encode(
        alt.X(
            "fraction_of_reads",
            title="fraction from another library strain",
            stack=True,
            scale=alt.Scale(nice=False, padding=3),
            axis=alt.Axis(tickCount=3),
        ),
        _y.axis(None),
        alt.Fill(
            "origin",
            title="came from",
            sort=origin_order,
            scale=alt.Scale(scheme="set1"),
            legend=alt.Legend(labelLimit=LABEL_LIMIT, orient="top", columns=1),
        ),
        alt.Order("origin_order", sort="ascending"),
        _stroke,
        strokeWidth=alt.condition(
            _strain_hover, alt.value(STROKE_WIDTH_HOVERED), alt.value(STROKE_WIDTH)
        ),
        tooltip=[
            "plate",
            "well",
            "strain",
            "origin",
            "reads",
            alt.Tooltip("fraction_of_reads", format=".3g"),
        ],
    )
    .mark_bar(height={"band": 0.85})
)
add_chart(
    alt.hconcat(
        *(
            panel.properties(height=alt.Step(10), width=width)
            for (panel, width) in [
                (_fates_panel, 170),
                (_composition_panel, 200),
                (_origins_panel, 130),
            ]
        ),
        spacing=20,
    )
    .add_params(_strain_hover)
    .resolve_scale(y="shared", fill="independent", stroke="shared")
    .configure_axis(grid=False)
    .configure_axisX(labelOverlap="greedy")
)

# --- wells holding another strain -------------------------------------------------------

flagged_other_strain = well_reads.query("flagged_other_strain").sort_values(
    "frac_other_strain_reads", ascending=False
)
add_markdown(f"""
    ## Wells holding another strain

    The wells where more than {single_well_config["max_other_strain_frac"]} of the parsed
    reads are barcodes of another library strain, with the {N_TOP_BARCODES} strains
    contributing most of them. Read the counts alongside the fraction: a well whose own virus
    failed to grow shows a large fraction of another strain from very little of it, and a
    shallow well reaches this threshold on very few reads, which is why `parsed reads` is
    here too.
    """)

if not len(flagged_other_strain):
    add_markdown(
        "No well is over that threshold, so no well holds an appreciable amount of "
        "another strain's material."
    )
else:
    add_table(
        flagged_other_strain[
            [
                "plate",
                "well",
                "strain",
                NON_NEUT_STANDARD_READS,
                "own_strain_reads",
                "other_strain_reads",
                "frac_other_strain_reads",
            ]
        ]
        .merge(
            top_other_strains(flagged_other_strain, N_TOP_BARCODES)
            .assign(
                described=lambda x: x["other_strain"]
                + " ("
                + x["fraction_of_reads"].map("{:.4g}".format)
                + ", "
                + x["came_from"]
                + ")"
            )
            .groupby(["plate", "well"])["described"]
            .apply("; ".join)
            .rename("top other strains")
            .reset_index(),
            on=["plate", "well"],
            how="left",
            validate="one_to_one",
        )
        .rename(
            columns={
                "strain": "own strain",
                NON_NEUT_STANDARD_READS: "non-neut-standard reads",
                "own_strain_reads": "own strain reads",
                "other_strain_reads": "other strain reads",
                "frac_other_strain_reads": "fraction another strain",
            }
        ),
        small_last_column=True,
        dropped=flagged_other_strain["strain"].isin(_dropped).reset_index(drop=True),
    )

# --- barcodes not matching any strain or neut standard ----------------------------------

add_markdown(f"""
    ## Barcodes not matching any strain or neut standard

    Reads that parsed as a barcode but match neither the viral library nor the
    neutralization standard. One within {MAX_SEQUENCING_ERROR_HAMMING} nucleotide of a
    known barcode is sequencing error off it; one further away is a genuinely different
    barcode, and so material that is not in the library at all.

    A well whose virus grew well but whose barcodes are not in the library being counted
    against would show up here and nowhere else, its reads all landing in this category
    while the well looks empty of virus everywhere else in the report.
    """)

flagged_unmatched = well_reads.query("flagged_unmatched").sort_values(
    "frac_unmatched", ascending=False
)
add_markdown(f"""
    The wells where more than {single_well_config["max_unmatched_frac"]} of the parsed reads
    matched neither set, one row each, with the {N_TOP_BARCODES} such barcodes the well holds
    most of. Each is given as its share of the well's non-neut-standard reads, followed by
    what it appears to be: the barcode of something known, misread, or a sequence far enough
    from everything known to be foreign to the library.
    """)

if not len(flagged_unmatched):
    add_markdown("No well is over that threshold.")
else:
    _top_unmatched = (
        invalid.merge(
            flagged_unmatched[["plate", "well"]],
            on=["plate", "well"],
            validate="many_to_one",
        )
        .sort_values("count", ascending=False)
        .groupby(["plate", "well"])
        .head(N_TOP_BARCODES)
        .merge(
            barcode_class.rename(columns={"barcode": "closest_valid_barcode"}),
            on="closest_valid_barcode",
            how="left",
            validate="many_to_one",
        )
        .merge(
            well_reads[["plate", "well", NON_NEUT_STANDARD_READS, "frac_unmatched"]],
            on=["plate", "well"],
            validate="many_to_one",
        )
        .assign(
            is_sequencing_error=lambda x: x[HAMMING_COL]
            <= MAX_SEQUENCING_ERROR_HAMMING,
            fraction_of_reads=lambda x: x["count"] / x[NON_NEUT_STANDARD_READS],
            # What a barcode appears to be. Close enough to a known barcode and it is that
            # barcode misread, so the thing it came off is named; further away the closest
            # match is a coincidence rather than an origin, and naming it would imply an
            # association the distance rules out, so only the distance is given.
            described=lambda x: x["barcode"]
            + " ("
            + x["fraction_of_reads"].map("{:.4g}".format)
            + ", "
            + (
                "sequencing error off "
                + x["strain"].where(~x["neut_standard"], "neut standard")
            ).where(
                x["is_sequencing_error"],
                x[HAMMING_COL].astype(str) + " nt from any known barcode",
            )
            + ")",
        )
        # one row per well, the barcodes of each collapsed into a single cell
        .groupby(
            ["plate", "well", "own_strain", NON_NEUT_STANDARD_READS, "frac_unmatched"],
            as_index=False,
            sort=False,
        )
        .aggregate(described=pd.NamedAgg("described", "; ".join))
        .rename(
            columns={
                "own_strain": "own strain",
                NON_NEUT_STANDARD_READS: "non-neut-standard reads",
                "frac_unmatched": "fraction of well unmatched",
                "described": "top barcodes matching nothing",
            }
        )
    )
    add_table(
        _top_unmatched,
        small_last_column=True,
        dropped=_top_unmatched["own strain"].isin(_dropped).reset_index(drop=True),
    )

# --- estimated titer of each strain ------------------------------------------------------

add_markdown(f"""
    ## Estimated titer of each strain

    An **estimate** of how well each strain grew on its own: its own reads against the
    neutralization standard spiked into its well, `own strain reads / neutralization standard
    reads`. Only those two are in the ratio, so a well holding material that does not belong
    does not get a higher titer for it; that is reported above instead.

    Read it as an estimate rather than a measurement. It is **relative**, in arbitrary units
    rather than infectious ones, saying how much of a strain there was for a fixed amount of
    standard at a dilution factor of {_dilutions[0]}. And the infection may not have been in
    the range over which barcode counts track how much virus a well held, in which case the
    ratio understates or overstates the difference between two strains rather than scaling with
    it. Reads are also summed over all of a strain's barcodes, of which strains have differing
    numbers, so a strain carrying more barcodes gives more reads at the same titer; the barcode
    count is in each point's tooltip.

    The standard's reads are this ratio's denominator, and so are what sets how precisely it
    is measured. A strain whose well holds fewer than
    {MIN_NEUT_STANDARD_READS_FOR_TITER} of them is named in a warning at the top of this
    report, its titer being imprecise however the point sits; this is the one place in the
    report where the amount of standard matters quantitatively. It only warns: the strain keeps
    its point, since a titer that is imprecise is still a measurement.
    """)

# The standard's reads set how precisely the titer is measured, so a well holding too few of
# them has a titer that should not be read closely. A count and not a share: with one strain
# per well the standard's *share* of the reads is `1 / (1 + titer)`, and so says nothing the
# titer does not already say.
well_reads["titer_imprecise"] = (
    well_reads["neut_standard_reads"] < MIN_NEUT_STANDARD_READS_FOR_TITER
)
if well_reads["titer_imprecise"].any():
    add_warning(
        f"these wells have fewer than {MIN_NEUT_STANDARD_READS_FOR_TITER} neutralization "
        "standard reads, which is the denominator of the titer, so their titers cannot be "
        "measured with any precision: "
        f"{well_reads.loc[well_reads['titer_imprecise'], ['plate', 'well', 'strain', 'own_strain_reads', 'neut_standard_reads']].to_dict(orient='records')}"
    )

# A well with no standard reads at all has no titer, there being nothing to divide by. That
# is arithmetic rather than a threshold: every well that has a denominator gets a titer,
# including the ones warned about just above.
well_reads["titer"] = well_reads["own_strain_reads"] / well_reads[
    "neut_standard_reads"
].where(well_reads["neut_standard_reads"] > 0)
no_titer = well_reads[well_reads["titer"].isna()]
if len(no_titer):
    add_warning(
        "these wells have no neutralization standard reads at all, so there is no "
        "denominator to divide by and hence no titer to draw for them: "
        f"{no_titer[['plate', 'well', 'strain']].to_dict(orient='records')}"
    )

titers = well_reads[well_reads["titer"].notna()]
_zero_titers = titers.query("titer == 0")
if len(_zero_titers):
    add_markdown(f"""
        {len(_zero_titers)} {"strain has" if len(_zero_titers) == 1 else "strains have"} a
        titer of exactly zero, having no reads of
        {"its" if len(_zero_titers) == 1 else "their"} own at all. A logarithmic axis cannot
        place zero, so {"it is" if len(_zero_titers) == 1 else "they are"} **absent from the
        right-hand panel** and
        {"appears" if len(_zero_titers) == 1 else "appear"} only on the left:
        {", ".join(f"`{strain}`" for strain in sorted(_zero_titers["strain"]))}.
        """)

if not len(titers):
    add_markdown("No well has a titer to report.")
else:
    # every strain of every plate on one axis, which the shared dilution factor checked
    # earlier is what makes comparable
    _titer_y = alt.Y(
        "strain",
        title=None,
        sort=alt.EncodingSortField("titer", op="min", order="descending"),
    )
    _titer_panels = []
    for _scale_type in ["linear", "log"]:
        _first = not _titer_panels  # only the first panel labels the shared axis
        _titer_panels.append(
            alt.Chart(titers)
            .encode(
                alt.X(
                    "titer",
                    title=f"estimated titer, {_scale_type} scale",
                    scale=alt.Scale(type=_scale_type),
                ),
                (
                    _titer_y.axis(alt.Axis(labelLimit=LABEL_LIMIT))
                    if _first
                    else _titer_y.axis(None)
                ),
                _stroke,
                strokeWidth=alt.condition(
                    _strain_hover,
                    alt.value(POINT_STROKE_WIDTH_HOVERED),
                    alt.value(POINT_STROKE_WIDTH),
                ),
                tooltip=[
                    "plate",
                    "well",
                    "strain",
                    alt.Tooltip("titer", format=".3g"),
                    "own_strain_reads",
                    "neut_standard_reads",
                    "parsed_reads",
                    "n_barcodes",
                    "n_barcodes_final",
                    "final_library",
                    "too_shallow",
                    "titer_imprecise",
                    "flagged_other_strain",
                    "flagged_unmatched",
                ],
            )
            .mark_point(filled=True, fill="gray", size=25)
            .properties(height=alt.Step(10), width=220)
        )
    add_chart(
        alt.hconcat(*_titer_panels, spacing=20)
        .add_params(_strain_hover)
        .resolve_scale(y="shared", stroke="shared")
        .configure_axis(grid=False)
        .configure_axisX(labelOverlap="greedy")
    )

# --- write the outputs -----------------------------------------------------------------

if warnings:
    report[warnings_index] = markdown.markdown(
        "## Warnings\n\n"
        + "\n".join(f" - {warning}" for warning in warnings)
        + "\n\nThese do not stop the analysis, but should be looked into."
    )
    print(f"Finished with {len(warnings)} warnings")

# the single other strain each well holds most of, which is what makes the CSV usable
# without the report; the rest of each well's other strains are in the report's table
_worst_other_strain = top_other_strains(well_reads, 1).rename(
    columns={
        "other_strain": "top_other_strain",
        "reads": "top_other_strain_reads",
        "source_well": "top_other_strain_well",
        "wells_away": "top_other_strain_wells_away",
    }
)[
    [
        "plate",
        "well",
        "top_other_strain",
        "top_other_strain_reads",
        "top_other_strain_well",
        "top_other_strain_wells_away",
    ]
]

well_composition = (
    well_reads.merge(
        _worst_other_strain, on=["plate", "well"], how="left", validate="one_to_one"
    )
    .sort_values(["plate", "position"])[
        [
            "plate",
            "well",
            "strain",
            "dilution_factor",
            "n_barcodes",
            "n_barcodes_final",
            "final_library",
            "parsed_reads",
            "neut_standard_reads",
            NON_NEUT_STANDARD_READS,
            *COMPOSITION_CATEGORIES,
            *(f"frac_{column}" for column in COMPOSITION_CATEGORIES),
            "frac_foreign",
            "frac_unmatched",
            "titer",
            "too_shallow",
            "titer_imprecise",
            "flagged_other_strain",
            "flagged_unmatched",
            "top_other_strain",
            "top_other_strain_reads",
            "top_other_strain_well",
            "top_other_strain_wells_away",
        ]
    ]
    .reset_index(drop=True)
)

well_composition.to_csv(
    snakemake.output.well_composition, index=False, float_format="%.4g"
)
Path(snakemake.output.html).write_text(
    PAGE_TEMPLATE.format(
        dropped_color=DROPPED_COLOR,
        title=f"Single virus per well infections {single_well}",
        body="\n".join(report),
    )
)

print(f"Wrote the report to {snakemake.output.html}")
