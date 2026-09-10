"""Freeze the stratified H1N1 titer chart into a static figure.

The chart splits sera by the ratio of their titers against two strains chosen from
dropdowns, at a threshold set by a slider, and reports the split in lines of text above
the plot. This keeps the strains, threshold, and split the chart opens with, folds what
those lines of text say into the color key, and lays the key and the age density beside
each other above the titer plot.

The key's labels have to name the numbers the chart computes -- how many sera fell on
each side of the threshold -- so they are read back out of the text the chart itself
draws rather than recomputed here, which would duplicate the pipeline's arithmetic.
"""

import json
import re
import sys

from interactive_charts import extract_spec, mark_type, render_svg

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# The chart's own readouts, which this figure reads its labels out of. The strains and
# the threshold are whichever the chart opens with, and the counts follow from them.
READOUT_SPLIT = re.compile(
    r"splitting sera at a (?P<threshold>[\d.]+)-fold titer ratio of "
    # lazy, so the clause below is what ends the reference rather than being read as
    # part of its name; a `derived_haplotype` can itself hold commas
    r"(?P<comparator>.+) to (?P<reference>.+?)"
    # the chart appends this when it drops the sera censored against either strain. The
    # key names the ratio, and its counts already come from the chart's own totals, so
    # the clause has nothing to add to the figure.
    r"(?:, dropping sera at the lower limit of detection against either strain)?"
)
READOUT_COUNTS = re.compile(
    r"n=(?P<above>\d+) above threshold, n=(?P<at_or_below>\d+) at or below, "
    r"n=(?P<unmeasured>\d+) not measured against both strains"
)

# What the color key becomes. The chart colors by the comparator's titer divided by the
# reference's, and puts a serum above the threshold when that ratio is strictly greater
# than it, so the labels carry the inequalities that way round.
KEY_TITLE = "ratio of {comparator} to {reference} titer"
KEY_LABELS = {
    "above threshold": "> {threshold} (n = {above})",
    "at or below threshold": "≤ {threshold} (n = {at_or_below})",
}

# the key, moved out of the gutter Vega gives a legend and into the row above the titer
# plot, where `LEGEND_ROW_LAYOUT` left-aligns it and the age density sits to its right
KEY_LEGEND = {"orient": "left", "columns": 1}

# identifies the color encoding this figure relabels, and the categories it must cover
KEY_ENCODING_TITLE = "comparator vs reference titer"
KEY_CATEGORIES = ["above threshold", "at or below threshold"]

# gap between the age density and the titer plot below it, which the chart sets to 1 for
# the rows of readout text this figure drops
ROW_SPACING = 12


def chart_parts(spec, html_path):
    """The `(text rows, titer plot and tree, age density)` of the stratified chart.

    The chart is rows of readout text above an `hconcat` of the titer plot over the
    tree, with the age density to their right. Identify each by shape and error on
    anything else, since that means the chart has changed and this figure, which
    rearranges those parts, has to be revisited.
    """
    text_rows = [row for row in spec["vconcat"] if mark_type(row) == "text"]
    plot_rows = [row for row in spec["vconcat"] if "hconcat" in row]
    if len(text_rows) + len(plot_rows) != len(spec["vconcat"]) or len(plot_rows) != 1:
        raise ValueError(
            f"{html_path}: expected rows of text above one row of plots, found "
            f"{len(text_rows)} text and {len(plot_rows)} plot rows"
        )
    if len(plot_rows[0]["hconcat"]) != 2:
        raise ValueError(f"{html_path}: plot row is not a pair of columns")
    panels, density = plot_rows[0]["hconcat"]
    if len(panels.get("vconcat", [])) != 2:
        raise ValueError(f"{html_path}: first column is not the titer plot and tree")
    if mark_type(density) != "area":
        raise ValueError(f"{html_path}: second column is not the age density")
    return (text_rows, panels, density)


def read_readouts(spec, text_rows, html_path):
    """What the chart's lines of readout text say, as one dict of named values.

    Rendering them is how their values are read: they are the output of the chart's own
    transforms over its sliders and dropdowns, so nothing short of running that
    dataflow knows what they say.
    """
    readouts = json.loads(json.dumps(spec))
    readouts["vconcat"] = text_rows
    readouts.pop("title", None)
    drawn = re.findall(r"<text[^>]*>([^<]*)</text>", render_svg(readouts))
    print(f"chart's readout text: {drawn}")

    values = {}
    for pattern in (READOUT_SPLIT, READOUT_COUNTS):
        matches = [m for m in map(pattern.fullmatch, drawn) if m]
        if len(matches) != 1:
            raise ValueError(
                f"{html_path}: {len(matches)} of its readouts match "
                f"{pattern.pattern!r}, expected 1"
            )
        values.update(matches[0].groupdict())
    return values


def relabel_key(panels, density, readouts, html_path):
    """Move the color key to the row above the plot and name the split it shows."""
    encodings = [
        layer["encoding"]["color"]
        for layer in panels["vconcat"][0]["layer"]
        if layer.get("encoding", {}).get("color", {}).get("title") == KEY_ENCODING_TITLE
    ]
    if not encodings:
        raise ValueError(f"{html_path}: found no {KEY_ENCODING_TITLE!r} encoding")
    for encoding in encodings:
        if encoding["scale"]["domain"] != KEY_CATEGORIES:
            raise ValueError(
                f"{html_path}: splits sera into "
                f"{encoding['scale']['domain']}, not the {KEY_CATEGORIES} "
                f"this figure has a label for"
            )
        # the key moves to the age density, which draws the same categories, so that the
        # two sit side by side in one row rather than in a row each
        encoding["legend"] = None

    if density["encoding"]["color"]["legend"] is not None:
        raise ValueError(f"{html_path}: the age density already draws a key")
    label_expr = " : ".join(
        f"datum.label === {json.dumps(category)} ? "
        f"{json.dumps(label.format(**readouts))}"
        for category, label in KEY_LABELS.items()
    )
    density["encoding"]["color"]["legend"] = {
        **KEY_LEGEND,
        "title": KEY_TITLE.format(**readouts),
        # every category is named above, so nothing should reach this fallback
        "labelExpr": f"{label_expr} : datum.label",
    }
    print(f"key title: {KEY_TITLE.format(**readouts)!r}")
    for category, label in KEY_LABELS.items():
        print(f"  {category!r} -> {label.format(**readouts)!r}")


spec = extract_spec(snakemake.input.chart_html)
text_rows, panels, density = chart_parts(spec, snakemake.input.chart_html)
readouts = read_readouts(spec, text_rows, snakemake.input.chart_html)
relabel_key(panels, density, readouts, snakemake.input.chart_html)

# The rows of readout text and the title go; the age density moves out of the column
# beside the plots and into a row above them, where it sits to the right of the key.
spec["vconcat"] = [density, panels]
spec["spacing"] = ROW_SPACING
del spec["title"]

# a static figure has no widgets
for param in spec["params"]:
    param.pop("bind", None)

svg = render_svg(spec)
with open(snakemake.output.figure_svg, "w", encoding="utf-8") as f:
    f.write(svg)

print(f"dropped {len(text_rows)} rows of readout text and the title")
print(
    "sera measured against only one of the two strains, which the chart excludes and "
    f"the key does not count: {readouts['unmeasured']}"
)
print("rendered " + re.search(r'width="\d+" height="\d+"', svg).group())
print(f"wrote {snakemake.output.figure_svg}")
