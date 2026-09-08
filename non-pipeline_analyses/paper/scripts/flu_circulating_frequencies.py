"""Combine the H3N2 and H1N1 circulating-sequence panels into one paper figure.

The panels are composed as SVG rather than re-plotted from the underlying counts, so
this figure always shows exactly what the `flu_circulating_frequencies` analysis
produced. Only the standard library is used, and nothing derived from the clock is
written, so an unchanged pair of panels always yields a byte-identical figure.
"""

import collections
import re
import sys
import xml.etree.ElementTree as ET

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

SVG = "http://www.w3.org/2000/svg"

# classes Vega gives the parts of a chart that this script has to find
LEGEND_CLASS = "mark-group role-legend"
FRAME_CLASS = "mark-group role-frame root"

# Layout of this figure: H3N2 above H1N1, stacked so each panel is as wide as the page
# allows. Both panels carry an identical color key, so the lower panel's is stripped and
# that panel cropped at the key's left edge, and the upper panel's key is moved down the
# empty right-hand column to sit vertically centered beside both panels.

# whitespace between the upper panel's x-axis title and the lower panel's title
PANEL_GAP = 15

Panel = collections.namedtuple(
    "Panel", ["path", "content", "offset", "width", "height", "plot_area"]
)


def parse_translate(transform, context):
    """(x, y) of a `translate(x,y)` transform."""
    match = re.fullmatch(r"translate\(([-\d.]+),([-\d.]+)\)", transform or "")
    if not match:
        raise ValueError(f"{context}: expected translate transform, got {transform!r}")
    return (float(match.group(1)), float(match.group(2)))


def parse_box(d, context):
    """(width, height) of a rectangular `M{x},{y}h{w}v{h}h-{w}Z` path."""
    match = re.fullmatch(r"M[\d.]+,[\d.]+h([\d.]+)v([\d.]+)h-\1Z", d or "")
    if not match:
        raise ValueError(f"{context}: expected a rectangular path, got {d!r}")
    return (float(match.group(1)), float(match.group(2)))


def fmt(x):
    """Format a coordinate, dropping a trailing `.0`."""
    return f"{x:g}"


def find_one(element, path, context, what):
    """The single match of `path` in `element`, erroring if there is not exactly one."""
    matches = element.findall(path)
    if len(matches) != 1:
        raise ValueError(f"{context}: expected 1 {what}, found {len(matches)}")
    return matches[0]


def load_panel(path):
    """Read a Vega-rendered chart SVG, checking it has the structure composed below."""
    root = ET.parse(path).getroot()
    # the chart's drawing, offset from the SVG origin to leave room for axis labels; the
    # panel's own white background rect is a sibling and so is dropped by taking this
    content = find_one(root, f"{{{SVG}}}g", path, "top-level <g>")
    frame = find_one(
        content, f"{{{SVG}}}g[@class='{FRAME_CLASS}']/{{{SVG}}}g", path, "chart frame"
    )
    plot_area = find_one(
        frame, f"{{{SVG}}}path[@class='background']", path, "plot area"
    ).get("d")
    return Panel(
        path=path,
        content=content,
        offset=parse_translate(content.get("transform"), path),
        width=float(root.get("width")),
        height=float(root.get("height")),
        plot_area=plot_area,
    )


def legend_parts(panel):
    """A panel's color key: its group, the translated group inside it, and its x."""
    legends = [e for e in panel.content.iter() if e.get("class") == LEGEND_CLASS]
    if len(legends) != 1:
        raise ValueError(f"{panel.path}: expected 1 color key, found {len(legends)}")
    inner = find_one(legends[0], f"{{{SVG}}}g", panel.path, "color key group")
    x, _ = parse_translate(inner.get("transform"), panel.path)
    return (legends[0], inner, x)


def strip_legend(panel):
    """Remove a panel's color key, returning the panel width left of where it was."""
    legend, _, legend_x = legend_parts(panel)
    parents = {child: parent for parent in panel.content.iter() for child in parent}
    parents[legend].remove(legend)
    return panel.offset[0] + legend_x


def center_legend(panel, panel_dy, figure_height):
    """Move a panel's color key to be vertically centered in the composed figure."""
    _, inner, legend_x = legend_parts(panel)
    box = find_one(
        inner, f"{{{SVG}}}path[@class='background']", panel.path, "color key box"
    )
    _, legend_height = parse_box(box.get("d"), panel.path)
    y = (figure_height - legend_height) / 2 - panel_dy - panel.offset[1]
    inner.set("transform", f"translate({fmt(legend_x)},{fmt(y)})")
    return y


top = load_panel(snakemake.input.h3n2_svg)
bottom = load_panel(snakemake.input.h1n1_svg)

if top.plot_area != bottom.plot_area:
    raise ValueError(
        "panels draw different-sized plot areas so cannot be aligned:\n"
        f"  {top.path}: {top.plot_area}\n"
        f"  {bottom.path}: {bottom.plot_area}"
    )

bottom_width = strip_legend(bottom)

# Align the plot frames: any difference in the panels' left offsets is the width of
# their y-axis labels, which is not something the reader should see.
base_x = max(top.offset[0], bottom.offset[0])
placements = [
    (top, base_x - top.offset[0], 0.0, top.width),
    (bottom, base_x - bottom.offset[0], top.height + PANEL_GAP, bottom_width),
]

width = max(dx + panel_width for _, dx, _, panel_width in placements)
height = max(dy + panel.height for panel, _, dy, _ in placements)
legend_y = center_legend(top, placements[0][2], height)

figure = ET.Element(
    f"{{{SVG}}}svg",
    {
        "version": "1.1",
        "class": "marks",
        "width": fmt(width),
        "height": fmt(height),
        "viewBox": f"0 0 {fmt(width)} {fmt(height)}",
    },
)
ET.SubElement(
    figure,
    f"{{{SVG}}}rect",
    {"width": fmt(width), "height": fmt(height), "fill": "white"},
)
for panel, dx, dy, panel_width in placements:
    placed = ET.SubElement(
        figure, f"{{{SVG}}}g", {"transform": f"translate({fmt(dx)},{fmt(dy)})"}
    )
    placed.append(panel.content)
    print(
        f"placed {panel.path} ({fmt(panel_width)} x {fmt(panel.height)}) "
        f"at ({fmt(dx)}, {fmt(dy)})"
    )

print(f"moved the color key to y={fmt(legend_y)} within {top.path}")

ET.register_namespace("", SVG)
with open(snakemake.output.figure_svg, "w", encoding="utf-8") as f:
    f.write(ET.tostring(figure, encoding="unicode"))
    f.write("\n")

print(f"wrote {fmt(width)} x {fmt(height)} figure to {snakemake.output.figure_svg}")
