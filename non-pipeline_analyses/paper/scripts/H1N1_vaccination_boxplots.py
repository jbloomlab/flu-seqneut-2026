"""Freeze the H1N1 titer and fold-change box plots into one two-panel figure.

Panel A is the box plots of the titers before and after vaccination and panel B those of
the fold changes, each for the two strains the charts start on. Each chart's spec is
pulled from its HTML, its overall title and widgets dropped so the strains it starts on
stay fixed, and the pair rendered the same way the vega-embed menu's "Save as SVG" would
-- but offline and reproducibly -- then stacked with their plot frames aligned. Nothing
derived from the clock is written, so unchanged charts always yield a byte-identical
figure.
"""

import sys
import xml.etree.ElementTree as ET

from interactive_charts import extract_spec, render_svg

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"

# the color key of the titer chart, identified by its title so that it can be moved to a
# row above the plot, which leaves the panels' left margins to the axis titles alone
VACCINATION_LEGEND_TITLE = "vaccination"
VACCINATION_LEGEND = {
    "orient": "top",
    "direction": "horizontal",
    "columns": 2,
    "titleOrient": "left",  # keeps the row one line tall
}

# how many color keys each panel is expected to move; the fold-change chart has none
EXPECTED_KEYS_MOVED = {"A": 1, "B": 0}

# whitespace between the two panels, and the panel letter's size and offset from the
# panel's top left
PANEL_GAP = 20  # px
LETTER_SIZE = 24  # px
LETTER_X = 4  # px
LETTER_Y = 22  # px


def freeze(spec, label):
    """Drop `spec`'s title and widgets, and raise its color key; return the keys raised."""
    del spec["title"]
    # a static figure has no widgets, and unbound params keep the values they start on
    for param in spec["params"]:
        param.pop("bind", None)
    moved = 0
    stack = [spec]
    while stack:
        node = stack.pop()
        if isinstance(node, dict):
            color = node.get("encoding", {}).get("color", {})
            if color.get("title") == VACCINATION_LEGEND_TITLE:
                color["legend"] = VACCINATION_LEGEND
                moved += 1
            stack.extend(node.values())
        elif isinstance(node, list):
            stack.extend(node)
    if moved != EXPECTED_KEYS_MOVED[label]:
        raise ValueError(
            f"panel {label}: moved {moved} color key(s), expected "
            f"{EXPECTED_KEYS_MOVED[label]}"
        )
    return moved


def fmt(x):
    """Format a coordinate, dropping a trailing `.0`."""
    return f"{x:g}"


def namespace_ids(element, prefix):
    """Prefix every id `element` defines, and the references to it, with `prefix`.

    The two panels are drawn by separate renders that number their clip paths alike, so
    in one document the second panel's would resolve to the first's.
    """
    ids = {e.get("id") for e in element.iter() if e.get("id")}
    for e in element.iter():
        for name, value in list(e.attrib.items()):
            if name == "id" and value in ids:
                e.set(name, prefix + value)
            elif value.startswith("#") and value[1:] in ids:
                e.set(name, "#" + prefix + value[1:])
            elif "url(#" in value:
                for old in ids:
                    value = value.replace(f"url(#{old})", f"url(#{prefix}{old})")
                e.set(name, value)


def load_panel(chart_html, label):
    """Render a chart HTML frozen as above: its size, offset, and drawing."""
    spec = extract_spec(chart_html)
    freeze(spec, label)
    root = ET.fromstring(render_svg(spec))
    content = root.find(f"{{{SVG}}}g")
    translate = content.get("transform")
    x, y = (float(v) for v in translate[len("translate(") : -1].split(","))
    namespace_ids(content, label.lower() + "-")
    return {
        "label": label,
        "content": content,
        "offset": (x, y),
        "width": float(root.get("width")),
        "height": float(root.get("height")),
    }


panels = [
    load_panel(snakemake.input.titer_chart_html, "A"),
    load_panel(snakemake.input.fold_change_chart_html, "B"),
]

# any difference in the panels' left offsets is the width of their y-axis titles, which is
# not something the reader should see, so the plot frames are aligned
base_x = max(panel["offset"][0] for panel in panels)
width = 0.0
height = 0.0
placements = []
for panel in panels:
    dx = base_x - panel["offset"][0]
    placements.append((panel, dx, height))
    width = max(width, dx + panel["width"])
    height += panel["height"] + PANEL_GAP
height -= PANEL_GAP

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
for panel, dx, dy in placements:
    placed = ET.SubElement(
        figure, f"{{{SVG}}}g", {"transform": f"translate({fmt(dx)},{fmt(dy)})"}
    )
    placed.append(panel["content"])
    letter = ET.SubElement(
        figure,
        f"{{{SVG}}}text",
        {
            "x": fmt(LETTER_X),
            "y": fmt(dy + LETTER_Y),
            "font-family": "sans-serif",
            "font-size": fmt(LETTER_SIZE),
            "font-weight": "bold",
        },
    )
    letter.text = panel["label"]
    print(
        f"placed panel {panel['label']} ({fmt(panel['width'])} x {fmt(panel['height'])}) "
        f"at ({fmt(dx)}, {fmt(dy)})"
    )

ET.register_namespace("", SVG)
ET.register_namespace("xlink", XLINK)
with open(snakemake.output.figure_svg, "w", encoding="utf-8") as f:
    f.write(ET.tostring(figure, encoding="unicode"))
    f.write("\n")

print(f"wrote {fmt(width)} x {fmt(height)} figure to {snakemake.output.figure_svg}")
