"""Helpers for freezing the pipeline's interactive Altair charts into static figures.

Imported by the figure scripts that do so, each of which decides for itself which parts
of its chart to keep; nothing here writes a file or looks at the clock, so an unchanged
chart always yields a byte-identical figure.
"""

import json
import re

import vl_convert

# Vega-Lite cannot align a legend within the row above a plot, so rendering goes through
# Vega, where this layout left-aligns one. A chart with no legend in that row renders
# byte-identically with or without it.
LEGEND_ROW_LAYOUT = {"top": {"anchor": "start", "direction": "horizontal"}}


def extract_spec(html_path):
    """The Vega-Lite spec embedded in a chart HTML written by Altair."""
    with open(html_path) as f:
        html = f.read()
    matches = re.findall(r"var spec = (\{.*?\});\n\s*var embedOpt", html, re.DOTALL)
    if len(matches) != 1:
        raise ValueError(
            f"{html_path}: found {len(matches)} embedded specs, expected 1"
        )
    return json.loads(matches[0])


def mark_type(node):
    """A spec node's mark type, or None if it is not a single-mark chart."""
    mark = node.get("mark")
    if isinstance(mark, dict):
        return mark.get("type")
    return mark


def render_svg(spec):
    """`spec` rendered to an SVG string."""
    vega = vl_convert.vegalite_to_vega(json.dumps(spec))
    vega["config"]["legend"]["layout"] = LEGEND_ROW_LAYOUT
    return vl_convert.vega_to_svg(json.dumps(vega))
