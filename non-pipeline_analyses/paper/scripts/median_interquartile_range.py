"""Freeze one of the interactive titer charts into a static figure.

The pipeline's charts are Vega-Lite `vconcat`s whose sera are filtered by bound sliders,
and, for the charts drawn over whole cohorts, whose cohorts are chosen by clicking a
legend. This pulls the spec out of the chart's HTML, fixes the cohort selection where
there is one, drops the parts that do not belong in print, and renders it the same way
the vega-embed menu's "Save as SVG" would -- but offline and reproducibly. Nothing
derived from the clock is written, so an unchanged chart always yields a byte-identical
figure.

The charts differ only in what they draw, so one script freezes each of them; the
`figure` wildcard says which.
"""

import json
import re
import sys

import vl_convert

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# Cohorts to show, one faceted row each, in the charts that offer a cohort selection:
# every cohort those charts hold except the pre- and post-vaccination ones, which are
# what the vaccination figures are about. `All` is the pooled row over all sera. The
# vaccination charts offer no cohort selection and show every panel they draw.
COHORTS = ["All", "CTS", "SCH", "UWMC", "VIDRL"]

# every cohort left out of the figure must end in one of these; anything else means the
# chart has gained a cohort someone has to decide about rather than silently drop
EXCLUDED_COHORT_SUFFIXES = ("_pre", "_post")

# the chart's clickable cohort legend, identified by its title so it can be dropped
COHORT_LEGEND_TITLE = "serum cohort (click to select)"

# The pre-/post-vaccination charts' color key, identified the same way. Vega gives a
# faceted chart's legend a gutter to the left of the whole plot, which for a two-entry
# key is mostly empty, so put it in a tight row above the plot instead. Vega-Lite cannot
# align a legend within that row, so `LEGEND_ROW_LAYOUT` aligns it at render time.
VACCINATION_LEGEND_TITLE = "vaccination"
VACCINATION_LEGEND = {
    "orient": "top",
    "direction": "horizontal",
    "columns": 2,
    "titleOrient": "left",  # keeps the row one line tall
}
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


def legend_title(node):
    """The title of a spec node's fill legend, or None if it has no fill encoding."""
    return node.get("encoding", {}).get("fill", {}).get("title")


def is_panel_pair(node):
    """Whether a spec node is the figure's panels: the faceted titer plot and tree."""
    children = node.get("vconcat", [])
    return len(children) == 2 and "facet" in children[0] and "layer" in children[1]


def figure_panels(spec, html_path):
    """The panels of `spec`, checking that whatever else it holds is droppable.

    Each chart wraps the panels this figure wants alongside, in some of the charts, a
    text readout of the slider settings and an invisible chart carrying the clickable
    cohort legend. Find the panels by shape and require everything else to be one of
    those two, so a chart that has gained a panel fails here rather than silently
    losing it.
    """
    panels = [child for child in spec["vconcat"] if is_panel_pair(child)]
    if len(panels) != 1:
        raise ValueError(f"{html_path}: found {len(panels)} panel pairs, expected 1")

    extras = [child for child in spec["vconcat"] if not is_panel_pair(child)]
    for extra in extras:
        # the readout of what the sliders are set to, and the cohort legend
        if mark_type(extra) != "text" and legend_title(extra) != COHORT_LEGEND_TITLE:
            raise ValueError(
                f"{html_path}: do not know whether to drop {json.dumps(extra)[:200]}"
            )
    return (panels[0], extras)


def raise_vaccination_legend(panels):
    """Move the pre-/post-vaccination color key to a row above the plot.

    Returns how many encodings were moved, which is none for the charts drawn over
    whole cohorts: their only color key is the cohort legend, which is dropped.
    """
    moved = 0
    for chart in panels["vconcat"]:
        for layer in chart.get("spec", chart).get("layer", []):
            color = layer.get("encoding", {}).get("color", {})
            if color.get("title") == VACCINATION_LEGEND_TITLE:
                color["legend"] = VACCINATION_LEGEND
                moved += 1
    return moved


def cohort_selections(spec):
    """Every param of `spec` that selects cohorts."""
    return [
        param
        for param in spec["params"]
        if param.get("select", {}).get("fields") == ["cohort"]
    ]


def cohort_legend(extras, html_path):
    """The dropped chart carrying the clickable cohort legend, if the chart has one."""
    legends = [extra for extra in extras if legend_title(extra) == COHORT_LEGEND_TITLE]
    if len(legends) > 1:
        raise ValueError(f"{html_path}: found {len(legends)} cohort legends")
    return legends[0] if legends else None


def check_cohorts(legend, html_path):
    """Check `COHORTS` against every cohort the chart offers."""
    in_chart = legend["encoding"]["fill"]["scale"]["domain"]
    missing = [cohort for cohort in COHORTS if cohort not in in_chart]
    if missing:
        raise ValueError(f"{html_path}: no such cohort {missing}, chart has {in_chart}")
    unexpected = [
        cohort
        for cohort in in_chart
        if cohort not in COHORTS and not cohort.endswith(EXCLUDED_COHORT_SUFFIXES)
    ]
    if unexpected:
        raise ValueError(
            f"{html_path}: cohort {unexpected} is neither shown in this figure nor a "
            f"pre-/post-vaccination cohort; decide whether the figure should show it"
        )


spec = extract_spec(snakemake.input.chart_html)
panels, extras = figure_panels(spec, snakemake.input.chart_html)
legend = cohort_legend(extras, snakemake.input.chart_html)
selections = cohort_selections(spec)

# A chart either offers a cohort selection, with the legend to click, or draws a fixed
# set of panels; requiring the legend and the param to agree is what keeps a chart that
# lost its legend from quietly having every cohort selected instead.
if legend is None:
    if selections:
        raise ValueError(
            f"{snakemake.input.chart_html}: has a cohort selection but no legend"
        )
    print("chart offers no cohort selection, so every panel it draws is shown")
else:
    if len(selections) != 1:
        raise ValueError(
            f"{snakemake.input.chart_html}: found {len(selections)} cohort "
            f"selections alongside its legend, expected 1"
        )
    check_cohorts(legend, snakemake.input.chart_html)
    selections[0]["value"] = [{"cohort": cohort} for cohort in COHORTS]
    print(f"cohorts shown: {', '.join(COHORTS)}")

# Keep the panels nested in their own `vconcat` rather than splicing them into the root,
# so the chart's `resolve`, `spacing`, and `center` settings still apply as they do in
# the browser. Dropping the title takes the text above the plot with it.
spec["vconcat"] = [panels]
del spec["title"]

# a static figure has no widgets; this also stops the cohort selection's
# `bind: "legend"` from referring to the legend chart dropped just above
for param in spec["params"]:
    param.pop("bind", None)

moved = raise_vaccination_legend(panels)
print(f"moved {moved} color key(s) to a row above the plot")

# Compile to Vega rather than rendering the Vega-Lite spec directly, so the legend row
# above the plot can be right-aligned; a chart with no legend in that row renders
# byte-identically either way.
vega = vl_convert.vegalite_to_vega(json.dumps(spec))
vega["config"]["legend"]["layout"] = LEGEND_ROW_LAYOUT
svg = vl_convert.vega_to_svg(json.dumps(vega))
with open(snakemake.output.figure_svg, "w", encoding="utf-8") as f:
    f.write(svg)

print("rendered " + re.search(r'width="\d+" height="\d+"', svg).group())
print(f"wrote {snakemake.output.figure_svg}")
