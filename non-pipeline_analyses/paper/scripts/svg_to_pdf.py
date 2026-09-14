"""Render a figure's SVG to the PDF that the manuscript includes.

LaTeX cannot read SVG, and converting the tracked SVG rather than re-rendering the chart
keeps the figure in the typeset manuscript identical to the figure the analysis drew.

`vl_convert` writes a PDF's resource dictionary in an order that varies between runs, so
converting the same SVG twice gives two different files even though they draw the same
thing. `qpdf` writes dictionary keys in sorted order, so rewriting the PDF through it
makes the figure -- and the manuscript that embeds it -- byte-reproducible.
"""

import subprocess
import sys
import tempfile

import vl_convert

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input.figure_svg, encoding="utf-8") as f:
    svg = f.read()
print(f"read {len(svg)} characters of SVG from {snakemake.input.figure_svg}")

with tempfile.NamedTemporaryFile(suffix=".pdf") as rendered:
    rendered.write(vl_convert.svg_to_pdf(svg))
    rendered.flush()
    subprocess.run(
        [
            "qpdf",
            "--deterministic-id",  # no random file id, which also varies per run
            "--object-streams=disable",
            rendered.name,
            snakemake.output.figure_pdf,
        ],
        check=True,
    )
print(f"wrote {snakemake.output.figure_pdf}")
