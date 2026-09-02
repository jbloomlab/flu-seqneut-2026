"""Protein-structure visualizations built by `prot-struct-viz`.

Each visualization is a directory under `data/prot-struct-viz_config/` holding a
`spec.yaml` and the CSVs and Markdown captions it names. Paths inside a spec resolve
relative to the spec file, so such a directory can be moved or copied intact.

"""

import pathlib

# may be absent or null: a project need not have any structure visualizations
prot_struct_viz = config.get("prot-struct-viz") or {}

# `build_docs` copies every HTML into `results/docs` by basename and rejects duplicates,
# so the directory name alone could collide with another analysis' chart
prot_struct_viz_html = "results/prot-struct-viz/{viz}_prot_struct_viz.html"


rule prot_struct_viz_page:
    """Render one self-contained HTML page showing data on a protein structure."""
    input:
        spec="data/prot-struct-viz_config/{viz}/spec.yaml",
        # The spec's other inputs sit beside it and are named only inside it, which cannot
        # be read here as `prot-struct-viz` is not in the environment running snakemake.
        # Everything else in the directory is therefore taken as an input, less the
        # scratch files that a leading `_` marks.
        spec_dir_files=lambda wc: sorted(
            str(f)
            for f in pathlib.Path("data/prot-struct-viz_config", wc.viz).iterdir()
            if f.is_file() and f.name != "spec.yaml" and not f.name.startswith("_")
        ),
    output:
        html=prot_struct_viz_html,
        report=prot_struct_viz_html.removesuffix(".html") + "_report.txt",
    log:
        "results/logs/prot_struct_viz_{viz}.txt",
    wildcard_constraints:
        viz="|".join(prot_struct_viz),
    conda:
        "../envs/prot-struct-viz.yml"
    shell:
        # `--out` rather than an `out` key in the spec, so no path out of `data/` and into
        # `results/` is written into a data file
        "prot-struct-viz {input.spec} --out {output.html} &> {log}"


# guarded because an empty section would render in the docs as a heading with no links
if prot_struct_viz:
    add_htmls_to_docs["Protein structure visualizations"] = {
        _title: prot_struct_viz_html.format(viz=_viz)
        for (_viz, _title) in prot_struct_viz.items()
    }

# the HTMLs are not listed as they are already inputs to `build_docs`
prot_struct_viz_outputs = expand(
    prot_struct_viz_html.removesuffix(".html") + "_report.txt",
    viz=prot_struct_viz,
)
