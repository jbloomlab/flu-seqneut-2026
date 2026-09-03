"""Narrative reports hand-written in Markdown and rendered to HTML pages in the docs.

Each report is one Markdown file in `data/reports/`, named for the key that configures
it. Its links to the rest of the project are rewritten and checked when it is rendered,
and the links it makes to other pages of the docs site are checked against that site
once `build_docs` has built it.

"""

import pathlib
import re

# may be absent or null: a project need not have any reports
reports = config.get("reports") or {}

#: A report's `figure:` targets. `scripts/render_report.py` parses the whole figure
#: syntax; only the target is needed here, to declare the file as an input so that
#: editing a figure rebuilds the report that inlines it.
FIGURE_RE = re.compile(r"\]\(figure:([^)\s]+)\)")

# `build_docs` copies every HTML into `results/docs` by basename and rejects duplicates,
# so the report name alone could collide with another analysis' chart
report_html = "results/reports/{report}.html"


rule render_report:
    """Render one hand-written Markdown report to a self-contained HTML page."""
    input:
        markdown=lambda wc: reports[wc.report]["markdown"],
        template="scripts/report_template.html",
        figures=lambda wc: sorted(
            {
                target
                for target in FIGURE_RE.findall(
                    pathlib.Path(reports[wc.report]["markdown"]).read_text()
                )
                if not target.startswith(("http://", "https://"))
            }
        ),
    output:
        html=report_html,
    log:
        "results/logs/render_report_{report}.txt",
    wildcard_constraints:
        report="|".join(reports),
    conda:
        "../seqneut-pipeline/environment.yml"
    params:
        repo_url=lambda wc: reports[wc.report]["repo_url"],
    script:
        "../scripts/render_report.py"


rule check_report_links:
    """Check a rendered report's local links against the docs site `build_docs` built."""
    input:
        docs=rules.build_docs.output.docs,
        html=report_html,
    output:
        check="results/reports/{report}_link_check.txt",
    log:
        "results/logs/check_report_links_{report}.txt",
    wildcard_constraints:
        report="|".join(reports),
    conda:
        "../seqneut-pipeline/environment.yml"
    script:
        "../scripts/check_report_links.py"


# guarded because an empty section would render in the docs as a heading with no links
if reports:
    add_htmls_to_docs["Reports"] = {
        _settings["title"]: report_html.format(report=_report)
        for (_report, _settings) in reports.items()
    }

# the HTMLs are not listed as they are already inputs to `build_docs`
reports_outputs = expand("results/reports/{report}_link_check.txt", report=reports)
