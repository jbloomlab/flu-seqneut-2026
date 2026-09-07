"""Render a hand-written Markdown report to a self-contained HTML page.

What a report may write -- the `docs:`, `repo:` and `embed:` link prefixes, and the
`{repo_url}` it can interpolate -- is described for whoever is writing one in
`data/reports/README.md`. Each link is rewritten to a real URL and checked here, except
that whether a `docs:` target is really a page of the site is left to
`scripts/check_report_links.py`, which runs once `build_docs` has built that site.

"""

import base64
import csv
import html
import os
import pathlib
import re
import string
import subprocess
import sys

import markdown
import markdown.extensions
import markdown.extensions.toc
import markdown.preprocessors
import markdown.treeprocessors

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

#: Heading levels shown in the sidebar. The `#` heading is the report's title, so the
#: sidebar starts at `##`; `####` and deeper render but would make the sidebar a wall.
TOC_DEPTH = "2-3"

#: Height in px given to an embed that cannot be measured in the browser and that names
#: no height of its own. Only reached by a cross-origin embed, such as a tree.
DEFAULT_EMBED_HEIGHT = 800

#: MIME type of each image format a `figure:` may name
FIGURE_MIME_TYPES = {
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
}

#: The option each kind takes after its target, as `{<option>=<number>}`
KIND_OPTIONS = {"embed": "height", "figure": None, "table": "tfoot"}

#: `![caption](embed:target)`, `![caption](figure:target)` or `![caption](table:target)`
#: alone on a line, optionally followed by the kind's option. Anywhere else these are an
#: error, since none of an iframe, a captioned figure or a table can sit inside the flow
#: of a paragraph.
INLINE_RE = re.compile(
    r"^!\[(?P<caption>[^\]]*)\]\((?P<kind>embed|figure|table):(?P<target>[^)\s]+)\)"
    r"(?:\{(?P<option>\w+)=(?P<value>\d+)\})?[ \t]*$",
    re.MULTILINE,
)


def tracked_paths():
    """Every path tracked by git, as a set of files plus the directories holding them."""
    files = subprocess.run(
        ["git", "ls-files"], capture_output=True, check=True, text=True
    ).stdout.splitlines()
    paths = set(files)
    for f in files:
        parent = pathlib.PurePosixPath(f).parent
        while str(parent) != ".":
            paths.add(str(parent))
            parent = parent.parent
    return paths


def resolve(href):
    """Turn one authored link target into the URL it stands for, or raise."""
    if href.startswith("#"):
        return href
    scheme, _, target = href.partition(":")
    if scheme == "docs":
        # `build_docs` copies every page into `results/docs` by basename, and this report
        # is copied in beside them, so the basename is the whole link
        page, sep, fragment = target.partition("#")
        return os.path.basename(page) + sep + fragment
    if scheme == "repo":
        if target not in tracked:
            raise ValueError(
                f"`repo:{target}` is not tracked by git, so it would 404 on GitHub. "
                "Note that most of `results/` is deliberately untracked."
            )
        # `HEAD` rather than a named branch, so the link follows the default branch
        return f"{repo_url}/blob/HEAD/{target}"
    if scheme in {"http", "https", "mailto"}:
        return href
    raise ValueError(
        f"link target {href!r} names no known scheme; write `docs:`, `repo:`, `embed:`, "
        "`https://`, `mailto:`, or `#anchor`"
    )


def embed_html(caption, target, height):
    """The raw HTML of one inline plot or tree: an iframe plus its caption."""
    src = resolve(target if ":" in target else f"docs:{target}")
    # A cross-origin frame cannot be measured from this page, so an absolute URL keeps
    # whatever height it is given; a page of our own docs site is measured and resized.
    attrs = ""
    if height or "//" in src:
        attrs = f' data-height="{height or DEFAULT_EMBED_HEIGHT}"'
    quoted = html.escape(src, quote=True)
    return (
        '<figure class="embed">\n'
        f'<iframe class="embed" src="{quoted}" loading="lazy" '
        f'title="{html.escape(caption, quote=True)}"{attrs}></iframe>\n'
        f"<figcaption>{html.escape(caption)} "
        f'(<a href="{quoted}" target="_blank" rel="noopener">open in a new tab</a>)'
        "</figcaption>\n"
        "</figure>"
    )


def figure_html(caption, target):
    """The raw HTML of one static figure: an image plus its caption."""
    if target.startswith(("http://", "https://")):
        src = target  # an image on another site, which the browser fetches
    else:
        suffix = pathlib.Path(target).suffix.lower()
        if suffix not in FIGURE_MIME_TYPES:
            raise ValueError(
                f"`figure:{target}` has extension {suffix!r}, but a figure must be one "
                f"of {sorted(FIGURE_MIME_TYPES)}"
            )
        if target not in figures:
            raise ValueError(
                f"`figure:{target}` is not among the figures `rules/reports.smk` found "
                "in this report, so a change to it would not rebuild the report"
            )
        # Inlined rather than linked, so the page carries its own figures. `build_docs`
        # copies only the pages it is given, so a linked file would not be beside it.
        data = base64.b64encode(pathlib.Path(target).read_bytes()).decode()
        src = f"data:{FIGURE_MIME_TYPES[suffix]};base64,{data}"
    return (
        '<figure class="figure">\n'
        f'<img src="{html.escape(src, quote=True)}" '
        f'alt="{html.escape(caption, quote=True)}">\n'
        f"<figcaption>{html.escape(caption)}</figcaption>\n"
        "</figure>"
    )


def table_html(caption, target, tfoot):
    """The raw HTML of one table read from a CSV, plus its caption.

    The last `tfoot` rows are put in the table's footer, which is how a row totalling
    the rows above it is marked; the CSV itself just holds that row last.

    """
    if pathlib.Path(target).suffix.lower() != ".csv":
        raise ValueError(f"`table:{target}` must name a `.csv` written by a rule")
    if target not in tables:
        raise ValueError(
            f"`table:{target}` is not among the tables `rules/reports.smk` found in "
            "this report, so a change to it would not rebuild the report"
        )
    with open(target, newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2 + tfoot:
        raise ValueError(f"`table:{target}` has too few rows: {len(rows)}")

    def cells(row, tag):
        return (
            "<tr>" + "".join(f"<{tag}>{html.escape(c)}</{tag}>" for c in row) + "</tr>"
        )

    body = rows[1 : len(rows) - tfoot]
    parts = [
        "<thead>" + cells(rows[0], "th") + "</thead>",
        "<tbody>" + "".join(cells(row, "td") for row in body) + "</tbody>",
    ]
    if tfoot:
        parts.append(
            "<tfoot>" + "".join(cells(row, "td") for row in rows[-tfoot:]) + "</tfoot>"
        )
    return (
        '<figure class="table">\n'
        f"<table>{''.join(parts)}</table>\n"
        f"<figcaption>{html.escape(caption)}</figcaption>\n"
        "</figure>"
    )


class EmbedPreprocessor(markdown.preprocessors.Preprocessor):
    """Replace each embed, figure or table line with the raw HTML it stands for."""

    def run(self, lines):
        def stash(match):
            caption = " ".join(match.group("caption").split())
            kind, target = match.group("kind"), match.group("target")
            option, value = match.group("option"), match.group("value")
            if option and option != KIND_OPTIONS[kind]:
                takes = KIND_OPTIONS[kind]
                raise ValueError(
                    f"`{kind}:{target}` names `{option}=`, but a {kind} takes "
                    + (f"only `{takes}=`" if takes else "no option")
                )
            if kind == "embed":
                html_text = embed_html(caption, target, value)
            elif kind == "figure":
                html_text = figure_html(caption, target)
            else:
                html_text = table_html(caption, target, int(value or 0))
            return self.md.htmlStash.store(html_text)

        return INLINE_RE.sub(stash, "\n".join(lines)).split("\n")


class LinkTreeprocessor(markdown.treeprocessors.Treeprocessor):
    """Rewrite and check every link once the inline processor has built them."""

    def run(self, root):
        for element in root.iter("a"):
            element.set("href", resolve(element.get("href")))
        for element in root.iter("img"):
            raise ValueError(
                f"{element.get('src')!r} is an image whose target is not an `embed:`, "
                "a `figure:` or a `table:`, or is one that is not alone on its line"
            )


class ReportExtension(markdown.extensions.Extension):
    """The report's link rewriting, as a Python-Markdown extension."""

    def extendMarkdown(self, md):
        # below `fenced_code` (25), so an embed shown inside a code fence stays literal
        md.preprocessors.register(EmbedPreprocessor(md), "report_embed", 24)
        md.treeprocessors.register(LinkTreeprocessor(md), "report_links", 5)


report = snakemake.wildcards.report
repo_url = snakemake.params.repo_url.rstrip("/")
figures = set(snakemake.input.figures)
tables = set(snakemake.input.tables)
tracked = tracked_paths()

text = pathlib.Path(snakemake.input.markdown).read_text()
# Replacing just this token, rather than formatting the whole text, so that the
# braces a report writes for other reasons -- an embed's `{height=}` -- are untouched
text = text.replace("{repo_url}", repo_url)

headings = [line for line in text.split("\n") if re.fullmatch(r"#[^#].*", line)]
if len(headings) != 1:
    raise ValueError(
        f"{snakemake.input.markdown} must have exactly one `#` heading, its title, but "
        f"has {len(headings)}: {headings}"
    )
if not text.lstrip().startswith("# "):
    raise ValueError(f"{snakemake.input.markdown} must open with its `#` title heading")
title = headings[0].removeprefix("# ").strip()
print(f"Rendering report {report!r} titled {title!r}")

md = markdown.Markdown(
    extensions=[
        "fenced_code",
        "tables",
        "sane_lists",
        markdown.extensions.toc.TocExtension(permalink=True, toc_depth=TOC_DEPTH),
        ReportExtension(),
    ]
)
body = md.convert(text)

page = string.Template(pathlib.Path(snakemake.input.template).read_text()).substitute(
    title=html.escape(title), toc=md.toc, body=body
)
pathlib.Path(snakemake.output.html).write_text(page)
print(f"Wrote {snakemake.output.html}")
