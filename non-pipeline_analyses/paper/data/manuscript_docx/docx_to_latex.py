"""Convert the Word manuscript in this directory to the LaTeX sources in `manuscript/`.

This is a one-time conversion, run by hand, whose output is committed as the manuscript
source; it is not part of the workflow. See README.md in this directory for how to run
it and what it needs.

The prose is never retyped. `pandoc` does the text conversion and everything here is a
mechanical transform of its output, so the words reaching `body.tex` are the words in the
`.docx`. Anything that cannot be resolved unambiguously is a hard error, not a guess.
"""

import html
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

HERE = pathlib.Path(__file__).resolve().parent
DOCX = HERE / "manuscript.docx"
OUT_DIR = HERE.parent.parent / "manuscript"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

# Sections of the Word document that do not survive as body text: the abstract moves to
# the front matter, Tables and Figures dissolve into floats placed near each first
# mention, and References becomes the bibliography.
DROPPED_SECTIONS = {"Abstract", "Tables", "Figures", "References"}

N_FIGURES = 8

# The running head needs short forms, which the Word document does not carry.
SHORT_AUTHORS = "Kikawa, Butler et al."
SHORT_TITLE = "Human neutralizing antibody landscape to influenza, summer 2026"

# Each figure is set across both columns. The height cap keeps a tall figure and its
# caption on one page; the figures here are up to 6.5 by 8.1 inches.
FIGURE_GRAPHIC = (
    r"\includegraphics[width=\textwidth,height=0.72\textheight,keepaspectratio]{%s}"
)

# The table the pipeline generates, relative to `manuscript/`.
TABLE_INPUT = "../tables/Table_1.tex"

# Crossref turns a DOI into BibTeX. Their etiquette policy asks that a request identify
# what is making it; this project's repository does that.
CROSSREF = "https://api.crossref.org/works/{doi}/transform/application/x-bibtex"
CROSSREF_AGENT = "flu-seqneut-2026 (https://github.com/jbloomlab/flu-seqneut-2026)"
CROSSREF_DELAY = 0.15  # seconds between requests, to stay inside their rate limit

# The only references Crossref does not index are WHO vaccine-composition
# recommendations, which are reports rather than articles.
WHO_ENTRY = """@techreport{{{key},
    author = {{{{World Health Organization}}}},
    title = {{{title}}},
    institution = {{World Health Organization}},
    year = {{{year}}},
    url = {{{url}}}
}}"""


# ---------------------------------------------------------------------------
# reading the Word document
# ---------------------------------------------------------------------------


def docx_paragraphs(path):
    """Every paragraph of the document, as (style, text)."""
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    out = []
    for p in root.find(W + "body").iter(W + "p"):
        style = p.find(W + "pPr/" + W + "pStyle")
        out.append(
            (
                style.get(W + "val") if style is not None else "",
                "".join(t.text or "" for t in p.iter(W + "t")),
            )
        )
    return out


def run_pandoc(workdir):
    """The document as a LaTeX fragment, one line per paragraph."""
    if shutil.which("pandoc") is None:
        sys.exit("pandoc is not on PATH; see README.md in this directory")
    out = workdir / "body_raw.tex"
    subprocess.run(
        ["pandoc", str(DOCX), "-t", "latex", "--wrap=none", "-o", str(out)], check=True
    )
    return out.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# structure of the pandoc output
# ---------------------------------------------------------------------------


def flatten_headings(text):
    r"""Collapse pandoc's `\hypertarget{..}{% \section{..}\label{..}}` to the heading.

    pandoc wraps every heading in an anchor whose closing brace lands on the heading
    line. Removing the wrapper line-wise is brace-safe, which a regex spanning the
    heading is not: these titles contain `\texorpdfstring{..}{..}`.
    """
    lines = text.split("\n")
    out = []
    i = 0
    while i < len(lines):
        if not re.fullmatch(r"\\hypertarget\{[^}]*\}\{%", lines[i]):
            out.append(lines[i])
            i += 1
            continue
        # gather until the brace the anchor opened closes; a heading the document left
        # blank spans two lines, as pandoc renders it `\texorpdfstring{\hfill\break`
        heading, depth = [], 1
        i += 1
        while i < len(lines) and depth:
            depth += lines[i].count("{") - lines[i].count("}")
            heading.append(lines[i])
            i += 1
        if depth:
            sys.exit(f"unterminated heading anchor at {heading[0]!r}")
        # kept to one line so the heading survives splitting the text by lines
        out.append(" ".join(x.strip() for x in heading)[:-1])
    return "\n".join(out)


def brace_group(s, start):
    """The contents of the brace group beginning at `s[start]`, and the index after it."""
    depth = 0
    for i in range(start, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return s[start + 1 : i], i + 1
    sys.exit(f"unbalanced braces in {s[start : start + 80]!r}")


def heading_of(line):
    r"""(depth, title) if the line is a `\section`/`\subsection`, else None."""
    m = re.match(r"\\(?:sub)*section(?:\[[^\]]*\])?\{", line)
    if not m:
        return None
    title, _ = brace_group(line, m.end() - 1)
    # pandoc wraps a title holding markup for the PDF bookmark; keep the visible half
    inner = re.match(r"\\texorpdfstring\{", title)
    if inner:
        title, _ = brace_group(title, inner.end() - 1)
    return line.count("sub", 0, m.end()), title.strip()


def drop_graphics(tex):
    """Remove embedded images.

    The Word document carries the figures as pasted images, and one of them sits inside
    a caption's heading. The manuscript takes its figures from the pipeline instead, so
    every one of these is dropped.
    """
    return re.sub(
        r"\\(?:protect\\?)?\\?includegraphics(?:\[[^\]]*\])?\{[^{}]*\}", "", tex
    )


def plain(tex):
    """Text with inline LaTeX markup removed, for matching and for headings."""
    out = re.sub(
        r"\\(?:emph|textbf|textit|underline|ul|hl|protect)\{([^{}]*)\}",
        r"\1",
        drop_graphics(tex),
    )
    out = re.sub(r"\\[a-zA-Z]+(?:\[[^\]]*\])?\s*", "", out)
    return out.replace("{", "").replace("}", "").strip()


def split_sections(text):
    r"""(front, sections) where sections is [(title, [lines])] split at each `\section`."""
    front, sections = [], []
    for line in text.split("\n"):
        head = heading_of(line)
        if head and head[0] == 0:
            sections.append((plain(head[1]), []))
        elif sections:
            sections[-1][1].append(line)
        else:
            front.append(line)
    return front, sections


def section(sections, name):
    """The one section with this title."""
    hits = [body for title, body in sections if title == name]
    if len(hits) != 1:
        sys.exit(f"expected exactly one {name!r} section, found {len(hits)}")
    return hits[0]


# ---------------------------------------------------------------------------
# cleaning up Google Docs artifacts
# ---------------------------------------------------------------------------

ZOTERO = (
    r"\\href\{https://www\.zotero\.org/google-docs/\?[^}]*\}\{((?:[^{}]|\{[^{}]*\})*)\}"
)


def strip_zotero(text):
    """Unwrap the Zotero links Google Docs leaves on citations and reference entries."""
    return re.sub(ZOTERO, r"\1", text)


def unlink(text):
    """Drop pandoc's underlining of hyperlinks, a Google Docs styling artifact."""
    while True:
        new = re.sub(r"\\underline\{((?:[^{}]|\{[^{}]*\})*)\}", r"\1", text)
        if new == text:
            return new
        text = new


# ---------------------------------------------------------------------------
# references
# ---------------------------------------------------------------------------


def parse_references(lines):
    """Each reference entry, as the fields needed to key and look it up."""
    refs = []
    for line in lines:
        entry = unlink(strip_zotero(line)).strip()
        if not entry:
            continue
        m = re.match(r"^(.*?)\s((?:19|20)\d\d[a-z]?)\.\s", entry)
        if not m:
            sys.exit(f"cannot find the author/year of reference: {entry[:100]!r}")
        authors, year = m.group(1), m.group(2)
        doi = re.search(r"https?://doi\.org/(\S+?)\.?\s*$", entry)
        url = re.search(r"(https?://\S+?)\.?\s*$", entry)
        refs.append(
            {
                "entry": entry,
                "year": year,
                "doi": doi.group(1) if doi else None,
                "url": url.group(1) if url else None,
                "title": re.search(r"``(.*?)''", entry).group(1),
                **name_parts(authors),
            }
        )
    return refs


def name_parts(authors):
    """The surnames a citation is disambiguated by, from a reference's author list.

    Entries read `Surname, Given, Given Surname, and Given Surname`, optionally cut
    short with `et al.`; an organization has no comma at all.
    """
    text = re.sub(r",?\s*et al\.?$", "", authors.rstrip(".").strip())
    if "," not in text:
        return {"surname": text, "given": "", "second": ""}
    parts = [p for p in re.split(r",\s*and\s+|\s+and\s+|,\s*", text) if p]
    rest = parts[2:]
    return {
        "surname": parts[0],
        "given": parts[1] if len(parts) > 1 else "",
        "second": rest[0].split()[-1] if rest else "",
    }


def assign_keys(refs):
    """A unique BibTeX key per reference, from its first author and year."""
    acronym = lambda s: "".join(w[0] for w in s.split() if w[:1].isupper()).lower()
    groups = {}
    for ref in refs:
        surname = ref["surname"]
        stem = acronym(surname) if " " in surname else surname.lower()
        stem = re.sub(r"[^a-z]", "", stem) + ref["year"]
        groups.setdefault(stem, []).append(ref)
    for stem, group in groups.items():
        for i, ref in enumerate(group):
            ref["key"] = stem if len(group) == 1 else f"{stem}{chr(ord('a') + i)}"


def clean_bibtex(bibtex):
    """Turn the HTML that Crossref puts in titles and journal names into LaTeX.

    Crossref serves fields as HTML, so an ampersand arrives as an entity and markup as
    tags. Left alone, a bare `&` is an alignment tab and stops the build.
    """
    bibtex = re.sub(r"<sup>(.*?)</sup>", r"\\textsuperscript{\1}", bibtex)
    bibtex = re.sub(r"<sub>(.*?)</sub>", r"\\textsubscript{\1}", bibtex)
    bibtex = re.sub(r"</?[a-z]+>", "", bibtex)
    bibtex = html.unescape(bibtex)
    return re.sub(r"(?<!\\)([&%#])", r"\\\1", bibtex)


def crossref_bibtex(doi):
    """The BibTeX record Crossref holds for a DOI."""
    request = urllib.request.Request(
        CROSSREF.format(doi=urllib.parse.quote(doi)),
        headers={"User-Agent": CROSSREF_AGENT},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8").strip()


def build_bibliography(refs):
    """`references.bib`, from Crossref for everything it indexes."""
    entries = []
    for ref in refs:
        if ref["doi"]:
            print(f"  fetching {ref['key']}: {ref['doi']}")
            bibtex = clean_bibtex(crossref_bibtex(ref["doi"]))
            time.sleep(CROSSREF_DELAY)
            bibtex = re.sub(r"^(\s*@\w+\{)[^,]*,", rf"\1{ref['key']},", bibtex, count=1)
            # The citation must read the year the manuscript cites. Where Crossref
            # records a different one -- an online-first date against the issue date
            # Zotero used -- the manuscript wins, or the in-text citation would change.
            found = re.search(r"year=\{?(\d{4})", bibtex)
            if found and found.group(1) != ref["year"]:
                print(f"    year {found.group(1)} -> {ref['year']} (as cited)")
                bibtex = re.sub(
                    r"year=\{?\d{4}\}?", f"year={{{ref['year']}}}", bibtex, count=1
                )
        else:
            bibtex = WHO_ENTRY.format(
                key=ref["key"],
                title=ref["title"],
                year=ref["year"],
                url=ref["url"],
            )
        entries.append(bibtex)
    return "\n\n".join(entries) + "\n"


# ---------------------------------------------------------------------------
# citations
# ---------------------------------------------------------------------------


def initials(given):
    return "".join(w[0] for w in given.replace(".", " ").split() if w[:1].isalpha())


def resolve_citation(name, year, refs):
    """The single reference an in-text `Author et al. Year` names."""
    candidates = [r for r in refs if r["year"] == year]
    name = name.strip()
    if name.endswith("et al."):
        name = name[: -len("et al.")].strip().rstrip(",")
    tokens = [t.strip() for t in name.split(",") if t.strip()]
    head, extra = tokens[0], tokens[1:]
    if " and " in head:
        first, second = (x.strip() for x in head.split(" and "))
        candidates = [
            r for r in candidates if r["surname"] == first and r["second"] == second
        ]
    else:
        # Zotero prefixes initials where two first authors share a surname
        given = re.match(r"^((?:[A-Z]\.\s*)+)(\S.*)$", head)
        if given:
            head = given.group(2).strip()
            wanted = given.group(1).replace(".", "").replace(" ", "")
        candidates = [r for r in candidates if r["surname"] == head]
        if given:
            candidates = [
                r for r in candidates if initials(r["given"]).startswith(wanted)
            ]
    # and adds further surnames where the first author and year still do not separate
    for token in extra:
        exact = [r for r in candidates if r["second"] == token]
        candidates = exact or [
            r
            for r in candidates
            if re.search(r"\b" + re.escape(token) + r"\b", r["entry"])
        ]
    if len(candidates) != 1:
        sys.exit(
            f"{name!r} {year} matches {len(candidates)} references"
            + "".join(f"\n    {c['entry'][:90]}" for c in candidates)
        )
    return candidates[0]


def replace_citations(text, refs):
    r"""Rewrite each Zotero citation group as a `\citep`, and record what was cited."""
    cited = set()

    def rewrite(match):
        group = match.group(1).strip()
        if not group.startswith("("):
            return match.group(0)  # a reference-list entry, not a citation
        keys = []
        for part in group[1:-1].split(";"):
            part = part.strip()
            split = re.match(
                r"^(.*?)\s((?:19|20)\d\d[a-z]?(?:,\s*(?:19|20)\d\d[a-z]?)*)$", part
            )
            if not split:
                sys.exit(f"cannot read the citation {part!r}")
            # one author cited for several years reads `Author 2026, 2025`
            for year in re.findall(r"(?:19|20)\d\d[a-z]?", split.group(2)):
                ref = resolve_citation(split.group(1), year, refs)
                keys.append(ref["key"])
                cited.add(ref["key"])
        return r"\citep{" + ",".join(keys) + "}"

    return re.sub(ZOTERO, rewrite, text), cited


# ---------------------------------------------------------------------------
# floats
# ---------------------------------------------------------------------------


def parse_captions(lines, pattern):
    r"""{number: caption} for headings matching `pattern`, plus the text following each.

    A caption in the Word document is a bold heading paragraph -- the sentence naming
    the figure or table -- followed by the rest of the caption as ordinary text.
    """
    captions = {}
    current = None
    for line in lines:
        head = heading_of(line)
        if head:
            found = re.match(pattern, plain(head[1]))
            if found:
                title, _ = brace_group(line, line.index("{", line.index("section")))
                inner = re.match(r"\\texorpdfstring\{", title)
                if inner:
                    title, _ = brace_group(title, inner.end() - 1)
                current = found.group(1)
                captions[current] = [drop_graphics(title).strip()]
            else:
                current = None  # a heading holding only the image
        elif current and line.strip():
            if line.lstrip().startswith(r"\begin{"):
                current = None  # the figure or table itself, not more of its caption
            elif "includegraphics" not in line:
                captions[current].append(line.strip())
    # the naming sentence is the heading, which the document sets bold, and the rest is
    # the caption proper
    return {n: (parts[0], " ".join(parts[1:])) for n, parts in captions.items()}


def caption_of(parts):
    """The caption, with its naming sentence bold as the Word document sets it."""
    head, rest = parts
    return rf"\textbf{{{head}}} {rest}".strip()


def figure_float(number, parts):
    """A figure set across both columns, captioned as the Word document captions it."""
    return "\n".join(
        [
            r"\begin{figure*}[tp]",
            r"\centering",
            FIGURE_GRAPHIC % f"Figure_{number}",
            rf"\caption*{{{caption_of(parts)}}}",
            rf"\label{{fig:{number}}}",
            r"\end{figure*}",
        ]
    )


def table_float(number, parts):
    """The table the pipeline generates, captioned as the Word document captions it."""
    return "\n".join(
        [
            r"\begin{table*}[tp]",
            r"\centering",
            rf"\caption*{{{caption_of(parts)}}}",
            rf"\label{{tab:{number}}}",
            rf"\input{{{TABLE_INPUT}}}",
            r"\end{table*}",
        ]
    )


def insert_floats(lines, floats):
    """Place each float just after the paragraph that first mentions it.

    Every mention is located in the prose as it stands before any float is added, so a
    paragraph naming two figures keeps them in numerical order rather than letting the
    second one land ahead of the first.
    """
    placed = {}
    for mention, block in floats:
        where = next(
            (
                i
                for i, line in enumerate(lines)
                if re.search(rf"\b{mention}\b", line) and not heading_of(line)
            ),
            None,
        )
        if where is None:
            sys.exit(f"no paragraph in the body mentions {mention!r}")
        placed.setdefault(where, []).append(block)
    out = []
    for i, line in enumerate(lines):
        out.append(line)
        for block in placed.get(i, []):
            out += ["", block]
    return out


# ---------------------------------------------------------------------------
# assembling the output
# ---------------------------------------------------------------------------


def build_frontmatter(title, front, abstract):
    """`frontmatter.tex`: the title block the template's preamble reads."""
    text = "\n".join(front)
    authors = unlink(text.split(r"\begin{enumerate}")[0]).strip()
    affiliations = [
        m.strip()
        for m in re.findall(
            r"\\item\s*\n\s*\\begin\{quote\}\s*\n(.*?)\n\s*\\end\{quote\}",
            text,
            re.DOTALL,
        )
    ]
    if not affiliations:
        sys.exit("could not find the numbered affiliations")
    notes = [
        line.strip()
        for line in text.rsplit(r"\begin{quote}", 1)[1]
        .split(r"\end{quote}")[0]
        .split("\n")
        if line.strip()
    ]
    block = [rf"\textsuperscript{{{i}}}{a}\\" for i, a in enumerate(affiliations, 1)]
    block += [rf"{unlink(n)}\\" for n in notes]
    block[-1] = block[-1][:-2]  # no line break after the last line
    return "\n".join(
        [
            "% Generated from manuscript.docx; see data/manuscript_docx/README.md.",
            rf"\newcommand{{\papertitle}}{{{title}}}",
            rf"\newcommand{{\paperauthors}}{{{authors}}}",
            r"\newcommand{\paperaffiliations}{%",
            *block,
            "}",
            r"\newcommand{\paperabstract}{%",
            abstract,
            "}",
            rf"\newcommand{{\papershortauthors}}{{{SHORT_AUTHORS}}}",
            rf"\newcommand{{\papershorttitle}}{{{SHORT_TITLE}}}",
            "",
        ]
    )


def build_body(sections, figures, tables, refs):
    """`body.tex`: the manuscript prose, with the floats placed."""
    lines = []
    for title, body in sections:
        if title in DROPPED_SECTIONS or not title:
            continue
        lines.append(rf"\section{{{title}}}")
        for line in body:
            head = heading_of(line)
            if head and not head[1]:
                continue  # an empty heading, from a blank heading paragraph in the doc
            if head:
                lines.append(rf"\subsection{{{head[1]}}}")
            else:
                lines.append(line)
    floats = [(f"Figure {n}", figure_float(n, c)) for n, c in sorted(figures.items())]
    floats += [(f"Table {n}", table_float(n, c)) for n, c in sorted(tables.items())]
    lines = insert_floats(lines, floats)
    text, cited = replace_citations("\n".join(lines), refs)
    uncited = {r["key"] for r in refs} - cited
    if uncited:
        sys.exit(f"references never cited in the text: {sorted(uncited)}")
    text = unlink(strip_zotero(text))
    return re.sub(r"\n{3,}", "\n\n", text).strip() + "\n"


def main():
    with tempfile.TemporaryDirectory() as tmp:
        raw = run_pandoc(pathlib.Path(tmp))
    text = flatten_headings(raw)
    front, sections = split_sections(text)

    titles = [t.strip() for style, t in docx_paragraphs(DOCX) if style == "Title"]
    if len(titles) != 1:
        sys.exit(f"expected one Title paragraph, found {len(titles)}")

    refs = parse_references(section(sections, "References"))
    assign_keys(refs)
    print(f"read {len(refs)} references")

    figures = parse_captions(section(sections, "Figures"), r"Figure (\d+)\.")
    if len(figures) != N_FIGURES:
        sys.exit(f"found {len(figures)} figure captions, expected {N_FIGURES}")
    tables = parse_captions(section(sections, "Tables"), r"Table (\d+)\.")
    print(f"read {len(figures)} figure and {len(tables)} table captions")

    abstract = unlink(
        strip_zotero("\n".join(x for x in section(sections, "Abstract") if x.strip()))
    )

    OUT_DIR.mkdir(exist_ok=True)
    (OUT_DIR / "frontmatter.tex").write_text(
        build_frontmatter(titles[0], front, abstract), encoding="utf-8"
    )
    (OUT_DIR / "body.tex").write_text(
        build_body(sections, figures, tables, refs), encoding="utf-8"
    )
    print("fetching the bibliography from Crossref")
    (OUT_DIR / "references.bib").write_text(build_bibliography(refs), encoding="utf-8")
    print(f"wrote frontmatter.tex, body.tex and references.bib to {OUT_DIR}")


if __name__ == "__main__":
    main()
