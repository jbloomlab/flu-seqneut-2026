"""Check that the converted LaTeX says exactly what the Word manuscript says.

Run by hand after `docx_to_latex.py`, to verify that the conversion moved the prose
across without changing it. It reported `PASS` for that conversion, which is what it was
for. It is **not** a check on the manuscript as it stands: the LaTeX is edited by hand
from here on and the Word document is frozen, so the two are expected to diverge -- and
already do, since the figure and table numbers the document spelled out are now generated
by LaTeX rather than written in the source. Both sides are reduced to the multiset of their words and of
their numbers and compared in both directions, so a dropped, added, or altered word
shows up wherever it happens.

Three things are left out of the comparison because they are deliberately not carried
across verbatim: the reference list, which LaTeX regenerates from `references.bib`; the
citations themselves, which became `\\citep` commands; and Table 1, whose cells the
workflow generates from the pipeline's summary rather than from the document.
"""

import argparse
import collections
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

# Headings after which the document stops being prose that the LaTeX carries verbatim.
STOP_AT = "References"
# Headings that are structure rather than prose: the template prints "Abstract" itself,
# and the Figures and Tables sections dissolve into the floats placed through the text.
EXPECTED_ABSENT = {"Abstract", "Figures", "Tables"}

# A word for this purpose is a run of letters, digits or the marks that sit inside one.
WORD = re.compile(r"[^\W_]+(?:['\u2019\u2010-]?[^\W_]+)*", re.UNICODE)
NUMBER = re.compile(r"\d(?:[\d,.]*\d)?")


def zotero_relations(archive):
    """The relationship ids that point at the Zotero links Google Docs leaves behind."""
    rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
    return {rel.get("Id") for rel in rels if "zotero.org" in (rel.get("Target") or "")}


def docx_text(path):
    """The document's prose, without its citations, reference list or table."""
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
        zotero = zotero_relations(archive)
    out = []
    section = ""
    for para in root.find(W + "body"):
        # the table's cells come from the pipeline's summary, not from the document
        if para.tag != W + "p":
            continue
        style = para.find(W + "pPr/" + W + "pStyle")
        style = style.get(W + "val") if style is not None else ""
        text = "".join(t.text or "" for t in para.iter(W + "t"))
        if style == "Heading1":
            section = text.strip()
        if section == STOP_AT:
            break
        runs = [
            "".join(t.text or "" for t in node.iter(W + "t"))
            for node in para
            # a citation is a \citep in the LaTeX, so it is not compared
            if not (node.tag == W + "hyperlink" and node.get(R + "id") in zotero)
        ]
        out.append("".join(runs))  # joined without a space: a run can split a word
    return "\n".join(out)


def latex_text(paths):
    """The manuscript sources with their markup, citations and table input removed."""
    text = "\n".join(p.read_text(encoding="utf-8") for p in paths)
    text = re.sub(r"(?m)^\s*%.*$", "", text)
    text = re.sub(r"(?m)^\\newcommand\{\\papershort(?:title|authors)\}\{.*$", "", text)
    text = re.sub(r"(?m)^\\textsuperscript\{[^}]*\}", "", text)
    text = re.sub(r"\\textsuperscript\{([^}]*)\}", r"\1", text)
    text = re.sub(r"\\citep\{[^}]*\}", "", text)
    text = re.sub(
        r"\\(?:input|includegraphics|label|graphicspath)\b\s*(\[[^\]]*\])?\{[^}]*\}",
        "",
        text,
    )
    text = re.sub(r"\\href\{[^}]*\}", " ", text)  # keep the visible half of a link
    text = re.sub(r"\\(?:begin|end)\{[^}]*\}(\[[^\]]*\])?", " ", text)
    text = re.sub(r"\\newcommand\{\\[a-zA-Z]+\}", " ", text)
    text = re.sub(r"\\[a-zA-Z]+\*?", " ", text)
    return text.replace("{", " ").replace("}", " ").replace("\\", " ")


def normalize(text):
    """Fold the typographic variants that LaTeX renders identically."""
    for fancy, plain_char in [("\u2019", "'"), ("\u2018", "'"), ("\u2010", "-")]:
        text = text.replace(fancy, plain_char)
    return text


def counts(text):
    text = normalize(text)
    return collections.Counter(WORD.findall(text)), collections.Counter(
        NUMBER.findall(text)
    )


def report(label, missing, extra):
    if not missing and not extra:
        print(f"  {label}: identical")
        return True
    for item, n in sorted(missing.items()):
        print(f"  {label} only in the Word document (x{n}): {item!r}")
    for item, n in sorted(extra.items()):
        print(f"  {label} only in the LaTeX (x{n}): {item!r}")
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docx", default="manuscript.docx", type=pathlib.Path)
    parser.add_argument(
        "--manuscript",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parent.parent.parent / "manuscript",
    )
    args = parser.parse_args()

    doc_words, doc_numbers = counts(docx_text(args.docx))
    tex_words, tex_numbers = counts(
        latex_text([args.manuscript / "frontmatter.tex", args.manuscript / "body.tex"])
    )
    print(
        f"Word document: {sum(doc_words.values())} words, {sum(doc_numbers.values())} numbers"
    )
    print(
        f"LaTeX source:  {sum(tex_words.values())} words, {sum(tex_numbers.values())} numbers"
    )

    missing = doc_words - tex_words
    for word in EXPECTED_ABSENT:
        del missing[word]
    ok = report("word", missing, tex_words - doc_words)
    ok &= report("number", doc_numbers - tex_numbers, tex_numbers - doc_numbers)
    print("PASS: the LaTeX says what the Word document says" if ok else "REVIEW NEEDED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
