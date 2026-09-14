# The Word manuscript the LaTeX was converted from

`manuscript.docx` is the manuscript as it stood in Google Docs on September 13, 2026,
exported from
https://docs.google.com/document/d/114EQK8n18Iq5AU1QisaOMNsDnG3JNO41A1PvDwvefDE. It is
kept as the record of where [../../manuscript/](../../manuscript) came from.

**This is history, not a step in the workflow.** The manuscript source is the LaTeX, and
that is where edits go; nothing here is re-run, and re-running it would overwrite the
LaTeX with the state of a document that is no longer current.

## What the conversion did

`docx_to_latex.py` wrote `frontmatter.tex`, `body.tex` and `references.bib`. `pandoc`
converted the prose and everything the script does is a mechanical transform of that
output, so the wording in the LaTeX is the wording in the document. On top of that it:

- split the title, authors, affiliations and abstract into the title block;
- moved each figure and table caption into a float placed at its first mention in the
  text, dropping the Word document's Tables and Figures sections;
- replaced the Zotero citation links with `\citep`, resolving each in-text author-year
  to exactly one reference and failing if any were ambiguous or uncited;
- built the bibliography by looking up each reference's DOI on Crossref, which is also
  how the author lists the document had shortened to "et al." were recovered. Where
  Crossref records a different year than the manuscript cites, the manuscript wins, so
  that no in-text citation changed.

Table 1 is not carried across: the workflow generates it from the pipeline's own summary
of the sera, and the manuscript inputs that.

`fidelity_check.py` then verified the result, comparing the multiset of words and of
numbers on each side in both directions. It reports `PASS` for this conversion.

## Running them

Both need the `seqneut-pipeline` conda environment; the conversion additionally needs
`pandoc` and network access for Crossref. On the Hutch cluster:

    module load Pandoc/2.13
    python docx_to_latex.py
    python fidelity_check.py
