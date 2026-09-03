"""Check that every local link a rendered report makes resolves in the built docs site.

`build_docs` copies each page into `results/docs` by basename, and the report is copied
in beside them, so a link the report makes to another page is just that page's basename.
Checking those against the directory `build_docs` actually built -- rather than against a
list of the pages it was expected to build -- is what keeps this from drifting as the
pipeline gains page types, and covers every page of the site.

"""

import pathlib
import re
import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

#: `href` and `src` attribute values in the rendered HTML
LINK_RE = re.compile(r"\b(?:href|src)=\"([^\"]+)\"")

docs = pathlib.Path(snakemake.input.docs)
report = pathlib.Path(snakemake.input.html)

targets = []
for link in dict.fromkeys(LINK_RE.findall(report.read_text())):
    # anything with a scheme, a protocol-relative URL, or an anchor within this page is
    # not the docs site's to resolve
    if link.startswith(("#", "//")) or re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", link):
        continue
    targets.append(link)

lines = [f"Local links made by {report}, checked against {docs}", ""]
missing = []
for target in sorted(targets):
    page = target.partition("#")[0]
    if (docs / page).is_file():
        lines.append(f"  ok       {target}")
    else:
        lines.append(f"  MISSING  {target}")
        missing.append(target)

lines += ["", f"{len(targets) - len(missing)} of {len(targets)} local links resolve."]
text = "\n".join(lines) + "\n"
print(text)

# written only once every link resolves, so a failed run leaves no report saying it did
if missing:
    raise ValueError(
        f"{report} links to {len(missing)} target(s) that are not pages of {docs}: "
        + ", ".join(missing)
    )
pathlib.Path(snakemake.output.check).write_text(text)
