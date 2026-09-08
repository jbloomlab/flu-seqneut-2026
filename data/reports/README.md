# Writing a report

A report is a narrative walk-through of the results, written by hand in Markdown and
rendered by `rules/reports.smk` to a page of the documentation site. This file describes
how to write one; what any particular report says is in the report itself.

Each report is one `.md` file here, named for its key in the `reports` block of
[config.yml](../../config.yml). Adding a report means adding that key and the file.

## Headings

Open with a single `#` heading, the title, on the first line; it also becomes the page's
browser title. `##` and `###` are sections and subsections, and appear in the sidebar.
`####` and deeper render as headings but stay out of the sidebar, which would otherwise
become a wall of links.

## Links

Write links with one of these prefixes. Anything else — including a plain relative link —
fails the build, so a dead link cannot ship.

| To link to | Write | Example |
| --- | --- | --- |
| another page of this site | `docs:` + the path the pipeline writes | `[the titers](docs:results/titer_plots/human_H1N1_recent_individual_sera_subclade.html)` |
| a file or directory on GitHub | `repo:` + its path in this repository | `[the plate layouts](repo:data/plates)` |
| anywhere else | the full URL | `[the H1N1 tree](https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1)` |
| a section of this report | `#` + its slug | `[above](#writing-a-report)` |

A `repo:` target must be tracked by git, since an untracked one would 404 on GitHub. Most
of `results/` is deliberately untracked, so link to a chart with `docs:`, not `repo:`.

`repo:` links point at the repository's default branch, so they follow it rather than
naming a branch that may later be renamed or deleted.

## The repository URL

`{repo_url}` is replaced anywhere in a report by the URL configured for it, so the URL is
written once in [config.yml](../../config.yml) rather than copied into each report:

```markdown
The code and data are at <{repo_url}>, and the open issues are at {repo_url}/issues.
```

Angle brackets around it make the URL itself the link text; `[the code]({repo_url})`
links other words to it instead. Only that exact token is replaced, so any other braces
a report writes are left alone.

## Embedding a plot or a tree

To show something inline rather than only linking to it, write an image on a line of its
own with an `embed:` target, which is either a `docs:`-style path or a full URL:

```markdown
![Median titer of each serum against recent H1N1 strains.](embed:results/titer_plots/human_H1N1_recent_individual_sera_subclade.html)

![The H1N1 tree.](embed:https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1?onlyPanels&d=tree,measurements&sidebar=closed)
```

The caption is shown beneath the frame, followed by a link that opens the same page in a
new tab.

Every frame is as wide as the page and is fitted to the height of the chart it holds, so
an embed takes no sizing of its own. A frame starts a window tall and, once its page has
loaded, shrinks or grows to that page's own height, refitting whenever the page changes
height -- a chart redrawn by one of its own controls, a window resized.

Fitting a frame means reading the page inside it, which a browser allows only for a page
from the same site. A frame keeps its window's height for an embed of another site, such
as a Nextstrain tree, and for every embed when a report is opened from disk rather than
served -- a further reason to preview by serving the site, as below. For a Nextstrain
tree, `onlyPanels` drops the site header and footer, and `d=` and `sidebar=closed` choose
what is shown.

Embeds are loaded lazily, but each chart is still 250 KB to 8 MB, so embed the few that
carry the argument and link to the rest.

## Showing a table

A CSV written by a rule is shown as a table with a `table:` target on a line of its own:

```markdown
![The human sera tested.](table:results/final_titer_data/human_sera_summary.csv){tfoot=1}
```

Write the path within the repository. The CSV is read at build time and rendered as a
real table in the page -- not an iframe -- so the report never recomputes anything it
shows, the report is rebuilt whenever the table changes, and the build fails if the path
is wrong. A CSV shown this way is written to be read: format its numbers and combine a
minimum and maximum into one range column in the rule that writes it, not here.

`{tfoot=<n>}` puts the last *n* rows in the table's footer, where they are shown in
bold. That is how a row totalling the rows above it is marked; the CSV itself simply
holds that row last.

## Showing a static figure

An SVG or PNG in this repository is shown the same way as a plot, with a `figure:`
target on a line of its own:

```markdown
![Public H3N2 HA1 sequences over time.](figure:non-pipeline_analyses/flu_circulating_frequencies/results/plots/H3N2_HA1_counts.svg)
```

Write the path within the repository, not a `raw.githubusercontent.com` URL to the same
file. The file is then read at build time and carried inside the report, so the figure
needs nothing fetched over the network, is rebuilt whenever the image changes, and fails
the build if the path is wrong. A full `https://` URL works for an image genuinely
hosted elsewhere, but gets none of that.

Figures are scaled to the width of the text, so they take no height; `.svg`, `.png`,
`.jpg` and `.jpeg` are understood.

## Building and previewing

```bash
snakemake -j1 --sdm conda -- results/reports/<report>_link_check.txt
```

`repo:` targets and link prefixes are checked as the report renders; `docs:` targets are
checked afterwards against the site that was actually built, and the result is written to
`results/reports/<report>_link_check.txt`.

To preview, serve the built site rather than opening the file directly, so that the
embeds and the links between pages behave as they will once published:

```bash
cd results/docs && python -m http.server
```
