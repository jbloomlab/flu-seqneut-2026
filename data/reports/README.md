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
links other words to it instead. Only that exact token is replaced, so braces a report
writes for other reasons -- an embed's `{height=}` -- are left alone.

## Embedding a plot or a tree

To show something inline rather than only linking to it, write an image on a line of its
own with an `embed:` target, which is either a `docs:`-style path or a full URL:

```markdown
![Median titer of each serum against recent H1N1 strains.](embed:results/titer_plots/human_H1N1_recent_individual_sera_subclade.html)

![The H1N1 tree.](embed:https://nextstrain.org/community/jbloomlab/flu-seqneut-2026@main/H1N1?onlyPanels&d=tree,measurements&sidebar=closed){height=900}
```

The caption is shown beneath the frame, followed by a link that opens the same page in a
new tab.

A page of this site is measured in the browser, so its frame follows the chart as it
grows and shrinks — which matters for the charts whose height changes when serum cohorts
are toggled in the legend. Anything on another site cannot be measured, so give it a
`{height=<px>}`; for a Nextstrain tree, `onlyPanels` drops the site header and footer,
and `d=` and `sidebar=closed` choose what is shown.

Embeds are loaded lazily, but each chart is still 250 KB to 8 MB, so embed the few that
carry the argument and link to the rest.

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

To preview, serve the built site rather than opening the file directly:

```bash
cd results/docs && python -m http.server
```

Opened over `file://`, every embed counts as another site and falls back to a fixed
height, so the resizing cannot be seen.
