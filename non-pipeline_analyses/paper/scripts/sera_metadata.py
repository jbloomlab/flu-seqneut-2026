"""Format the pipeline's summary of the human sera as a table for the manuscript.

The table is written as HTML so it can be pasted straight into the manuscript with its
formatting intact. Only what the pipeline already summarized is shown, restyled for
print: months rather than `YYYY-MM` dates, and a range given only where there is one.
"""

import csv
import html
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

MONTHS = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]
EN_DASH = "–"

# Every collection date falls in this year, so the header states it once and the cells
# give the month alone. A date from another year is an error, not a silent mislabel.
COLLECTION_YEAR = "2026"

# Dropped from the descriptions, since the vaccine names already say which substrate.
SUBSTRATE_NOTES = [" (cell-based)", " (egg-based)"]

# What the row summarizing every serum is called, as it carries no cohort name.
TOTAL_NAME = "All sera"


def month(iso):
    """'2026-06' -> 'Jun', checking the year the header claims."""
    year, mon = iso.split("-")
    if year != COLLECTION_YEAR:
        raise ValueError(f"{iso} is not in {COLLECTION_YEAR}, which the header states")
    return MONTHS[int(mon) - 1]


def description(text):
    """A cohort description with the substrate parentheticals removed."""
    for note in SUBSTRATE_NOTES:
        text = text.replace(note, "")
    return text


def date_cell(median, span):
    """A month, with its range in parentheses only where the range is wider."""
    if not median:
        return ""
    low, high = span.split(" to ")
    if low == high == median:
        return month(median)
    return f"{month(median)} ({month(low)}{EN_DASH}{month(high)})"


def range_cell(median, span):
    """A median with its range in parentheses."""
    if not median:
        return ""
    low, high = span.split(" to ")
    return f"{median} ({low}{EN_DASH}{high})"


def percent_cell(percent):
    """A percentage, written with its sign."""
    return f"{percent}%" if percent else ""


with open(snakemake.input.summary_csv, newline="") as f:
    rows = list(csv.DictReader(f))
print(f"read {len(rows)} rows from {snakemake.input.summary_csv}")

HEADERS = [
    "Cohort",
    "Description",
    "Individuals",
    "Sera",
    f"Collection month in {COLLECTION_YEAR}, median (range)",
    "Age in years, median (range)",
    "Days post-vaccination, median (range)",
    "Vaccinated in prior year",
]
NUMERIC = {2, 3, 7}  # columns right-aligned

body = []
for row in rows:
    # The summary row carries no cohort name; it is set off as a total instead.
    total = not row["cohort"]
    cells = [
        TOTAL_NAME if total else row["cohort"],
        "" if total else description(row["description"]),
        row["n_individuals"],
        row["n_sera"],
        date_cell(row["collection_date_median"], row["collection_date_range"]),
        range_cell(row["age_median"], row["age_range"]),
        range_cell(
            row["days_post_vaccination_median"], row["days_post_vaccination_range"]
        ),
        # Left empty for the total, where a percentage over every serum would pool the
        # cohorts that record prior-year vaccination with those that record none, and
        # for a cohort that records none, which is not the same as none of its donors
        # having been vaccinated.
        "" if total else percent_cell(row["percent_vaccinated_in_prior_year"]),
    ]
    style = "total" if total else "row"
    tds = []
    for i, cell in enumerate(cells):
        classes = [style] + (["num"] if i in NUMERIC else [])
        tds.append(f'<td class="{" ".join(classes)}">{html.escape(cell)}</td>')
    body.append("      <tr>\n        " + "\n        ".join(tds) + "\n      </tr>")

ths = "\n        ".join(
    f'<th class="head{" num" if i in NUMERIC else ""}">{html.escape(h)}</th>'
    for i, h in enumerate(HEADERS)
)

with open(snakemake.output.table_html, "w", encoding="utf-8") as f:
    f.write(f"""<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <title>Human sera used in this study</title>
    <style>
      body {{ font-family: Arial, Helvetica, sans-serif; font-size: 9pt; }}
      table {{ border-collapse: collapse; }}
      th, td {{ padding: 4pt 7pt; vertical-align: top; }}
      .head {{
        font-weight: bold;
        text-align: left;
        border-top: 1.5pt solid #000000;
        border-bottom: 0.75pt solid #000000;
        background-color: #f2f2f2;
      }}
      .row {{ border-bottom: 0.5pt solid #d9d9d9; }}
      .total {{ font-weight: bold; border-top: 0.75pt solid #000000; }}
      .num {{ text-align: right; }}
    </style>
  </head>
  <body>
    <table>
      <tr>
        {ths}
      </tr>
{chr(10).join(body)}
    </table>
  </body>
</html>
""")

print(f"wrote {len(rows)} rows to {snakemake.output.table_html}")
