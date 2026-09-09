"""Summarize the sera in each sera set as a table, for a report to inline with `table:`.

The sera sets are cohorts of the multicohort CSV, named and ordered by
`summarize_sera.cohort_names` in the config, and must hold each serum exactly once so
that the final row really is the total of the rows above it.

"""

import sys

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

#: How many offending sera an error lists before trailing off
N_SERA_IN_ERROR = 5

#: Decimal places the numbers in the table are rounded to. The ages of the cohorts that
#: release an age range rather than an exact age are the midpoints of those ranges, so
#: any more digits than this would be spurious precision.
DECIMALS = 1


def fmt(number):
    """One number as it is written in the table."""
    return f"{round(number, DECIMALS):g}"


def numeric_stats(values):
    """Median and range of a numeric series, as strings."""
    return {
        "median": fmt(values.median()),
        "range": f"{fmt(values.min())} to {fmt(values.max())}",
    }


def month_stats(values):
    """Median and range of a series of `YYYY-MM` dates, as strings.

    The median is taken as the nearest quantile of the months themselves, so it is
    always a month that was actually observed rather than one halfway between two.

    """
    months = pd.PeriodIndex(values, freq="M")
    ordinals = pd.Series(months.astype("int64"))
    median = pd.Period(
        ordinal=int(ordinals.quantile(0.5, interpolation="nearest")), freq="M"
    )
    return {"median": str(median), "range": f"{months.min()} to {months.max()}"}


def summarize(description, cohort, sera):
    """One row of the table, summarizing a frame holding one row per serum."""
    row = {
        "description": description,
        "cohort": cohort,
        # a serum with no `subject_id` counts as its own individual: the cohorts that
        # record none have one serum per donor
        "n_individuals": sera["subject_id"].nunique() + sera["subject_id"].isna().sum(),
        "n_sera": len(sera),
    }
    for key, value in month_stats(sera["serum_collection_date"]).items():
        row[f"collection_date_{key}"] = value
    for key, value in numeric_stats(sera["age_numeric"]).items():
        row[f"age_{key}"] = value
    # Left empty unless every serum records it and at least one was drawn after a
    # vaccination: a median over a mix of pre-vaccination sera (all 0 days) and
    # post-vaccination sera would mean nothing, and for the pre-vaccination sera on
    # their own the number just restates that they are pre-vaccination.
    days = sera["days_post_vaccination"]
    stats = numeric_stats(days) if days.notna().all() and days.any() else {}
    for key in ["median", "range"]:
        row[f"days_post_vaccination_{key}"] = stats.get(key, "")
    # over just the sera that record it, so a set recording none is left empty rather
    # than reading as zero percent
    vaccinated = sera[prior_year_vaccination_column].dropna()
    row["percent_vaccinated_in_prior_year"] = (
        fmt(100 * (vaccinated == "Yes").mean()) if len(vaccinated) else ""
    )
    return row


cohort_names = snakemake.params.config["cohort_names"]
total_name = snakemake.params.config["total_name"]
prior_year_vaccination_column = snakemake.params.config["prior_year_vaccination_column"]

sera = pd.read_csv(
    snakemake.input.sera_multicohort, dtype={"serum": str, "subject_id": str}
)
print(f"Read {len(sera)} rows from {snakemake.input.sera_multicohort}")

if prior_year_vaccination_column not in sera.columns:
    raise ValueError(
        "`summarize_sera.prior_year_vaccination_column` names a column not in "
        f"{snakemake.input.sera_multicohort}: {prior_year_vaccination_column}"
    )
invalid = set(sera[prior_year_vaccination_column].dropna()) - {"Yes", "No"}
if invalid:
    raise ValueError(
        f"`{prior_year_vaccination_column}` in {snakemake.input.sera_multicohort} must "
        f"be `Yes` or `No`, but also has {sorted(invalid)}"
    )

unknown = [cohort for cohort in cohort_names if cohort not in set(sera["cohort"])]
if unknown:
    raise ValueError(
        f"`summarize_sera.cohort_names` names cohorts not in "
        f"{snakemake.input.sera_multicohort}: {unknown}"
    )

listed = sera[sera["cohort"].isin(cohort_names)]
n_sets = listed.groupby("serum").size().reindex(sera["serum"].unique(), fill_value=0)
if not (n_sets == 1).all():
    unassigned = sorted(n_sets.index[n_sets == 0])
    repeated = sorted(n_sets.index[n_sets > 1])
    raise ValueError(
        "the cohorts in `summarize_sera.cohort_names` must hold each serum exactly "
        f"once, but {len(unassigned)} sera are in none of them "
        f"({unassigned[:N_SERA_IN_ERROR]}) and {len(repeated)} are in more than one "
        f"({repeated[:N_SERA_IN_ERROR]})"
    )

rows = [
    summarize(name, cohort, sera[sera["cohort"] == cohort])
    for (cohort, name) in cohort_names.items()
]
# the listed cohorts hold each serum exactly once, so this really is their total
rows.append(summarize(total_name, "", sera.drop_duplicates(subset="serum")))

summary = pd.DataFrame(rows)
summary.to_csv(snakemake.output.csv, index=False)
print(f"\nWrote {snakemake.output.csv}:\n{summary.to_string(index=False)}")
