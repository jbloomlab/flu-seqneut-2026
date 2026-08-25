"""Shared category taxonomy for HA sequence vs. library mutation-count binning.

A category is either "identical", "{n}_mutation(s)" for 1 <= n <= max_n, or
"more_than_{max_n}_mutations" -- used both for the full-resolution categories
computed by count_matches.py (max_n = max_mutations_computed) and, reused with
max_n = plot_mutation_threshold, for the coarser categories shown in the plots.
"""


def category_label(n, max_n):
    """n is 0 (identical), 1..max_n, or None/greater than max_n (catch-all)."""
    if n is None or n > max_n:
        return f"more_than_{max_n}_mutations"
    if n == 0:
        return "identical"
    return f"{n}_mutation" if n == 1 else f"{n}_mutations"


def mutation_count(category, max_n):
    """Inverse of category_label: 0, 1..max_n, or max_n + 1 for the catch-all."""
    if category == "identical":
        return 0
    if category == f"more_than_{max_n}_mutations":
        return max_n + 1
    return int(category.split("_", 1)[0])


def pretty_label(category):
    """Human-readable text for legends/tooltips, e.g. '1 mutation different'."""
    if category == "identical":
        return "identical"
    if category.startswith("more_than_"):
        n = category.removeprefix("more_than_").removesuffix("_mutations")
        return f">{n} mutations different"
    n = category.split("_", 1)[0]
    return f"{n} mutation different" if n == "1" else f"{n} mutations different"


def all_categories(max_n):
    """Full ordered category list from identical to the >max_n catch-all."""
    return [category_label(n, max_n) for n in range(max_n + 1)] + [
        category_label(None, max_n)
    ]
