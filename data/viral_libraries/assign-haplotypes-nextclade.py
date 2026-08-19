"""Script to assign `subclade` and `derived_haplotype` in the designed library CSV.

Both columns are derived from `nt_sequence_HA_ectodomain` with `nextclade`, rather
than carried over from the library-design analysis, which supplied them from
hand-curated haplotype tables and so left some strains labelled with a bare subclade.

`subclade` is nextclade's subclade call. `derived_haplotype` is that subclade plus the
amino-acid mutations separating the strain from its subclade founder, HA1 mutations
first and bare, then HA2 mutations prefixed `HA2_` (e.g. `K:F192V,HA2_R32K`).

Run by hand from this directory, after `create-designed-library-csv.py` and before
`create-actual-library-csv.py`; see README.md.
"""

import datetime
import json
import pathlib
import re
import subprocess
import tempfile

import pandas as pd

input_csv = "flu-seqneut-2026-barcode-to-strain-designed.csv"
output_csv = "flu-seqneut-2026-barcode-to-strain-designed.csv"
provenance_txt = "assign-haplotypes-nextclade_provenance.txt"

# nextclade dataset per `subtype`. Each is the HA dataset for its subtype and provides
# the `subclade` calls whose nomenclature this library uses. The default references are
# used rather than the `_broad` variants (A/Wisconsin/67/2005, A/California/7/2009)
# because the oldest strains here are from 2019 and assign cleanly against these. A new
# subtype or lineage needs a new entry; note pre-2009 seasonal H1N1 is not h1n1pdm and
# its dataset offers no clade calls at all, so it would need a different approach.
datasets = {
    "H3N2": "nextstrain/flu/h3n2/ha/EPI1857216",
    "H1N1": "nextstrain/flu/h1n1pdm/ha/MW626062",
}

# columns this script rewrites
assigned_columns = ["subclade", "derived_haplotype"]

# genes the haplotype is built from, in the order their mutations are listed. Mutations
# in any other gene are an error, as the naming convention has no place to put them.
haplotype_genes = ["HA1", "HA2"]

# prefix applied to mutations of each gene in `derived_haplotype`; HA1 mutations are
# bare, so that the labels this produces match the pre-existing HA1-only convention.
gene_prefixes = {"HA1": "", "HA2": "HA2_"}

# Minimum fraction of each gene the sequences must cover, below which the assignment is
# not trustworthy. HA1 is short of 1 for H1N1 because the constructs omit the first three
# ectodomain residues. HA2 is short of 1 for both because the constructs stop before the
# transmembrane and cytoplasmic tail, which is constant across strains and so carries no
# strain-derived differences to detect.
min_cds_coverage = {"HA1": 0.98, "HA2": 0.75}

# HA1 sites the H1N1 constructs do not cover, so no mutation may be called in them
uncovered_ha1_sites = 3

mutation_regex = re.compile(r"^([A-Za-z0-9]+):([A-Z*-])(\d+)([A-Z*-])$")


def read_strains(csv_path):
    """Read `csv_path` down to one row per strain, validating strain consistency."""
    df = pd.read_csv(csv_path, dtype=str)
    strains = df.drop_duplicates(["strain", "subtype", "nt_sequence_HA_ectodomain"])
    if len(strains) != strains["strain"].nunique():
        raise ValueError(
            "Strains with more than one subtype or sequence: "
            f"{strains.loc[strains['strain'].duplicated(), 'strain'].tolist()}"
        )
    unknown = sorted(set(strains["subtype"]) - set(datasets))
    if unknown:
        raise ValueError(f"No nextclade dataset configured for subtype(s): {unknown}")
    return df, strains


def run_nextclade(strains, subtype, workdir):
    """Return nextclade's output frame and dataset tag for one subtype's strains.

    The dataset is fetched without a tag, so the latest is always used; it is fetched
    to disk rather than passed to `run --dataset-name` so the tag it resolved to can be
    read back and recorded.
    """
    subtype_strains = strains[strains["subtype"] == subtype]
    fasta = workdir / f"{subtype}.fasta"
    fasta.write_text(
        "".join(
            f">{row['strain']}\n{row['nt_sequence_HA_ectodomain']}\n"
            for _, row in subtype_strains.iterrows()
        )
    )

    dataset_dir = workdir / f"dataset_{subtype}"
    subprocess.run(
        ["nextclade", "dataset", "get", "--name", datasets[subtype]]
        + ["--output-dir", str(dataset_dir)],
        check=True,
    )
    tag = json.loads((dataset_dir / "pathogen.json").read_text())["version"]["tag"]

    tsv = workdir / f"nextclade_{subtype}.tsv"
    subprocess.run(
        ["nextclade", "run", "--input-dataset", str(dataset_dir)]
        + ["--output-tsv", str(tsv), "--in-order", "--silent", str(fasta)],
        check=True,
    )

    nextclade = pd.read_csv(tsv, sep="\t", dtype=str)
    missing = set(subtype_strains["strain"]) - set(nextclade["seqName"])
    if missing:
        raise ValueError(
            f"nextclade returned no result for {subtype}: {sorted(missing)}"
        )
    return nextclade, tag


def check_nextclade_output(nextclade, subtype):
    """Fail fast on anything that would make the assignment untrustworthy."""
    failed = nextclade[nextclade["errors"].notnull()]
    if len(failed):
        raise ValueError(
            f"nextclade errors for {subtype}:\n{failed[['seqName', 'errors']]}"
        )

    bad_qc = nextclade[nextclade["qc.overallStatus"] != "good"]
    if len(bad_qc):
        raise ValueError(
            f"nextclade QC not good for {subtype}:\n"
            f"{bad_qc[['seqName', 'qc.overallStatus']]}"
        )

    # these datasets make `subclade` the default clade, so the two must agree; if they
    # ever diverge the dataset has changed in a way this script's naming assumes it has not
    disagree = nextclade[nextclade["clade"] != nextclade["subclade"]]
    if len(disagree):
        raise ValueError(
            f"nextclade `clade` and `subclade` disagree for {subtype}:\n"
            f"{disagree[['seqName', 'clade', 'subclade']]}"
        )

    deletions = nextclade[nextclade["founderMuts['clade'].aaDeletions"].notnull()]
    if len(deletions):
        raise ValueError(
            f"{subtype} strains have amino-acid deletions from their subclade founder, "
            "which `derived_haplotype` has no notation for:\n"
            f"{deletions['seqName'].tolist()}"
        )

    for gene, minimum in min_cds_coverage.items():
        coverage = (
            nextclade["cdsCoverage"]
            .str.extract(rf"{gene}:([\d.]+)", expand=False)
            .astype(float)
        )
        if coverage.isnull().any():
            raise ValueError(f"No {gene} coverage reported for {subtype}")
        if (coverage < minimum).any():
            below = nextclade.loc[coverage < minimum, "seqName"].tolist()
            raise ValueError(
                f"{subtype} strains below the {gene} coverage floor of {minimum}: {below}"
            )


def parse_founder_mutations(value, strain):
    """Return {gene: [(site, mutation)]} parsed from a `founderMuts` cell."""
    by_gene = {gene: [] for gene in haplotype_genes}
    if pd.isnull(value):
        return by_gene
    for mutation in value.split(","):
        match = mutation_regex.match(mutation)
        if not match:
            raise ValueError(f"Cannot parse mutation {mutation!r} for {strain}")
        gene, wildtype, site, mutant = match.groups()
        if gene not in by_gene:
            raise ValueError(
                f"{strain} has a mutation in {gene} ({mutation}), but `derived_haplotype` "
                f"only names mutations in {haplotype_genes}"
            )
        site = int(site)
        if gene == "HA1" and site <= uncovered_ha1_sites:
            raise ValueError(
                f"{strain} has a mutation at HA1 site {site}, which the constructs do "
                "not cover, so it cannot have been observed"
            )
        by_gene[gene].append((site, f"{wildtype}{site}{mutant}"))
    return by_gene


def build_haplotype(subclade, by_gene):
    """Return the `derived_haplotype` string for one strain."""
    mutations = [
        gene_prefixes[gene] + mutation
        for gene in haplotype_genes
        for _, mutation in sorted(by_gene[gene])
    ]
    return subclade + (":" + ",".join(mutations) if mutations else "")


def assign(strains):
    """Return a frame of strain, subclade, derived_haplotype, plus the dataset tags."""
    assignments = []
    tags = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = pathlib.Path(tmpdir)
        for subtype in sorted(strains["subtype"].unique()):
            nextclade, tags[subtype] = run_nextclade(strains, subtype, workdir)
            check_nextclade_output(nextclade, subtype)
            for _, row in nextclade.iterrows():
                by_gene = parse_founder_mutations(
                    row["founderMuts['clade'].aaSubstitutions"], row["seqName"]
                )
                assignments.append(
                    {
                        "strain": row["seqName"],
                        "subclade": row["subclade"],
                        "derived_haplotype": build_haplotype(row["subclade"], by_gene),
                    }
                )

    assigned = pd.DataFrame(assignments)
    implied = assigned["derived_haplotype"].str.split(":").str[0]
    if not (assigned["subclade"] == implied).all():
        raise ValueError("`subclade` does not match the `derived_haplotype` prefix")
    return assigned, tags


def report(strains, assigned, tags):
    """Print what nextclade assigned and how it differs from the current values."""
    print(f"Read {len(strains)} strains from {input_csv}\n")
    for subtype, tag in sorted(tags.items()):
        print(f"{subtype}: {datasets[subtype]} (tag {tag})")
    version = subprocess.run(
        ["nextclade", "--version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    print(f"nextclade version: {version}\n")

    comparison = strains[["strain", *assigned_columns]].merge(
        assigned, on="strain", how="left", validate="one_to_one", suffixes=("_old", "")
    )

    changed_subclade = comparison[comparison["subclade_old"] != comparison["subclade"]]
    print(f"=== subclade: {len(changed_subclade)} of {len(comparison)} changed ===")
    if len(changed_subclade):
        print(
            changed_subclade[["strain", "subclade_old", "subclade"]].to_string(
                index=False
            )
        )
    print()

    changed_haplotype = comparison[
        comparison["derived_haplotype_old"] != comparison["derived_haplotype"]
    ]
    print(
        f"=== derived_haplotype: {len(changed_haplotype)} of {len(comparison)} changed ==="
    )
    if len(changed_haplotype):
        print(
            changed_haplotype[
                ["strain", "derived_haplotype_old", "derived_haplotype"]
            ].to_string(index=False)
        )
    print()

    # A duplicate now means two strains whose strain-derived ectodomain is identical
    # across both HA1 and HA2, which is a stronger statement than under the previous
    # HA1-only convention. Reported rather than raised, for a human to judge.
    counts = assigned["derived_haplotype"].value_counts()
    duplicated = counts[counts > 1]
    print(
        f"=== derived_haplotype shared by more than one strain: {len(duplicated)} ==="
    )
    for haplotype in duplicated.index:
        sharing = assigned.loc[
            assigned["derived_haplotype"] == haplotype, "strain"
        ].tolist()
        print(f"  {haplotype}: {sharing}")
    print()


def write_provenance(tags):
    """Record what produced the assignment, since the dataset is always the latest."""
    version = subprocess.run(
        ["nextclade", "--version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    lines = [
        f"# written by {pathlib.Path(__file__).name}",
        f"date: {datetime.datetime.now(tz=datetime.UTC).date().isoformat()}",
        f"nextclade: {version}",
    ]
    lines += [
        f"{subtype}: {datasets[subtype]} {tags[subtype]}" for subtype in sorted(tags)
    ]
    pathlib.Path(provenance_txt).write_text("\n".join(lines) + "\n")
    print(f"Wrote {provenance_txt}")


if __name__ == "__main__":
    df, strains = read_strains(input_csv)
    assigned, tags = assign(strains)
    report(strains, assigned, tags)

    updated = df.drop(columns=assigned_columns).merge(
        assigned, on="strain", how="left", validate="many_to_one"
    )
    if updated[assigned_columns].isnull().any(axis=None):
        raise ValueError("Some rows were not assigned a subclade or derived_haplotype")
    updated = updated[df.columns]

    write_provenance(tags)
    updated.to_csv(output_csv, index=False)
    print(f"Wrote {len(updated)} rows to {output_csv}")
