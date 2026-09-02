"""Regenerate the five coloring CSVs for the 8FAW antigenic-regions example.

The CSVs are committed and no rule runs this script -- ``spec.yaml`` reads the
committed files. Run it by hand from the repository root when the inputs below
change::

    python data/prot-struct-viz_config/8faw_antigenic_regions/make_coloring_csv.py

It writes one CSV per view of ``spec.yaml``, each named after the view that
reads it -- so a name here that no view claims is a name that has gone stale:

* ``antigenic-regions-w-glycans.csv`` -- every residue, HA1 colored by antigenic
  region, plus the host N-glycans and the LSTc receptor analogue;
* ``antigenic-regions.csv`` -- the same without the host-glycan rows, so that the
  view's ``glycans: hide`` can take them away;
* ``perth-2009-to-subclade-k.csv`` and ``2025-26-to-2026-27-vaccine.csv`` -- the
  sites that differ between two HAs, in red, plus the same LSTc rows;
* ``subclade-k-with-region-d-mutations.csv`` -- the second of those two lists again,
  with HA1 222 and 223 added in a color of their own.

It exists because the ~490 rows of the first two are derived rather than typed,
and a table that large is only auditable if its derivation ships with it. Two
things it reads rather than hard-codes:

* this repository's H3N2 site-numbering map, which is what turns a residue
  number into an HA1 or HA2 site number;
* PDB 8FAW, fetched from RCSB, so that only *modeled* residues get a row.

The two mutation lists are the exception: they are transcribed below rather than
recomputed from the sequence libraries they came from. That is inherited from
the upstream example this directory was copied from, where those libraries were
not available; here they are, so recomputing the lists from
``data/viral_libraries/`` would drop the transcription step entirely. Until that
is done, what guards them is this script's own checks: the mapping from an
HA1/HA2 site to a residue of 8FAW, the expected length of each list, and that
every site it names is actually modeled.

The numbering frame is the load-bearing assumption throughout, so it is asserted
rather than trusted: if a revised map or a different entry ever shifts it, this
script fails instead of writing a plausible-looking but wrong CSV.
"""

import csv
import pathlib
import urllib.request

import gemmi

#: This repository's H3N2 site-numbering map, located from this script rather
#: than from the working directory. Author numbering in 8FAW is a single chain
#: running across HA1 and HA2, which is the same frame as this map's
#: ``sequential_site``.
NUMBERING_MAP = (
    pathlib.Path(__file__).parents[3]
    / "data/nextstrain-prot-titers-tree_data/H3N2_site_numbering_map.tsv"
)

PDB_ID = "8FAW"
STRUCTURE_URL = f"https://files.rcsb.org/download/{PDB_ID}.cif"

#: The polymer. 8FAW's asymmetric unit is one protomer; the trimer is assembly 1.
POLYMER_CHAIN = "A"

#: The LSTc receptor analogue, Neu5Ac-a2,6-Gal-b1,4-GlcNAc-b1,3-Gal-b1,4-Glc,
#: reducing end first. Naming a sugar in the CSV is what replaces its 3D-SNFG
#: symbols with a plain colored ball-and-stick. Every sugar in the entry is named
#: here -- the receptor analogue and the host N-glycans alike -- so nothing is
#: left to a view's `glycans` setting.
RECEPTOR_CHAIN = "D"

#: Antigenic sites A-E of H3 HA1, from Table 2 ("Amino acids assigned to
#: antigenic sites") of Stray & Pittman, Virol J 2012;9:91, PMC3499391.
#: Footnote h of that table: "H3 HA residues numbered as for mature HA1 of
#: A/Aichi/2/1968". Footnote g: residues the table encloses in parentheses were
#: assigned in previous studies but did not meet the authors' own inclusion
#: criteria. Both kinds are included here -- this is the conventional full
#: definition of the five sites, and the distinction is about how the 2012
#: authors scored evolutionary rates, not about whether antibodies bind there.
SITES = {
    "A": [
        121,
        122,
        123,
        124,
        125,
        126,
        127,
        129,
        131,
        132,
        133,
        134,
        135,
        136,
        137,
        138,
        140,
        142,
        143,
        144,
        145,
        146,
    ],
    "B": [
        155,
        156,
        157,
        158,
        159,
        160,
        186,
        188,
        189,
        190,
        192,
        193,
        194,
        196,
        197,
        198,
        199,
        246,
        247,
    ],
    "C": [49, 50, 53, 54, 271, 273, 275, 276, 278],
    "D": [
        167,
        201,
        202,
        203,
        204,
        205,
        206,
        207,
        214,
        216,
        217,
        218,
        219,
        220,
        222,
        223,
        225,
        226,
        227,
        242,
    ],
    "E": [62, 63, 75, 78, 79, 80, 81, 82, 83, 91, 92, 94],
}

#: Colorblind-safe qualitative palette. Deliberately not the colors of the
#: paper's own figure, whose yellow site B would be illegible against the pale
#: gray used for the rest of the structure.
# Paul Tol's "muted" scheme, cool subset. Every one of these is cool, so on the
# structure warm means "changed or added" -- the red of a mutated site, the gold
# of a glycan -- and nothing else is close to MUTATED_COLOR.
SITE_COLORS = {
    "A": "#332288",
    "B": "#88CCEE",
    "C": "#117733",
    "D": "#AA4499",
    "E": "#44AA99",
}

#: The two grays are not decoration: they draw the HA1/HA2 boundary, which is
#: what the _HA1 and _HA2 label suffixes exist to tell you about.
HA1_COLOR = "#e8e8e8"
HA2_COLOR = "#bdbdbd"
RECEPTOR_COLOR = "#000000"
GLYCAN_COLOR = "#ffd700"

#: Sites that differ between the two HAs of a comparison. Warm, where every
#: `SITE_COLORS` entry is cool: the last three views drop the antigenic-region
#: coloring, and the one of them that puts a region color back has to hold it
#: apart from this.
MUTATED_COLOR = "#e41a1c"

#: HA1 222 and 223, drawn on top of a mutation list in one view. They are two of
#: the twenty residues `SITES` assigns to antigenic region D, and they sit in the
#: 220-loop, one of the three elements lining the receptor-binding site.
REGION_D_SITES = [222, 223]

#: The purple that the first two views already give antigenic region D, so the
#: pair reads as the same thing here as it does there. That is the whole reason
#: this is not a new color: the other views' key is the key for this one too.
#: It is a cool color against `MUTATED_COLOR`, which is the other requirement --
#: this is the only view drawing two colored classes of site at once.
REGION_D_COLOR = SITE_COLORS["D"]

CITATION = "Stray & Pittman 2012 Virol J 9:91"

#: How each HA is named in a tooltip and in the CSV's ``notes`` column, with the
#: GenBank accession of the sequence actually compared.
PERTH = "A/Perth/16/2009 (GQ293081)"
SUBCLADE_K = "subclade K (A/Darwin/1415/2025, PX422923)"
DC_2023 = "A/District_Of_Columbia/27/2023 (PV280355)"
DARWIN_2025 = "A/Darwin/1415/2025 (PX422923)"

#: Sites differing between two HAs, as ``(protein, site, from, to)``.
#:
#: Both lists were computed once from the three HA ectodomain protein sequences
#: named above -- as they appear in the Bloom lab's sequencing libraries, which
#: is why each is pinned to an accession here -- and then transcribed. Subclade K
#: is represented by the one strain the 2026 library assigns ``derived_haplotype
#: == "K"``, A/Darwin/1415/2025, which is also the 2026-2027 vaccine strain: that
#: is why the two comparisons share a target and differ only in the comparator.
#:
#: All three ectodomains are 501 residues with no indels, so they compare
#: position by position with no alignment, and ectodomain position is HA1 site
#: 1-329 then HA2 site (position - 329) -- the same frame `FRAME_CHECKS` pins.
#: Substitutions outside what 8FAW models are dropped rather than listed:
#: ``K2N`` and ``L3I`` fall before HA1 11, and ``T328A`` after HA1 325.
PERTH_TO_SUBCLADE_K = [
    ("HA1", 33, "Q", "R"),
    ("HA1", 45, "S", "N"),
    ("HA1", 48, "T", "I"),
    ("HA1", 50, "E", "K"),
    ("HA1", 53, "D", "N"),
    ("HA1", 62, "K", "G"),
    ("HA1", 83, "K", "E"),
    ("HA1", 92, "K", "R"),
    ("HA1", 94, "Y", "N"),
    ("HA1", 96, "N", "S"),
    ("HA1", 121, "N", "K"),
    ("HA1", 122, "N", "D"),
    ("HA1", 131, "T", "K"),
    ("HA1", 135, "T", "K"),
    ("HA1", 140, "I", "K"),
    ("HA1", 142, "R", "G"),
    ("HA1", 144, "K", "N"),
    ("HA1", 145, "N", "S"),
    ("HA1", 156, "H", "S"),
    ("HA1", 158, "N", "D"),
    ("HA1", 159, "F", "N"),
    ("HA1", 164, "L", "Q"),
    ("HA1", 171, "N", "K"),
    ("HA1", 173, "Q", "R"),
    ("HA1", 186, "G", "D"),
    ("HA1", 189, "K", "R"),
    ("HA1", 190, "D", "N"),
    ("HA1", 192, "I", "F"),
    ("HA1", 193, "F", "S"),
    ("HA1", 195, "Y", "F"),
    ("HA1", 198, "A", "S"),
    ("HA1", 212, "T", "A"),
    ("HA1", 214, "S", "I"),
    ("HA1", 225, "N", "D"),
    ("HA1", 276, "K", "E"),
    ("HA1", 278, "N", "K"),
    ("HA1", 312, "N", "S"),
    ("HA2", 77, "I", "V"),
    ("HA2", 155, "G", "E"),
    ("HA2", 160, "D", "N"),
]

DC_2023_TO_DARWIN_2025 = [
    ("HA1", 135, "T", "K"),
    ("HA1", 144, "S", "N"),
    ("HA1", 145, "N", "S"),
    ("HA1", 158, "N", "D"),
    ("HA1", 160, "I", "K"),
    ("HA1", 173, "Q", "R"),
    ("HA1", 189, "K", "R"),
    ("HA2", 49, "S", "N"),
]

#: Expected lengths of the two lists above, asserted for the same reason
#: `check_sites` asserts the antigenic sites' size: a hand-transcribed list is
#: worth nothing if a dropped line goes unnoticed.
MUTATION_COUNTS = {"PERTH_TO_SUBCLADE_K": 40, "DC_2023_TO_DARWIN_2025": 8}

#: Residue number of HA2 site 1 in 8FAW's author numbering. HA1 sites are author
#: numbers already, so HA1 needs no offset.
HA2_OFFSET = 329

#: Residues whose identity fixes the numbering frame. The four aromatics line
#: the receptor-binding site and the nine cysteines form HA1's disulfides; all
#: are invariant across H3 and all are quoted in H3 HA1 numbering.
FRAME_CHECKS = {
    98: "TYR",
    153: "TRP",
    183: "HIS",
    195: "TYR",
    14: "CYS",
    52: "CYS",
    64: "CYS",
    76: "CYS",
    97: "CYS",
    139: "CYS",
    277: "CYS",
    281: "CYS",
    305: "CYS",
}

#: No view of this example draws labels into the scene, so the three columns
#: that style a drawn label -- ``show_label``, ``label_color``, ``label_size``
#: -- are left out rather than written empty. ``label`` stays: it is the
#: mouseover tooltip, which is how a reader identifies a residue here.
HEADER = [
    "chain",
    "residue",
    "color",
    "label",
    "representation",
    "notes",
]


def load_numbering_map():
    """Return ``{sequential_site: (protein, protein_site)}`` from the lab's map."""
    with NUMBERING_MAP.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            int(row["sequential_site"]): (row["protein"], int(row["protein_site"]))
            for row in reader
        }


def load_structure():
    """Return the 8FAW model, fetched from RCSB."""
    with urllib.request.urlopen(STRUCTURE_URL) as response:
        text = response.read().decode()
    structure = gemmi.make_structure_from_block(
        gemmi.cif.read_string(text).sole_block()
    )
    structure.setup_entities()
    return structure[0]


def partition_residues(model):
    """Split the modeled residues into polymer, host glycans and receptor.

    gemmi splits waters and heteroatoms into a second chain reusing the author
    ID, so the polymer and the two N-glycans that hang off it both appear under
    chain A. Classify by component rather than by chain: amino acids on the
    polymer chain are the protein, everything on the receptor chain is the
    receptor, and every other non-water residue is a host glycan.

    Returns three ``[(chain, residue_number, component_name), ...]`` lists.
    """
    polymer, glycans, receptor = [], [], []
    for chain in model:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info is not None and info.is_water():
                continue
            entry = (chain.name, residue.seqid.num, residue.name)
            if chain.name == RECEPTOR_CHAIN:
                receptor.append(entry)
            elif (
                chain.name == POLYMER_CHAIN
                and info is not None
                and info.is_amino_acid()
            ):
                polymer.append(entry)
            else:
                glycans.append(entry)
    if not polymer:
        raise SystemExit(f"{PDB_ID} has no polymer on chain {POLYMER_CHAIN!r}")
    if not receptor:
        raise SystemExit(f"{PDB_ID} has no receptor on chain {RECEPTOR_CHAIN!r}")
    return polymer, glycans, receptor


def check_frame(polymer):
    """Fail loudly if author numbering is not H3 HA1 numbering after all."""
    seen = {number: component for _, number, component in polymer}
    wrong = {
        num: (expected, seen.get(num))
        for num, expected in FRAME_CHECKS.items()
        if seen.get(num) != expected
    }
    if wrong:
        raise SystemExit(
            f"{PDB_ID} chain {POLYMER_CHAIN} is not in H3 HA1 numbering: "
            + "; ".join(
                f"expected {exp} at {num}, found {got}"
                for num, (exp, got) in sorted(wrong.items())
            )
        )


def check_sites():
    """Fail loudly if the transcribed site definitions have drifted."""
    flat = [site for residues in SITES.values() for site in residues]
    if len(flat) != len(set(flat)):
        raise SystemExit("antigenic sites overlap; they are disjoint in Table 2")
    if len(flat) != 82:
        raise SystemExit(f"expected 82 antigenic-site residues, got {len(flat)}")


def site_of(residue_number):
    """Return the antigenic site containing this HA1 site, or ``None``."""
    for site, residues in SITES.items():
        if residue_number in residues:
            return site
    return None


def map_ha2_offset(numbering):
    """Return the numbering map's own HA2 offset, checked against `HA2_OFFSET`.

    The map derives it from its own rows; `HA2_OFFSET` is what the hard-coded
    mutation lists assume. They are the same number for the same reason -- HA1
    is 329 residues -- but nothing forces that, so it is checked once here
    rather than trusted in two places.
    """
    offset = max(
        seq for seq, (protein, _) in numbering.items() if protein == "HA2"
    ) - max(site for protein, site in numbering.values() if protein == "HA2")
    if offset != HA2_OFFSET:
        raise SystemExit(
            f"numbering map puts HA2 site 1 at residue {offset + 1}, but the "
            f"mutation lists assume {HA2_OFFSET + 1}"
        )
    return offset


def polymer_rows(polymer, numbering):
    """One row per modeled polymer residue, colored by antigenic site."""
    ha2_offset = map_ha2_offset(numbering)
    rows = []
    for _, number, _component in polymer:
        if number in numbering:
            protein, site = numbering[number]
        else:
            # 8FAW models one residue past the end of the map: author 502, the
            # HA2 C-terminus, which is HA2 site 173. The map's HA2 offset is
            # constant, so applying it here gives that site rather than guessing
            # -- and dropping the residue instead would leave the C-terminus
            # uncolored and unlabeled.
            protein, site = "HA2", number - ha2_offset
            print(
                f"note: residue {number} is past the end of the numbering map; "
                f"it is in HA2, and the map's own HA2 offset of {ha2_offset} "
                f"puts it at {site}_{protein}"
            )
        region = site_of(site) if protein == "HA1" else None
        # One label for every residue, antigenic site or not: the site number
        # alone is ambiguous across HA1 and HA2, and this is the mouseover text.
        label = f"{site}_{protein}"
        if region is not None:
            rows.append(
                [
                    POLYMER_CHAIN,
                    number,
                    SITE_COLORS[region],
                    label,
                    "",
                    f"HA1 site {site}; antigenic site {region} of {CITATION}",
                ]
            )
        else:
            rows.append(
                [
                    POLYMER_CHAIN,
                    number,
                    HA1_COLOR if protein == "HA1" else HA2_COLOR,
                    label,
                    "",
                    f"{protein} site {site}; not in a defined antigenic site",
                ]
            )
    return rows


def receptor_rows(receptor):
    """One row per LSTc sugar, in every CSV.

    Every view keeps the receptor analogue: it is the landmark that says which
    end of a protomer is membrane-distal. Four of the five views set ``glycans:
    hide``, and 8FAW's sugars all classify as glycans -- the host N-glycans and
    these five alike -- so naming them here is the only reason they survive.
    """
    rows = []
    for index, (_, number, component) in enumerate(receptor, start=1):
        last = index == len(receptor)
        rows.append(
            [
                RECEPTOR_CHAIN,
                number,
                RECEPTOR_COLOR,
                "LSTc" if last else component,
                "ball-and-stick",
                (
                    f"LSTc sugar {index} of {len(receptor)} ({component}); "
                    "named here so it is drawn in this color rather than as an "
                    "SNFG symbol, and so that `glycans: hide` does not take it "
                    "away"
                ),
            ]
        )
    return rows


def glycan_rows(glycans):
    """One row per host N-glycan sugar, drawn in one color rather than as SNFG.

    Only the first view's CSV gets these. Leaving them out of the other four is
    what lets those views' ``glycans: hide`` take them away, since a residue
    named in the CSV is drawn whatever the heteroatom options say.
    """
    return [
        [
            chain,
            number,
            GLYCAN_COLOR,
            component,
            "ball-and-stick",
            (
                f"Host N-glycan sugar ({component}); named here so it is drawn "
                "in this color rather than as an SNFG symbol"
            ),
        ]
        for chain, number, component in sorted(glycans)
    ]


def residue_number(protein, site):
    """Return the 8FAW author residue number of an HA1 or HA2 site."""
    if protein == "HA1":
        return site
    if protein == "HA2":
        return site + HA2_OFFSET
    raise SystemExit(f"unknown protein {protein!r}; expected 'HA1' or 'HA2'")


def check_mutations(name, mutations, modeled):
    """Fail loudly if a hard-coded mutation list has drifted or gone unmodeled.

    Both failures would otherwise be silent: a dropped line just loses a red
    site, and a site 8FAW does not model would be reported by ``on_mismatch``
    but only as one line among the rest of the build's output.
    """
    expected = MUTATION_COUNTS[name]
    if len(mutations) != expected:
        raise SystemExit(f"{name} has {len(mutations)} entries, expected {expected}")
    seen = {(protein, site) for protein, site, _, _ in mutations}
    if len(seen) != len(mutations):
        raise SystemExit(f"{name} names the same site twice")
    absent = [
        f"{site}_{protein}"
        for protein, site, _, _ in mutations
        if residue_number(protein, site) not in modeled
    ]
    if absent:
        raise SystemExit(f"{name} names sites 8FAW does not model: {absent}")


def mutation_rows(mutations, from_name, to_name):
    """One red row per differing site, and nothing for the rest of the molecule.

    The antigenic-region coloring is gone from these views: every residue
    without a row here takes the view's ``default_color``, so the only thing
    painted is what changed. The tooltip carries the substitution itself, which
    is the only way to read the view -- no label is drawn into the scene.
    """
    return [
        [
            POLYMER_CHAIN,
            residue_number(protein, site),
            MUTATED_COLOR,
            f"{old}{site}{new}_{protein}",
            "",
            (
                f"{protein} site {site} differs: {old} in {from_name}, "
                f"{new} in {to_name}"
            ),
        ]
        for protein, site, old, new in mutations
    ]


def region_d_rows(mutations, modeled):
    """One row per site of `REGION_D_SITES`, for the view that adds them.

    These are written alongside `mutation_rows` into a single CSV, so all three
    of the things that could make that CSV lie are checked here: that the sites
    really are in region D, that 8FAW models them, and that the mutation list
    does not already claim one -- two rows for one residue would leave the file
    saying nothing about which color wins.
    """
    stray = [site for site in REGION_D_SITES if site not in SITES["D"]]
    if stray:
        raise SystemExit(f"HA1 sites {stray} are not in antigenic region D")
    absent = [
        site for site in REGION_D_SITES if residue_number("HA1", site) not in modeled
    ]
    if absent:
        raise SystemExit(f"{PDB_ID} does not model HA1 sites {absent}")
    claimed = sorted(
        {site for protein, site, _, _ in mutations if protein == "HA1"}
        & set(REGION_D_SITES)
    )
    if claimed:
        raise SystemExit(
            f"HA1 sites {claimed} are painted by the mutation list already"
        )
    return [
        [
            POLYMER_CHAIN,
            residue_number("HA1", site),
            REGION_D_COLOR,
            f"{site}_HA1",
            "",
            (
                f"HA1 site {site}; antigenic region D of {CITATION}, in the "
                "220-loop of the receptor-binding site"
            ),
        ]
        for site in REGION_D_SITES
    ]


def write_csv(name, rows, note):
    """Write one CSV beside this script and say what went into it."""
    out_path = pathlib.Path(__file__).parent / name
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)
    print(f"wrote {out_path} ({len(rows)} rows: {note})")


def main():
    check_sites()
    numbering = load_numbering_map()
    model = load_structure()

    polymer, glycans, receptor = partition_residues(model)
    check_frame(polymer)
    modeled = {number for _, number, _ in polymer}
    check_mutations("PERTH_TO_SUBCLADE_K", PERTH_TO_SUBCLADE_K, modeled)
    check_mutations("DC_2023_TO_DARWIN_2025", DC_2023_TO_DARWIN_2025, modeled)

    # Each list is ordered so the file sorts by chain then residue.
    protein = polymer_rows(polymer, numbering)
    host = glycan_rows(glycans)
    lstc = receptor_rows(receptor)

    write_csv(
        "antigenic-regions-w-glycans.csv",
        protein + host + lstc,
        f"{len(protein)} polymer, {len(host)} host glycan, {len(lstc)} receptor",
    )
    write_csv(
        "antigenic-regions.csv",
        protein + lstc,
        f"{len(protein)} polymer, {len(lstc)} receptor; host glycans left out so "
        "the view's `glycans: hide` can take them away",
    )
    write_csv(
        "perth-2009-to-subclade-k.csv",
        mutation_rows(PERTH_TO_SUBCLADE_K, PERTH, SUBCLADE_K) + lstc,
        f"{len(PERTH_TO_SUBCLADE_K)} differing sites, {len(lstc)} receptor",
    )
    vaccine = mutation_rows(DC_2023_TO_DARWIN_2025, DC_2023, DARWIN_2025)
    write_csv(
        "2025-26-to-2026-27-vaccine.csv",
        vaccine + lstc,
        f"{len(DC_2023_TO_DARWIN_2025)} differing sites, {len(lstc)} receptor",
    )
    region_d = region_d_rows(DC_2023_TO_DARWIN_2025, modeled)
    write_csv(
        "subclade-k-with-region-d-mutations.csv",
        # Re-sorted rather than concatenated: 222 and 223 fall between the last
        # HA1 site of the mutation list and its HA2 one.
        sorted(vaccine + region_d, key=lambda row: row[1]) + lstc,
        f"{len(vaccine)} differing sites, {len(region_d)} region D, "
        f"{len(lstc)} receptor",
    )


if __name__ == "__main__":
    main()
