"""Regenerate the coloring CSVs of this repository's H1 structure page.

The CSVs are committed and no rule runs this script -- ``spec.yaml`` reads the
committed files. Run it by hand from the repository root when the inputs below
change::

    python data/prot-struct-viz_config/H1/make_coloring_csv.py

It writes one CSV per view of ``spec.yaml``, each named after the view that
reads it -- so a name here that no view claims is a name that has gone stale:

* ``antigenic-regions.csv`` -- every modeled residue of all three protomers, HA1
  colored by antigenic region. The host N-glycans get no row, so the view's
  ``glycans: hide`` can take them away;
* ``d-3-1-1-mutations.csv`` -- the sites that differ between the two subclades,
  in red;
* ``d-3-1-1-with-sites-155-157.csv`` -- that list again, with every site of
  `HIGHLIGHT_SITES` added in one color.

It exists because the ~1500 rows of the first are derived rather than typed, and
a table that large is only auditable if its derivation ships with it. Everything
it needs is read rather than hard-coded:

* this repository's H1N1 site-numbering map, which is what turns a residue number
  into an HA1 or HA2 site number;
* this repository's viral library, whose ``protein_sequence_HA_ectodomain`` is
  what the differing sites are computed from -- there is no transcribed mutation
  list here, only the two ``derived_haplotype`` values that pick the strains;
* ``config.yml``, for the residues the library's ectodomains leave off the start
  of HA1;
* PDB 9GSP, fetched from RCSB, so that only *modeled* residues get a row.

The antigenic-region definitions are the exception: they are transcribed from the
table cited below, and `check_sites` is what guards the transcription.

The numbering frame is the load-bearing assumption throughout, so it is asserted
rather than trusted: if a revised map or a different entry ever shifts it, this
script fails instead of writing a plausible-looking but wrong CSV. For the same
reason it reads the captions back and holds what they spell out -- the sites of
each antigenic region, and each substitution named in prose -- to what it draws.
"""

import csv
import pathlib
import re
import urllib.request

import gemmi
import yaml

#: This repository, located from this script rather than from the working
#: directory, so the run path above is the only one a reader needs.
REPO = pathlib.Path(__file__).parents[3]

#: Author numbering in 9GSP is a single chain running across HA1 and HA2, which
#: is the same frame as this map's ``sequential_site``.
NUMBERING_MAP = (
    REPO / "data/nextstrain-prot-titers-tree_data/H1N1_site_numbering_map.tsv"
)

#: The library the differing sites are computed from, and the subtype whose rows
#: of it are read.
LIBRARY = REPO / "data/viral_libraries/flu-seqneut-2026-barcode-to-strain-actual.csv"
SUBTYPE = "H1N1"

#: Read for one key: the residues prepended to the library's ectodomain
#: sequences to bring them back to HA1 site 1. Its length is the offset between
#: an ectodomain position and a sequential site, and this is where the project
#: already records it.
CONFIG = REPO / "config.yml"
PREFIX_KEY = "nextstrain-prot-titers-tree_prefix_alignment"

PDB_ID = "9GSP"
STRUCTURE_URL = f"https://files.rcsb.org/download/{PDB_ID}.cif"

#: The three protomers. 9GSP's deposited coordinates are already the trimer, so
#: every protomer needs its own rows -- there is no symmetry expansion to
#: replicate one chain's annotation onto the other two.
POLYMER_CHAINS = ("A", "B", "C")

#: ``derived_haplotype`` values of the two strains compared. This is the one
#: experimental choice the script makes: each subclade is represented by the one
#: strain the library assigns that subclade as its whole haplotype.
FROM_HAPLOTYPE = "D.3.1"
TO_HAPLOTYPE = "D.3.1.1"

#: HA1 sites drawn on top of the mutation list by the view that adds them. These
#: are named rather than derived, because what picks them is which sites are of
#: interest rather than any property the library could be asked for. What is
#: derived is everything then said about them: `background_substitutions` reads
#: the substitutions the library carries at each, and fails if a site has none --
#: which is what catches a number typed wrong here.
HIGHLIGHT_SITES = (155, 157)

#: Antigenic sites of H1 HA1, from Table 2 ("Amino acid sequence of antigenic
#: sites for historic H1N1 and H1N1pdm natural isolate viruses") of Wilson et
#: al., Virology 2015;485:252-62, PMC5737639. That table's footnote: "Antigenic
#: sites are based on those determined for A/PR/8/34 by Caton et al. (1982)."
#: The residue numbers themselves are those of the table's reference row,
#: A/California/07/2009, so this is Caton's set carried into the numbering of
#: the pandemic 2009 lineage rather than into A/PR/8/34's -- which is why it
#: transfers to 9GSP without an alignment step. Every residue the table shows is
#: included: unlike the H3 table it encloses none of them in parentheses.
SITES = {
    "Sa": [124, 125, 153, 154, 155, 156, 157, 159, 160, 161, 162, 163, 164],
    "Sb": [184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195],
    "Ca1": [166, 167, 168, 169, 170, 203, 204, 205, 235, 236, 237],
    "Ca2": [137, 138, 139, 140, 141, 142, 221, 222],
    "Cb": [70, 71, 72, 73, 74, 75],
}

#: The H3 page's colorblind-safe qualitative palette, reused in a declared order.
#: The order is deliberately not an attempt to match each H1 site to the H3 site
#: it sits nearest: two of them would then share a color, and the reader would be
#: invited to read an equivalence off the two pages that the two sets of
#: definitions do not actually assert.
# Paul Tol's "muted" scheme, cool subset. Every one of these is cool, so on the
# structure warm means "changed" -- the red of a mutated site -- and nothing else
# is close to MUTATED_COLOR.
SITE_COLORS = {
    "Sa": "#332288",
    "Sb": "#88CCEE",
    "Ca1": "#117733",
    "Ca2": "#AA4499",
    "Cb": "#44AA99",
}

#: The two grays are not decoration: they draw the HA1/HA2 boundary, which is
#: what the _HA1 and _HA2 label suffixes exist to tell you about. 9GSP is an
#: uncleaved HA0, so that boundary is a position in one polypeptide rather than
#: a break between two -- and the loop across it is disordered.
HA1_COLOR = "#e8e8e8"
HA2_COLOR = "#bdbdbd"

#: Sites that differ between the two subclades. Warm, where every `SITE_COLORS`
#: entry is cool: the two mutation views drop the antigenic-region coloring, and
#: the one of them that puts a region color back has to hold it apart from this.
MUTATED_COLOR = "#e41a1c"

#: Antigenic region every one of `HIGHLIGHT_SITES` has to be in. The caption names
#: that region in prose, so this is what holds the two together -- and it is where
#: `HIGHLIGHT_SITES_COLOR` comes from.
HIGHLIGHT_REGION = "Sa"

#: One color for every one of `HIGHLIGHT_SITES`, the indigo the antigenic-regions
#: view gives `HIGHLIGHT_REGION`. `named_site_rows` holds every highlighted site
#: to that region, so the shared color states something true rather than being
#: only a choice of ink.
HIGHLIGHT_SITES_COLOR = SITE_COLORS[HIGHLIGHT_REGION]

CITATION = "Wilson et al. 2015 Virology 485:252-62"

#: Caption of the antigenic-regions view. It spells out the sites of each
#: region, so `check_caption_sites` holds it to `SITES`.
ANTIGENIC_CAPTION = "antigenic-regions.md"

#: Caption of the view that adds `HIGHLIGHT_SITES`. Which of the substitutions at
#: those sites it singles out is a scientific choice, so `check_highlight_caption`
#: does not require it to name them all -- only that whatever it does name is a
#: substitution the library actually carries there.
HIGHLIGHT_SITES_CAPTION = "d-3-1-1-with-sites-155-157.md"

#: Captions of the mutation views. They name each substitution between the two
#: subclades in prose, so `check_mutation_captions` holds them to what the library
#: says.
MUTATION_CAPTIONS = ("d-3-1-1-mutations.md", HIGHLIGHT_SITES_CAPTION)

#: CSV of the view that adds `HIGHLIGHT_SITES`. Unlike the other names it is a
#: constant, because it carries those sites' numbers and `named_site_rows` checks
#: the name still agrees with what is drawn into it.
HIGHLIGHT_SITES_CSV = "d-3-1-1-with-sites-155-157.csv"

#: Residue number of HA2 site 1 in 9GSP's author numbering. HA1 sites are author
#: numbers already, so HA1 needs no offset.
HA2_OFFSET = 327

#: Residues whose identity fixes the numbering frame, quoted in the same H1
#: numbering as the site definitions. The nine cysteines form HA1's disulfides
#: and the four aromatics line the receptor-binding site; the three HA2 entries
#: are there to pin the HA1/HA2 offset, and are quoted in author numbering like
#: the rest (HA2 14, 21 and 26).
FRAME_CHECKS = {
    4: "CYS",
    42: "CYS",
    55: "CYS",
    67: "CYS",
    90: "CYS",
    98: "TYR",
    136: "CYS",
    150: "TRP",
    180: "HIS",
    192: "TYR",
    275: "CYS",
    279: "CYS",
    303: "CYS",
    341: "TRP",
    348: "TRP",
    353: "HIS",
}

#: No view of this page draws labels into the scene, so the three columns that
#: style a drawn label -- ``show_label``, ``label_color``, ``label_size`` -- are
#: left out rather than written empty. ``label`` stays: it is the mouseover
#: tooltip, which is how a reader identifies a residue here.
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


def load_prefix_length():
    """Return how many HA1 residues the library's ectodomain sequences omit."""
    with CONFIG.open() as handle:
        return len(yaml.safe_load(handle)[PREFIX_KEY][SUBTYPE])


def load_library():
    """Return ``{derived_haplotype: (strain, accession, sequence)}``.

    Every strain of `SUBTYPE`, not only the ones compared: the haplotypes on the
    `TO_HAPLOTYPE` background are what `background_substitutions` reads.

    One strain per haplotype and one sequence per strain, both required rather
    than assumed: the library holds several barcodes of each strain, and a
    haplotype naming two strains would leave a comparison ambiguous.
    """
    wanted = (FROM_HAPLOTYPE, TO_HAPLOTYPE)
    found = {}
    with LIBRARY.open() as handle:
        for row in csv.DictReader(handle):
            haplotype = row["derived_haplotype"]
            if row["subtype"] != SUBTYPE:
                continue
            entry = (
                row["strain"],
                row["genbank_accession"],
                row["protein_sequence_HA_ectodomain"],
            )
            if haplotype in found and found[haplotype] != entry:
                raise SystemExit(
                    f"{LIBRARY.name}: derived_haplotype {haplotype!r} does not "
                    f"name one strain and sequence; found {found[haplotype][0]} "
                    f"and {entry[0]}"
                )
            found[haplotype] = entry
    absent = [haplotype for haplotype in wanted if haplotype not in found]
    if absent:
        raise SystemExit(
            f"{LIBRARY.name}: no {SUBTYPE} strain has derived_haplotype {absent}"
        )
    return found


def load_structure():
    """Return the 9GSP model, fetched from RCSB."""
    with urllib.request.urlopen(STRUCTURE_URL) as response:
        text = response.read().decode()
    structure = gemmi.make_structure_from_block(
        gemmi.cif.read_string(text).sole_block()
    )
    structure.setup_entities()
    return structure[0]


def polymer_residues(model):
    """Return the modeled polymer residues of every protomer.

    gemmi splits waters and heteroatoms into a second chain reusing the author
    ID, so the N-acetylglucosamines that hang off each protomer appear under that
    protomer's chain. Classify by component rather than by chain: amino acids on
    a protomer chain are the protein, and anything else non-water is a host
    N-glycan, which gets no row -- the views' ``glycans: hide`` is what takes
    those away.

    Returns a ``[(chain, residue_number, component_name), ...]`` list.
    """
    polymer = []
    for chain in model:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            if info is not None and info.is_water():
                continue
            if (
                chain.name in POLYMER_CHAINS
                and info is not None
                and info.is_amino_acid()
            ):
                polymer.append((chain.name, residue.seqid.num, residue.name))
    absent = set(POLYMER_CHAINS) - {chain for chain, _, _ in polymer}
    if absent:
        raise SystemExit(f"{PDB_ID} has no polymer on chain(s) {sorted(absent)}")
    return polymer


def check_frame(polymer):
    """Fail loudly if author numbering is not H1 numbering after all."""
    for want_chain in POLYMER_CHAINS:
        seen = {
            number: component
            for chain, number, component in polymer
            if chain == want_chain
        }
        wrong = {
            num: (expected, seen.get(num))
            for num, expected in FRAME_CHECKS.items()
            if seen.get(num) != expected
        }
        if wrong:
            raise SystemExit(
                f"{PDB_ID} chain {want_chain} is not in H1 numbering: "
                + "; ".join(
                    f"expected {exp} at {num}, found {got}"
                    for num, (exp, got) in sorted(wrong.items())
                )
            )


def check_ha2_offset(numbering):
    """Fail loudly if the numbering map's own HA2 offset is not `HA2_OFFSET`.

    The map derives it from its own rows; `HA2_OFFSET` is what `residue_number`
    uses to turn an HA2 site into an author number. They are the same number for
    the same reason -- HA1 is 327 residues -- but nothing forces that, so it is
    checked once here rather than trusted in two places.
    """
    offset = max(
        seq for seq, (protein, _) in numbering.items() if protein == "HA2"
    ) - max(site for protein, site in numbering.values() if protein == "HA2")
    if offset != HA2_OFFSET:
        raise SystemExit(
            f"numbering map puts HA2 site 1 at residue {offset + 1}, but "
            f"{HA2_OFFSET + 1} is what this script assumes"
        )


def check_sites():
    """Fail loudly if the transcribed site definitions have drifted."""
    flat = [site for residues in SITES.values() for site in residues]
    if len(flat) != len(set(flat)):
        raise SystemExit("antigenic sites overlap; they are disjoint in Table 2")
    if len(flat) != 50:
        raise SystemExit(f"expected 50 antigenic-site residues, got {len(flat)}")


def check_caption_sites():
    """Fail loudly if `ANTIGENIC_CAPTION`'s site lists have drifted from `SITES`.

    The caption spells the sites out for a reader, which puts the same 50 numbers
    in two places. Rather than trust them to stay in step, read them back and
    compare: a revised Table 2 transcription then fails here instead of leaving
    the page's key quietly describing the previous one.
    """
    text = (pathlib.Path(__file__).parent / ANTIGENIC_CAPTION).read_text()
    for region, sites in SITES.items():
        matches = re.findall(rf"^- .*\bsite {region}\b[^:]*:(.*)$", text, re.MULTILINE)
        if len(matches) != 1:
            raise SystemExit(
                f"{ANTIGENIC_CAPTION}: expected one bullet listing site "
                f"{region}, found {len(matches)}"
            )
        listed = set()
        for token in matches[0].split(","):
            bounds = [int(bound) for bound in token.strip().split("-")]
            listed |= set(range(bounds[0], bounds[-1] + 1))
        if listed != set(sites):
            raise SystemExit(
                f"{ANTIGENIC_CAPTION}: site {region} lists {sorted(listed)}, "
                f"but SITES has {sorted(sites)}"
            )


def check_mutation_captions(mutations):
    """Fail loudly if a mutation caption does not name every site it draws.

    The captions group the substitutions by antigenic region in prose, which is
    the only place they are typed now that the list itself is derived. Reading
    them back is what keeps a revised library from leaving them describing the
    previous one.
    """
    named = [f"{old}{site}{new}" for _, site, old, new in mutations]
    count = f"differ at {len(named)} sites"
    for caption in MUTATION_CAPTIONS:
        text = (pathlib.Path(__file__).parent / caption).read_text()
        absent = [mutation for mutation in named if mutation not in text]
        if absent:
            raise SystemExit(f"{caption} does not name {absent}")
        if count not in text:
            raise SystemExit(f"{caption} does not say the subclades {count}")


def check_highlight_caption(substitutions):
    """Fail loudly if `HIGHLIGHT_SITES_CAPTION` names a substitution nothing has.

    Which substitutions at the highlighted sites to single out is a scientific
    choice, so the caption is not held to naming all of them -- the CSV's tooltips
    carry the full derived list. What it is held to is that every substitution it
    does name at one of those sites is one the library carries there, so a typo or
    a revised library fails here rather than leaving prose the data contradicts.
    """
    text = (pathlib.Path(__file__).parent / HIGHLIGHT_SITES_CAPTION).read_text()
    wrong = [
        substitution
        for substitution, site in re.findall(r"\b([A-Z](\d+)[A-Z])\b", text)
        if int(site) in substitutions and substitution not in substitutions[int(site)]
    ]
    if wrong:
        raise SystemExit(
            f"{HIGHLIGHT_SITES_CAPTION} names {wrong}, which no {TO_HAPLOTYPE}"
            f"-background strain of {LIBRARY.name} carries"
        )


def site_of(residue_number_):
    """Return the antigenic site containing this HA1 site, or ``None``."""
    for site, residues in SITES.items():
        if residue_number_ in residues:
            return site
    return None


def residue_number(protein, site):
    """Return the 9GSP author residue number of an HA1 or HA2 site."""
    if protein == "HA1":
        return site
    if protein == "HA2":
        return site + HA2_OFFSET
    raise SystemExit(f"unknown protein {protein!r}; expected 'HA1' or 'HA2'")


def differing_sites(from_sequence, to_sequence, numbering, prefix_length):
    """Return ``[(protein, site, from, to), ...]`` for two ectodomain sequences.

    The library's ectodomains leave off the first `prefix_length` residues of HA1
    and hold no indels, so ectodomain position ``p`` is sequential site
    ``p + prefix_length`` and the two sequences compare position by position with
    no alignment. Both of those are asserted rather than assumed, since either
    one being false would shift every site silently.
    """
    if len(from_sequence) != len(to_sequence):
        raise SystemExit(
            f"the two ectodomains are {len(from_sequence)} and "
            f"{len(to_sequence)} residues, so they need an alignment"
        )
    if len(from_sequence) + prefix_length != len(numbering):
        raise SystemExit(
            f"a {len(from_sequence)}-residue ectodomain plus {prefix_length} "
            f"prepended residues is not the numbering map's {len(numbering)} sites"
        )
    return [
        (*numbering[position + prefix_length], old, new)
        for position, (old, new) in enumerate(zip(from_sequence, to_sequence), start=1)
        if old != new
    ]


def check_mutations(mutations, modeled):
    """Fail loudly if a differing site is unmodeled or named twice.

    A site 9GSP does not model would be reported by ``on_mismatch``, but only as
    one line among the rest of the build's output, so it fails here instead. Two
    rows for one residue would leave the CSV saying nothing about which wins.
    """
    seen = {(protein, site) for protein, site, _, _ in mutations}
    if len(seen) != len(mutations):
        raise SystemExit("the same site differs twice; the numbering map is not 1:1")
    absent = [
        f"{site}_{protein}"
        for protein, site, _, _ in mutations
        if residue_number(protein, site) not in modeled
    ]
    if absent:
        raise SystemExit(f"{PDB_ID} does not model differing sites {absent}")


def background_substitutions(library, numbering, prefix_length):
    """Return ``{site: [substitution, ...]}`` at each of `HIGHLIGHT_SITES`.

    Every library strain whose haplotype is `TO_HAPLOTYPE` or an extension of it
    is compared to `TO_HAPLOTYPE` itself, so what is reported at a highlighted
    site is read off the sequences rather than off the haplotype names. A site
    nothing mutates is a site typed wrong in `HIGHLIGHT_SITES`, and fails here
    rather than reaching the page as an unexplained patch of color.
    """
    found = {site: set() for site in HIGHLIGHT_SITES}
    for haplotype, (_, _, sequence) in library.items():
        if haplotype != TO_HAPLOTYPE and not haplotype.startswith(f"{TO_HAPLOTYPE}:"):
            continue
        for protein, site, old, new in differing_sites(
            library[TO_HAPLOTYPE][2], sequence, numbering, prefix_length
        ):
            if protein == "HA1" and site in found:
                found[site].add(f"{old}{site}{new}")
    absent = sorted(site for site, seen in found.items() if not seen)
    if absent:
        raise SystemExit(
            f"no {TO_HAPLOTYPE}-background strain of {LIBRARY.name} mutates HA1 "
            f"site(s) {absent}"
        )
    return {site: sorted(seen) for site, seen in found.items()}


def polymer_rows(polymer, numbering):
    """One row per modeled polymer residue, colored by antigenic site."""
    rows = []
    for chain, number, _component in polymer:
        if number not in numbering:
            raise SystemExit(
                f"{PDB_ID} chain {chain} models residue {number}, which is past "
                "the end of the numbering map"
            )
        protein, site = numbering[number]
        region = site_of(site) if protein == "HA1" else None
        # One label for every residue, antigenic site or not: the site number
        # alone is ambiguous across HA1 and HA2, and this is the mouseover text.
        label = f"{site}_{protein}"
        if region is not None:
            rows.append(
                [
                    chain,
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
                    chain,
                    number,
                    HA1_COLOR if protein == "HA1" else HA2_COLOR,
                    label,
                    "",
                    f"{protein} site {site}; not in a defined antigenic site",
                ]
            )
    return rows


def mutation_rows(mutations, from_name, to_name):
    """Red rows for the differing sites, on every protomer.

    The antigenic-site coloring is gone from these views: every residue without
    a row here takes the view's ``default_color``, so the only thing painted is
    what changed. The tooltip carries the substitution itself, which is the only
    way to read the view -- no label is drawn into the scene.
    """
    return [
        [
            chain,
            residue_number(protein, site),
            MUTATED_COLOR,
            f"{old}{site}{new}_{protein}",
            "",
            (
                f"{protein} site {site} differs: {old} in {from_name}, "
                f"{new} in {to_name}"
            ),
        ]
        for chain in POLYMER_CHAINS
        for protein, site, old, new in mutations
    ]


def site_rows(site, color, note, mutations, modeled):
    """A row per protomer for one HA1 site drawn on top of the mutation list.

    These are written alongside `mutation_rows` into a single CSV, so the two
    things that could make that CSV lie are checked here: that 9GSP models the
    site, and that the mutation list does not already claim it -- two rows for
    one residue would leave the file saying nothing about which color wins.
    """
    number = residue_number("HA1", site)
    if number not in modeled:
        raise SystemExit(f"{PDB_ID} does not model {site}_HA1")
    if ("HA1", site) in {(p, s) for p, s, _, _ in mutations}:
        raise SystemExit(f"{site}_HA1 is painted by the mutation list already")
    return [[chain, number, color, f"{site}_HA1", "", note] for chain in POLYMER_CHAINS]


def named_site_rows(substitutions, mutations, modeled):
    """Rows for every one of `HIGHLIGHT_SITES`, all in `HIGHLIGHT_SITES_COLOR`.

    They share one color, and that color is `HIGHLIGHT_REGION`'s, so every site
    drawn here has to be in that region: a site outside it would be painted in
    the ink the antigenic-regions view reserves for one that is.
    """
    outside = [site for site in HIGHLIGHT_SITES if site_of(site) != HIGHLIGHT_REGION]
    if outside:
        raise SystemExit(
            f"{HIGHLIGHT_SITES_CSV} would paint site(s) {outside} in the color of "
            f"antigenic region {HIGHLIGHT_REGION}, which they are not in"
        )
    absent = [site for site in HIGHLIGHT_SITES if str(site) not in HIGHLIGHT_SITES_CSV]
    if absent:
        raise SystemExit(f"{HIGHLIGHT_SITES_CSV} is not named for site(s) {absent}")
    rows = []
    for site in sorted(HIGHLIGHT_SITES):
        rows += site_rows(
            site,
            HIGHLIGHT_SITES_COLOR,
            (
                f"HA1 site {site}; antigenic region {HIGHLIGHT_REGION} of "
                f"{CITATION}; the site of {', '.join(substitutions[site])} on the "
                f"{TO_HAPLOTYPE} background"
            ),
            mutations,
            modeled,
        )
    return rows


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
    check_caption_sites()
    numbering = load_numbering_map()
    check_ha2_offset(numbering)
    prefix_length = load_prefix_length()
    library = load_library()
    model = load_structure()

    polymer = polymer_residues(model)
    check_frame(polymer)
    modeled = {number for _, number, _ in polymer}

    from_strain, from_accession, from_sequence = library[FROM_HAPLOTYPE]
    to_strain, to_accession, to_sequence = library[TO_HAPLOTYPE]
    mutations = differing_sites(from_sequence, to_sequence, numbering, prefix_length)
    check_mutations(mutations, modeled)
    check_mutation_captions(mutations)
    substitutions = background_substitutions(library, numbering, prefix_length)
    check_highlight_caption(substitutions)
    print(
        f"{FROM_HAPLOTYPE} ({from_strain}) to {TO_HAPLOTYPE} ({to_strain}) differs at "
        + ", ".join(
            f"{old}{site}{new}_{protein}" for protein, site, old, new in mutations
        )
    )
    print(
        f"{TO_HAPLOTYPE}-background strains of the library mutate the highlighted "
        "sites with "
        + "; ".join(", ".join(substitutions[site]) for site in sorted(substitutions))
    )

    # Every CSV is sorted by chain then residue, so a reader can find a residue
    # in it the same way in all of them.
    def by_position(rows):
        return sorted(rows, key=lambda row: (row[0], row[1]))

    write_csv(
        "antigenic-regions.csv",
        by_position(polymer_rows(polymer, numbering)),
        f"{len(polymer)} modeled residues across {len(POLYMER_CHAINS)} protomers; "
        "host glycans left out so the view's `glycans: hide` can take them away",
    )
    subclade = mutation_rows(
        mutations,
        f"subclade {FROM_HAPLOTYPE} ({from_strain}, {from_accession})",
        f"subclade {TO_HAPLOTYPE} ({to_strain}, {to_accession})",
    )
    write_csv(
        "d-3-1-1-mutations.csv",
        by_position(subclade),
        f"{len(mutations)} differing sites on each of "
        f"{len(POLYMER_CHAINS)} protomers",
    )
    write_csv(
        HIGHLIGHT_SITES_CSV,
        by_position(subclade + named_site_rows(substitutions, mutations, modeled)),
        f"{len(mutations)} differing sites and HA1 "
        + ", ".join(str(site) for site in sorted(HIGHLIGHT_SITES))
        + f" on each of {len(POLYMER_CHAINS)} protomers",
    )


if __name__ == "__main__":
    main()
