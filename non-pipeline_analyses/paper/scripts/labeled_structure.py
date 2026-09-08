"""Compose a structure figure: the rendered HA surface with its colored sites labeled.

`rule molstar_render` rendered the Mol* page to an image and reported where on it each
labeled residue sits. This lays that image into an SVG and draws the labels over it as
vector text on leader lines, which is both crisp in print and impossible for the surface
to bury -- Mol*'s own on-structure labels are depth-tested and sink into a surface this
bumpy.

The label text and color are the coloring CSV's, so nothing about what a site is called
or how it is colored is restated here. What is here is only where each label sits, which
is a matter of the silhouette a particular camera gives, and so is per figure: the
`figure` wildcard says which.
"""

import base64
import colorsys
import csv
import struct
import sys

import numpy
from PIL import Image

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# Where each site's label goes: which side of the structure it sits on, the x of the text
# edge nearest the structure, and the y of the text's centre -- all in the rendered
# image's pixels, so x may be negative or past the image where a label reaches into the
# margin. Read down each side and the values follow the silhouette: a wide head takes its
# labels close in, a narrow stalk lets them come well inside. The y values are ordered to
# match the anchors they point at, which is what keeps the leader lines from crossing.
LABEL_PLACEMENTS = {
    "H3_subclade_K_mutations_structure": {
        "158": ("left", 310, 280),
        "160": ("left", 350, 440),
        "135": ("left", 125, 710),
        "378": ("left", 440, 2025),
        "189": ("right", 1725, 260),
        "145": ("right", 1725, 490),
        "144": ("right", 1725, 615),
        "222": ("right", 1725, 760),
        "223": ("right", 1725, 895),
        "173": ("right", 1725, 1075),
    },
    "H1_D_3_1_1_mutations_structure": {
        "155": ("left", 380, 250),
        "139": ("left", 180, 900),
        "157": ("right", 1600, 200),
        "113": ("right", 1750, 700),
        "283": ("right", 1760, 1300),
        "302": ("right", 1700, 1620),
    },
}

# A leader points at the centre of the colored patch the render draws for its site,
# found in the image rather than taken from the projected position `rule molstar_render`
# reports. That position is where the residue's nearest atom lands, and the surface
# bulges over it, so it can sit beside the patch rather than in it -- which is what put
# H3's 144 leader on bare grey. The projected positions are the seeds: a colored pixel is
# claimed by whichever seed is nearest, and that is what stops three adjacent sites in
# one unbroken red region from all pointing at the same place. A site is matched on hue
# alone, because the surface shades its color over most of the brightness range.
PATCH_SEARCH_RADIUS = 85
PATCH_MIN_SATURATION = 0.3
PATCH_HUE_TOLERANCE = 0.055

# Each site is drawn once per protomer, and the render reports where every copy projects
# to. The copy labeled is the nearest one the render actually shows: nearest alone is not
# enough, because a copy can face the camera and still have its patch hidden behind what
# is in front of it, which is where H1's 283 sits on the front protomer. A patch has to
# reach this many pixels to count as shown -- the ones really on screen run to thousands,
# so this only rules out a sliver around an edge.
MIN_PATCH_PIXELS = 400

# The image the placements above were tuned against. They are pixel positions over one
# particular camera at one particular size, so a differently sized render would leave
# every label pointing somewhere plausible but wrong; this makes that a failure instead.
TUNED_FOR_IMAGE = (1996, 2996)

# The figure being composed, which is all this script needs the wildcard for.
FIGURE = snakemake.wildcards.figure
if FIGURE not in LABEL_PLACEMENTS:
    raise ValueError(f"no label placements for figure {FIGURE!r}")
PLACEMENTS = LABEL_PLACEMENTS[FIGURE]

# Blank space kept either side of the image for the labels that reach outside it.
MARGIN = 210

# Type and leader lines, in image pixels. The image is about 2000 x 3000.
FONT_SIZE = 82
FONT_FAMILY = "Helvetica, Arial, sans-serif"
LEADER_WIDTH = 5
ANCHOR_RADIUS = 13

# Gap between a leader line and the text it points to, so the line does not touch a glyph.
TEXT_GAP = 17

# Baseline offset that centres a line of text on its given y. Set here rather than with
# `dominant-baseline`, which not every SVG renderer honours.
BASELINE_SHIFT = 0.36

# The CSV writes a site's compartment into the label after an underscore (`T135K_HA1`)
# because that label is also a Mol* tooltip, which has no surrounding text to say which
# numbering it is in. A figure does: HA1 numbering is what all but one of these sites are
# in and what the caption will say, so it goes unwritten, and HA2 is spelled out as the
# exception it is.
COMPARTMENT_SEPARATOR = "_"
IMPLIED_COMPARTMENT = "HA1"


def png_size(png_path):
    """The `(width, height)` in the PNG's header."""
    with open(png_path, "rb") as f:
        header = f.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError(f"{png_path} is not a PNG")
    return struct.unpack(">II", header[16:24])


def read_projected_copies(path):
    """The projected copies in `path`, nearest the camera first within each site."""
    with open(path, newline="") as f:
        copies = [
            {**row, "x": float(row["x"]), "y": float(row["y"])}
            for row in csv.DictReader(f)
        ]
    return sorted(copies, key=lambda copy: float(copy["camera_distance"]))


def read_csv(path, key_columns):
    """`path` as a dict keyed on `key_columns`, whose combination must be unique."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    keyed = {tuple(row[column] for column in key_columns): row for row in rows}
    if len(keyed) != len(rows):
        raise ValueError(f"{path}: {key_columns} is not unique")
    return keyed


def hue_of(hex_color):
    """The hue, in turns, of a `#rrggbb` color."""
    channels = hex_color.lstrip("#")
    if len(channels) != 6:
        raise ValueError(f"{hex_color!r} is not a #rrggbb color")
    return colorsys.rgb_to_hsv(
        *(int(channels[i : i + 2], 16) / 255 for i in (0, 2, 4))
    )[0]


def patch_centres(structure_png, copies, coloring):
    """For each site, the centre of the patch of the nearest copy the render shows.

    `copies` is every projected copy of every site, nearest the camera first. Colored
    pixels are claimed by whichever copy's projected position is nearest, so a copy on a
    protomer facing away competes for nothing, and each site is then taken from the first
    of its copies whose claim is big enough to be a patch rather than a sliver.
    """
    rgb = numpy.asarray(Image.open(structure_png).convert("RGB"), dtype=numpy.float32)
    rgb /= 255.0
    red, green, blue = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    high, low = rgb.max(axis=2), rgb.min(axis=2)
    spread = numpy.maximum(high - low, 1e-6)
    saturation = numpy.where(high > 0, (high - low) / numpy.maximum(high, 1e-6), 0.0)
    hue = (
        numpy.select(
            [high == red, high == green],
            [((green - blue) / spread) % 6, (blue - red) / spread + 2],
            (red - green) / spread + 4,
        )
        / 6.0
    )
    colored = saturation >= PATCH_MIN_SATURATION

    rows, columns = numpy.mgrid[0 : rgb.shape[0], 0 : rgb.shape[1]].astype(
        numpy.float32
    )
    distances = numpy.stack(
        [(columns - copy["x"]) ** 2 + (rows - copy["y"]) ** 2 for copy in copies]
    )
    claimed_by = distances.argmin(axis=0)
    within_reach = distances.min(axis=0) <= PATCH_SEARCH_RADIUS**2

    by_site = {}
    for index, copy in enumerate(copies):
        from_hue = numpy.abs(
            (hue - hue_of(coloring[copy["residue"]]["color"]) + 0.5) % 1.0 - 0.5
        )
        patch = (
            colored
            & (from_hue <= PATCH_HUE_TOLERANCE)
            & (claimed_by == index)
            & within_reach
        )
        by_site.setdefault(copy["residue"], []).append((copy, int(patch.sum()), patch))

    centres = {}
    for site, found in sorted(by_site.items()):
        shown = [entry for entry in found if entry[1] >= MIN_PATCH_PIXELS]
        sizes = ", ".join(
            f"{copy['chain']}/{copy['unit']} {size} px" for copy, size, _ in found
        )
        if not shown:
            raise ValueError(
                f"{structure_png} shows no copy of site {site} in "
                f"{coloring[site]['color']}: patches were {sizes}, against the "
                f"{MIN_PATCH_PIXELS} px that counts as shown"
            )
        copy, _, patch = shown[0]
        centres[site] = (float(columns[patch].mean()), float(rows[patch].mean()))
        moved = (
            (centres[site][0] - copy["x"]) ** 2 + (centres[site][1] - copy["y"]) ** 2
        ) ** 0.5
        print(
            f"  {site}: patches {sizes}; labeling {copy['chain']}/{copy['unit']} "
            f"centred ({centres[site][0]:.5g}, {centres[site][1]:.5g}), "
            f"{moved:.5g} px from where it projects"
        )
    return centres


def label_text(coloring_row):
    """What a site's label reads in the figure."""
    label = coloring_row["label"]
    site, separator, compartment = label.partition(COMPARTMENT_SEPARATOR)
    if not separator:
        raise ValueError(
            f"label {label!r} names no compartment after "
            f"{COMPARTMENT_SEPARATOR!r}, so there is no telling whether it is in "
            f"{IMPLIED_COMPARTMENT} numbering"
        )
    return site if compartment == IMPLIED_COMPARTMENT else f"{site} {compartment}"


def leader_and_text(site, anchor, coloring_row):
    """The SVG for one label: its leader line, the dot on the site, and the text."""
    side, edge_x, label_y = PLACEMENTS[site]
    if side not in ("left", "right"):
        raise ValueError(f"site {site} is placed {side!r}, not 'left' or 'right'")
    anchor_x, anchor_y = anchor
    gap = TEXT_GAP if side == "left" else -TEXT_GAP
    color = coloring_row["color"]
    text = label_text(coloring_row)
    print(
        f"  {site}: {text!r} {side} of the structure, text edge at "
        f"({edge_x}, {label_y}), leader to ({anchor_x:.5g}, {anchor_y:.5g})"
    )
    return [
        (
            f'<line x1="{fmt(anchor_x + MARGIN)}" y1="{fmt(anchor_y)}" '
            f'x2="{fmt(edge_x + MARGIN + gap)}" y2="{fmt(label_y)}" '
            f'stroke="{color}" stroke-width="{LEADER_WIDTH}"/>'
        ),
        (
            f'<circle cx="{fmt(anchor_x + MARGIN)}" cy="{fmt(anchor_y)}" '
            f'r="{ANCHOR_RADIUS}" fill="{color}"/>'
        ),
        (
            f'<text x="{fmt(edge_x + MARGIN)}" '
            f'y="{fmt(label_y + FONT_SIZE * BASELINE_SHIFT)}" '
            f'text-anchor="{"end" if side == "left" else "start"}" '
            f'font-family="{FONT_FAMILY}" font-size="{FONT_SIZE}" fill="{color}">'
            f"{text}</text>"
        ),
    ]


def fmt(value):
    """A number written for an SVG attribute, without a trailing `.0`."""
    return f"{value:.5g}"


# The anchors name the chain each site was labeled on, which is how the right row of a
# coloring CSV that records the same site once per protomer is found.
copies = read_projected_copies(snakemake.input.anchors_csv)
coloring_by_chain = read_csv(snakemake.input.coloring_csv, ["chain", "residue"])
coloring = {
    copy["residue"]: coloring_by_chain[(copy["chain"], copy["residue"])]
    for copy in copies
}

unplaced = sorted(set(coloring) - set(PLACEMENTS))
unanchored = sorted(set(PLACEMENTS) - set(coloring))
if unplaced or unanchored:
    raise ValueError(
        f"{snakemake.input.anchors_csv} and this figure's placements disagree: "
        f"no placement for {unplaced}, no anchor for {unanchored}"
    )

width, height = png_size(snakemake.input.structure_png)
if (width, height) != TUNED_FOR_IMAGE:
    raise ValueError(
        f"{snakemake.input.structure_png} is {width} x {height}, but this figure's "
        f"label placements were tuned against {TUNED_FOR_IMAGE[0]} x "
        f"{TUNED_FOR_IMAGE[1]}; re-tune them against the new render"
    )
print(f"{snakemake.input.structure_png} is {width} x {height}")

anchors = patch_centres(snakemake.input.structure_png, copies, coloring)
png = base64.b64encode(open(snakemake.input.structure_png, "rb").read()).decode()

labels = []
for site in sorted(anchors, key=lambda s: PLACEMENTS[s][2]):
    labels += leader_and_text(site, anchors[site], coloring[site])

svg = "\n".join(
    [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{width + 2 * MARGIN}" height="{height}" '
            f'viewBox="0 0 {width + 2 * MARGIN} {height}">'
        ),
        f'<rect width="{width + 2 * MARGIN}" height="{height}" fill="white"/>',
        (
            f'<image x="{MARGIN}" y="0" width="{width}" height="{height}" '
            f'xlink:href="data:image/png;base64,{png}"/>'
        ),
        *labels,
        "</svg>",
        "",
    ]
)
with open(snakemake.output.figure_svg, "w", encoding="utf-8") as f:
    f.write(svg)

print(
    f"wrote {snakemake.output.figure_svg}: {width + 2 * MARGIN} x {height} "
    f"with {len(anchors)} labels ({len(svg) / 1e6:.3g} MB)"
)
