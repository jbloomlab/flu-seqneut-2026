"""Render a structure page to a high-resolution image and locate its labeled sites.

The page is a Mol* viewer, so rendering it means running it: a headless browser loads a
served copy with `scripts/molstar_figure.js` appended, and that script hands back Mol*'s
own offscreen render together with the image coordinates of the residues the figure
labels. Composing the labels over that image is `scripts/labeled_structure.py`; the two
are separate rules because the render takes minutes and the label placement is what gets
adjusted.

Nothing here re-creates what `prot-struct-viz` did. The page's coloring, camera, and
geometry are used as generated; only the surface type is changed, and only by updating
the representation nodes the page already built.

The `figure` wildcard says which page; the differences between them are the view to
render and which chains the sites to label are recorded on.
"""

import base64
import csv
import functools
import http.server
import json
import pathlib
import shutil
import subprocess
import sys
import threading

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# Per figure: the view it shows, and the chains whose coloring-CSV rows get an anchor.
# The page opens on the first view in its spec, and the injected script checks that this
# is the one rather than switching to it.
#
# The chains differ because the two entries reach a trimer differently. 8FAW deposits one
# protomer and Mol* generates the other two, so the H3 sites are recorded on chain A
# alone and its other rows are the bound LSTc analogue, which the figure draws but does
# not label. 9GSP deposits all three protomers, so the H1 sites are recorded once per
# protomer and all three chains are named. Either way a site ends up with three copies to
# choose between, and the chains named for one figure hold equivalent copies of the same
# sites -- which is checked below.
FIGURE_VIEWS = {
    "H3_subclade_K_mutations_structure": (
        "subclade-k-mutations-also-indicating-sites-222-and-223",
        ["A"],
    ),
    "H1_D_3_1_1_mutations_structure": (
        "d-3-1-1-mutations-also-indicating-sites-155-and-157",
        ["A", "B", "C"],
    ),
}

# The surface the figure draws in place of the one the interactive page draws: the
# solvent-excluded molecular surface rather than the smoother Gaussian one.
SURFACE_FROM = "gaussian-surface"
SURFACE_TO = "molecular-surface"

# Live canvas in CSS pixels, which fixes the figure's aspect ratio, and the whole-number
# factor the export is scaled up by. The export has to be an exact multiple of the live
# drawing buffer: that is what makes a projected anchor land on the exported image under
# a pure scale. This works out at roughly 2000 x 3000, which is more than 300 dpi at any
# size a page will print this figure.
CANVAS_CSS = (1000, 1500)
EXPORT_SCALE = 2

# White, to match the other figures.
TRANSPARENT_BACKGROUND = False

# Chromium renders WebGL in software here, which it will only do when told that a major
# performance caveat is acceptable; without the third flag Mol* reports no WebGL at all.
# The debugging port is what keeps the browser alive: headless Chromium given only a URL
# exits before the page's scripts run.
CHROMIUM = "chromium-browser"
CHROMIUM_FLAGS = [
    "--headless",
    "--no-sandbox",
    "--disable-dev-shm-usage",
    "--use-gl=angle",
    "--use-angle=swiftshader",
    "--override-use-software-gl-for-tests",
    "--remote-debugging-port=0",
]

# How long to wait for the browser, and how patiently the injected script waits on Mol*.
RENDER_TIMEOUT_SECONDS = 1800
POLL_MS = 250
LOAD_POLLS = 480
SETTLE_POLLS = 8


def labeled_sites(coloring_csv, chains):
    """The `[{chain, residue}]` of the coloring CSV's rows on `chains`.

    Rows for the same residue on different chains must agree on what the site is called
    and how it is colored, since they are meant to be copies of one another and the
    figure will label whichever copy faces the camera.
    """
    with open(coloring_csv, newline="") as f:
        rows = list(csv.DictReader(f))
    sites = [row for row in rows if row["chain"] in chains]
    if not sites:
        raise ValueError(f"{coloring_csv}: no rows on chain {chains}")

    by_residue = {}
    for row in sites:
        shown = (row["label"], row["color"])
        if by_residue.setdefault(row["residue"], shown) != shown:
            raise ValueError(
                f"{coloring_csv}: residue {row['residue']} is drawn as "
                f"{by_residue[row['residue']]} on one of {chains} and as {shown} on "
                f"another, so its chains are not copies of one another"
            )
    print(
        f"{coloring_csv}: {len(sites)} of {len(rows)} rows are on chain(s) "
        f"{', '.join(chains)}, covering {len(by_residue)} sites"
    )
    return [{"chain": row["chain"], "residue": row["residue"]} for row in sites]


def served_page(page_html, injected_js, view_slug, sites, directory):
    """Write `page_html` into `directory` with the figure script appended."""
    html = pathlib.Path(page_html).read_text()
    settings = {
        "view_slug": view_slug,
        "sites": sites,
        "surface_from": SURFACE_FROM,
        "surface_to": SURFACE_TO,
        "canvas_css": list(CANVAS_CSS),
        "export_scale": EXPORT_SCALE,
        "transparent": TRANSPARENT_BACKGROUND,
        "poll_ms": POLL_MS,
        "load_polls": LOAD_POLLS,
        "settle_polls": SETTLE_POLLS,
    }
    closing = "</body>"
    if html.count(closing) != 1:
        raise ValueError(f"{page_html}: found {html.count(closing)} {closing} tags")
    script = (
        f"<script>var PSV_FIGURE = {json.dumps(settings)};</script>\n"
        f"<script>\n{pathlib.Path(injected_js).read_text()}</script>\n"
    )
    page = directory / "page.html"
    page.write_text(html.replace(closing, script + closing))
    return page


def run_browser(directory, page_name):
    """Serve `directory` over loopback, run the page in Chromium, return what it posts.

    The page posts its result back to the server it was served from, so this returns
    once the figure is in hand -- no polling the browser and no guessing at how long
    Mol* needs for a surface this size.
    """
    if shutil.which(CHROMIUM) is None:
        raise RuntimeError(
            f"{CHROMIUM} is not on PATH; this figure renders a Mol* page in a browser"
        )
    posted = {}
    finished = threading.Event()

    class Handler(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass  # the browser's request log says nothing this rule's log needs

        def do_POST(self):
            body = self.rfile.read(int(self.headers["Content-Length"]))
            posted[self.path] = json.loads(body)
            self.send_response(204)
            self.end_headers()
            finished.set()

    handler = functools.partial(Handler, directory=str(directory))
    with http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler) as server:
        port = server.server_address[1]
        threading.Thread(target=server.serve_forever, daemon=True).start()
        command = [
            CHROMIUM,
            *CHROMIUM_FLAGS,
            f"--user-data-dir={directory / 'chromium'}",
            f"--window-size={CANVAS_CSS[0] + 400},{CANVAS_CSS[1] + 400}",
            f"http://127.0.0.1:{port}/{page_name}",
        ]
        print(f"running {' '.join(command)}")
        browser = subprocess.Popen(
            command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        try:
            if not finished.wait(timeout=RENDER_TIMEOUT_SECONDS):
                raise RuntimeError(
                    f"the page posted nothing back within {RENDER_TIMEOUT_SECONDS} s"
                )
        finally:
            browser.kill()
            browser.wait()

    if "/error" in posted:
        for line in posted["/error"]["progress"]:
            print(f"  {line}")
        raise RuntimeError(f"the page reported: {posted['/error']['error']}")
    if "/figure" not in posted:
        raise RuntimeError(f"the page posted to {sorted(posted)}, expected /figure")
    return posted["/figure"]


def projected_copies(anchors, sites):
    """Every copy of every site, where it projects to and how far away it is.

    Both specs show a trimer, so each site comes back three times over -- as symmetry
    copies of one chain, or as one row per deposited protomer. All of them are reported:
    which to label is a question of which one the render actually shows in color, and
    that is answered against the image by `scripts/labeled_structure.py`.
    """
    by_residue = {}
    for anchor in anchors:
        by_residue.setdefault(anchor["residue"], []).append(anchor)
    missing = sorted({site["residue"] for site in sites} - set(by_residue))
    if missing:
        raise ValueError(f"the page draws no residue {missing}")

    copies = []
    for residue, found in sorted(by_residue.items(), key=lambda kv: kv[0]):
        found = sorted(found, key=lambda copy: copy["distance"])
        where = ", ".join(
            f"{copy['chain']}/{copy['unit']} at ({copy['x']:.5g}, {copy['y']:.5g}), "
            f"{copy['distance']:.5g} from the camera"
            for copy in found
        )
        print(f"  {residue}: {len(found)} copies -- {where}")
        copies += [
            {
                "chain": copy["chain"],
                "residue": residue,
                "unit": copy["unit"],
                "x": round(copy["x"]),
                "y": round(copy["y"]),
                "camera_distance": f"{copy['distance']:.5g}",
            }
            for copy in found
        ]
    return copies


view_slug, chains = FIGURE_VIEWS[snakemake.wildcards.figure]
sites = labeled_sites(snakemake.input.coloring_csv, chains)
directory = (
    pathlib.Path(snakemake.output.structure_png).parent
    / f"_molstar_render_{snakemake.wildcards.figure}"
)
directory.mkdir(parents=True, exist_ok=True)
page = served_page(
    snakemake.input.page_html,
    snakemake.input.injected_js,
    view_slug,
    sites,
    directory,
)
figure = run_browser(directory, page.name)
for line in figure["progress"]:
    print(line)

anchors = projected_copies(figure["anchors"], sites)
prefix = "data:image/png;base64,"
if not figure["png"].startswith(prefix):
    raise ValueError(f"the page posted back {figure['png'][:40]!r}, expected a PNG")
png = base64.b64decode(figure["png"][len(prefix) :])
pathlib.Path(snakemake.output.structure_png).write_bytes(png)

with open(snakemake.output.anchors_csv, "w", newline="") as f:
    writer = csv.DictWriter(
        f, ["chain", "residue", "unit", "x", "y", "camera_distance"]
    )
    writer.writeheader()
    writer.writerows(anchors)

print(
    f"wrote {snakemake.output.structure_png} "
    f"({figure['image']['width']} x {figure['image']['height']}, {len(png) / 1e6:.3g} MB) "
    f"and {len(anchors)} projected copies to {snakemake.output.anchors_csv}"
)
