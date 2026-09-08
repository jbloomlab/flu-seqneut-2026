/* Drive the Mol* page a `prot-struct-viz` spec produced, and hand back a figure.
 *
 * Appended to a served copy of the page by `scripts/molstar_render.py`, which defines
 * `PSV_FIGURE` ahead of it; every choice this makes comes from there, so the mechanics
 * here are the same for any such page. Two things go back: the high-resolution image
 * Mol* renders offscreen, and where on that image each requested residue sits, so that
 * labels can be drawn over it as vector text rather than buried in the surface.
 *
 * Results are POSTed back rather than written anywhere, because the page is served over
 * loopback and the caller is waiting on the request; nothing here has to guess how long
 * Mol* needs.
 */
(async function () {
    var post = function (path, body) {
        return fetch(path, {method: "POST", body: body});
    };
    var progress = [];
    var note = function (message) {
        progress.push(message);
    };
    var sleep = function (ms) {
        return new Promise(function (resolve) { setTimeout(resolve, ms); });
    };

    try {
        // The page clears #status when its own load chain resolves, and puts the failure
        // there if it does not.
        var status = document.getElementById("status");
        for (var i = 0; i < PSV_FIGURE.load_polls && status.textContent !== ""; i++) {
            await sleep(PSV_FIGURE.poll_ms);
        }
        if (status.textContent !== "") {
            throw new Error("page did not load: " + status.textContent);
        }
        var plugin = window.viewer.plugin;

        // The page opens on the first view in the spec. Check that is the view the figure
        // is about rather than selecting it, so a reordered spec fails here instead of
        // quietly changing which view the paper shows.
        var select = document.getElementById("view-select");
        var shown = select ? select.value : PSV_FIGURE.view_slug;
        if (shown !== PSV_FIGURE.view_slug) {
            throw new Error("page opens on view '" + shown + "', not '" + PSV_FIGURE.view_slug + "'");
        }

        // Update the surface representations in place. Adding one instead would arrive in
        // Mol*'s own element coloring, because MVS hangs the annotation colors off each
        // representation node it creates; updating keeps that node and its color theme.
        var build = plugin.build();
        var swapped = [];
        plugin.state.data.cells.forEach(function (cell, ref) {
            var values = cell.params && cell.params.values;
            if (values && values.type && values.type.name === PSV_FIGURE.surface_from) {
                swapped.push(String(ref));
                build.to(ref).update(function (old) {
                    old.type = {name: PSV_FIGURE.surface_to, params: {}};
                });
            }
        });
        if (swapped.length === 0) {
            throw new Error("page draws no '" + PSV_FIGURE.surface_from + "' representation");
        }
        await build.commit();
        plugin.state.data.cells.forEach(function (cell) {
            var values = cell.params && cell.params.values;
            if (values && values.type && values.type.name === PSV_FIGURE.surface_from) {
                throw new Error("a '" + PSV_FIGURE.surface_from + "' representation survived the swap");
            }
        });
        note("swapped " + swapped.length + " " + PSV_FIGURE.surface_from
             + " representations for " + PSV_FIGURE.surface_to);

        // Give the live canvas the figure's aspect ratio. The export is then an exact
        // multiple of it, which is what lets a projected anchor be scaled onto the
        // exported image instead of re-derived against a different camera.
        var box = document.getElementById("viewer");
        box.style.width = PSV_FIGURE.canvas_css[0] + "px";
        box.style.height = PSV_FIGURE.canvas_css[1] + "px";
        plugin.canvas3d.handleResize();
        var gl = plugin.canvas3d.webgl.gl;
        var size = null;
        var settled = 0;
        for (var j = 0; j < PSV_FIGURE.load_polls && settled < PSV_FIGURE.settle_polls; j++) {
            await sleep(PSV_FIGURE.poll_ms);
            var now = [gl.drawingBufferWidth, gl.drawingBufferHeight];
            settled = (size && size[0] === now[0] && size[1] === now[1]) ? settled + 1 : 0;
            size = now;
        }
        if (settled < PSV_FIGURE.settle_polls) {
            throw new Error("canvas size never settled, last " + size);
        }
        var scale = PSV_FIGURE.export_scale;
        var image = [size[0] * scale, size[1] * scale];
        note("canvas " + size[0] + "x" + size[1] + ", exporting " + image[0] + "x" + image[1]);

        // Where each requested residue is on the image. Coordinates come from the loaded
        // assembly through Mol*'s own accessors -- `unit.conformation.position` applies
        // the assembly operator -- so each symmetry copy of a residue answers separately
        // and the caller decides which copy to label.
        //
        // Each copy is placed at its atom nearest the camera rather than at the centroid
        // of its atoms. A residue reaching into the protein has a centroid under the
        // surface, which projects to a point beside the colored patch instead of on it;
        // the nearest atom is the part of the residue the surface actually shows.
        var lib = molstar.lib.structure;
        var props = lib.StructureProperties;
        var structure = null;
        plugin.state.data.cells.forEach(function (cell) {
            if (!structure && cell.transform.transformer.id === "ms-plugin.structure-from-model"
                && cell.obj) {
                structure = cell.obj.data;
            }
        });
        if (!structure) {
            throw new Error("page holds no structure node");
        }
        var wanted = new Set(PSV_FIGURE.sites.map(function (site) {
            return site.chain + "/" + site.residue;
        }));
        var found = new Set();
        var camera = plugin.canvas3d.camera;
        var eye = camera.position;
        var location = lib.StructureElement.Location.create(structure);
        var nearest = new Map();
        structure.units.forEach(function (unit) {
            location.unit = unit;
            var elements = unit.elements;
            for (var k = 0; k < elements.length; k++) {
                location.element = elements[k];
                var residue = String(props.residue.auth_seq_id(location))
                    + (props.residue.pdbx_PDB_ins_code(location) || "");
                var chain = props.chain.auth_asym_id(location);
                var site = chain + "/" + residue;
                if (!wanted.has(site)) continue;
                found.add(site);
                var point = [0, 0, 0];
                unit.conformation.position(elements[k], point);
                var dx = point[0] - eye[0];
                var dy = point[1] - eye[1];
                var dz = point[2] - eye[2];
                var distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                var key = site + "@" + unit.id;
                var best = nearest.get(key);
                if (!best) {
                    best = {chain: chain, residue: residue, unit: unit.id, atoms: 0,
                            point: point, distance: distance};
                    nearest.set(key, best);
                } else if (distance < best.distance) {
                    best.point = point;
                    best.distance = distance;
                }
                best.atoms += 1;
            }
        });
        var unseen = [];
        wanted.forEach(function (site) {
            if (!found.has(site)) unseen.push(site);
        });
        if (unseen.length) {
            throw new Error("the page draws no residue " + unseen.join(", "));
        }

        var anchors = [];
        nearest.forEach(function (best) {
            var projected = [0, 0, 0, 0];
            camera.project(projected, best.point);
            anchors.push({
                chain: best.chain,
                residue: best.residue,
                unit: best.unit,
                atoms: best.atoms,
                // `camera.project` puts the origin bottom left, as WebGL does; an image
                // has it top left.
                x: projected[0] * scale,
                y: image[1] - projected[1] * scale,
                distance: best.distance,
            });
        });

        // Mol*'s exporter re-frames the image to fit its content by default, which would
        // throw away the camera the spec's `orientation:` block authored and move every
        // anchor. Off means the export is the live camera at a larger size.
        var shot = plugin.helpers.viewportScreenshot;
        shot.behaviors.cropParams.next({auto: false, relativePadding: 0});
        shot.resetCrop();
        shot.behaviors.values.next(Object.assign({}, shot.values, {
            transparent: PSV_FIGURE.transparent,
            resolution: {name: "custom", params: {width: image[0], height: image[1]}},
        }));
        var uri = await shot.getImageDataUri();

        await post("/figure", JSON.stringify({
            png: uri,
            image: {width: image[0], height: image[1]},
            anchors: anchors,
            progress: progress,
        }));
    } catch (err) {
        await post("/error", JSON.stringify({
            error: String((err && err.stack) || err),
            progress: progress,
        }));
    }
})();
