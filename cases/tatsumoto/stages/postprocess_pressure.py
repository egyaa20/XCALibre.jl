#!/usr/bin/env python3
"""Render the pressure field from a tatsumoto heated-run VTU snapshot.

Writes two things into <outdir>:
  <field>_xsection.png  cross-section perpendicular to flow (Z-normal, at
                         the domain's axial midpoint) -- the O-grid
                         quarter-circle, showing radial/BL structure
  <field>_scene.html     self-contained interactive scene of the full
                         (decimated) outer surface (rotate/zoom/pan)

Two whole-pipe PNG views were tried and abandoned first (a longitudinal
slice, then an isometric surface view): both collapse to an unreadable
sliver, because the pipe genuinely is that shape (300 mm long by 2.7 mm
radius, over 100:1) -- no camera angle fixes an aspect ratio that
extreme. A cross-section sidesteps it entirely (both in-plane dimensions
are ~R) and is the more physically informative view besides. The HTML
covers whole-pipe inspection instead, where rotate/zoom actually works.

Field defaults to "p" (full reconstructed pressure) rather than "p_rgh"
(the solver's internal reduced/hydrodynamic pressure) -- both exist in the
VTU, "p" is the physically meaningful one for inspection.

Same PNG-first-then-HTML-best-effort pattern as hpc/postprocess_vtk.py: a
missing trame backend costs the interactive scene, never the static image.
Decimation/resampling logic is identical to that script (see it for the
"why" on the vtkQuadricDecimation point-data-drop fix) -- kept as a
separate copy here rather than a shared import because the two scripts
assume different geometries (flat 2D vs this 3D quarter-pipe) and have no
other code in common.

Usage: postprocess_pressure.py <input.vtu> <outdir> [--field p] [--decimate 0.75]
"""
import argparse
import sys

import pyvista as pv

pv.OFF_SCREEN = True

MIN_CELLS_FOR_DECIMATION = 5000


def render_png(mesh, field: str, png_path: str) -> None:
    """Cross-section perpendicular to flow (Z-normal), at the domain's
    axial midpoint.

    Two whole-pipe views were tried and abandoned first: a longitudinal
    slice and an isometric surface view both collapse to an unreadable
    sliver, because the pipe genuinely IS that shape (300 mm long by
    2.7 mm radius, over 100:1) -- no camera angle fixes an aspect ratio
    that extreme. A cross-section has none of that problem (its two
    in-plane dimensions are both ~R), and it's the physically informative
    view anyway: this is exactly the O-grid quarter-circle the mesh
    generator builds, showing the radial/BL pressure distribution at one
    axial station. The interactive HTML covers the whole-pipe view instead
    (rotate/zoom lets you inspect any region there).
    """
    center = mesh.center
    z_min, z_max = mesh.bounds[4], mesh.bounds[5]
    # Not exactly at the midpoint: the axial extrusion is on 1 mm layer
    # boundaries (n_axial segments over the domain length), and the exact
    # midpoint lands ON one of those boundaries whenever the domain length
    # in mm is even -- cutting exactly along existing mesh geometry rather
    # than through cell interiors, which produced scattered degenerate
    # polygons (dark blotches, near-zero interpolated values) instead of a
    # clean cross-section. An off-integer fractional offset avoids it.
    z_cut = z_min + 0.5013 * (z_max - z_min)
    plane = mesh.slice(normal=(0, 0, 1), origin=(center[0], center[1], z_cut))
    if plane.n_points == 0:
        raise RuntimeError(f"Z-normal slice at z={z_cut} produced no points "
                            f"(domain z range: {z_min}..{z_max})")

    plotter = pv.Plotter(off_screen=True, window_size=(900, 900))
    plotter.add_mesh(plane, scalars=field, cmap="viridis", show_edges=False)
    plotter.view_xy()
    plotter.add_scalar_bar(title=field)
    plotter.screenshot(png_path)
    plotter.close()


def prepare_surface(mesh, field: str, reduction: float):
    """Cell-data volume grid -> decimated triangulated surface carrying `field`."""
    surf = (mesh.cell_data_to_point_data()
                .extract_surface(algorithm="dataset_surface")
                .triangulate())

    if reduction <= 0.0:
        return surf
    if surf.n_cells < MIN_CELLS_FOR_DECIMATION:
        print(f"skipping decimation: surface has {surf.n_cells} triangles "
              f"(< {MIN_CELLS_FOR_DECIMATION})")
        return surf
    try:
        decimated = surf.decimate(reduction)
    except Exception as exc:  # noqa: BLE001
        print(f"WARN decimation failed ({exc}), exporting full surface", file=sys.stderr)
        return surf
    if decimated.n_cells == 0:
        print("WARN decimation produced an empty mesh, exporting full surface", file=sys.stderr)
        return surf

    # vtkQuadricDecimation rebuilds the point set and drops point data with
    # it, so the field has to be interpolated back onto the reduced mesh.
    if field not in decimated.point_data:
        decimated = decimated.sample(surf)
        if field not in decimated.point_data:
            print(f"WARN could not carry '{field}' onto the decimated mesh, "
                  "exporting full surface", file=sys.stderr)
            return surf

    print(f"decimated {surf.n_cells} -> {decimated.n_cells} triangles "
          f"(reduction={reduction})")
    return decimated


def export_html(surf, field: str, html_path: str) -> None:
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(surf, scalars=field, cmap="viridis", show_edges=False)
    plotter.add_scalar_bar(title=field)
    plotter.view_isometric()
    plotter.export_html(html_path)
    plotter.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vtu_path")
    ap.add_argument("outdir")
    ap.add_argument("--field", default="p")
    ap.add_argument("--decimate", type=float, default=0.75)
    args = ap.parse_args()

    mesh = pv.read(args.vtu_path)
    field = args.field

    if field not in mesh.cell_data and field not in mesh.point_data:
        raise SystemExit(
            f"no '{field}' field in {args.vtu_path}; "
            f"cell_data={list(mesh.cell_data.keys())} "
            f"point_data={list(mesh.point_data.keys())}")

    png_path = f"{args.outdir}/{field}_xsection.png"
    render_png(mesh, field, png_path)
    print(f"wrote {png_path}")

    surf = prepare_surface(mesh, field, args.decimate)
    html_path = f"{args.outdir}/{field}_scene.html"
    try:
        export_html(surf, field, html_path)
        print(f"wrote {html_path}")
    except Exception as exc:  # noqa: BLE001
        print(f"WARN interactive HTML export failed: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
