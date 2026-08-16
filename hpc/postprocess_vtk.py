#!/usr/bin/env python3
"""Post-process a single legacy .vtk file written by XCALibre.jl.

Usage: postprocess_vtk.py <input.vtk> <outdir> [--field p] [--decimate 0.75]

Writes into <outdir>:
  <field>_xy.png      static view of the field at the X-Y plane
  <field>_scene.html  self-contained interactive scene (rotate/zoom/pan)

Runs off-screen (no display needed) on a compute node. The PNG is the
proven path and is written first; the HTML export is best-effort, so a
missing trame backend degrades to "no interactive scene" rather than
losing the static image too.

Decimation notes: .vtk fields from XCALibre are CELL data on an
unstructured volume grid, which vtkQuadricDecimation cannot consume
directly. The surface has to be extracted and triangulated first, and
cell data converted to point data BEFORE decimation so the field
interpolates across merged triangles instead of being dropped.
"""
import argparse
import sys

import pyvista as pv

pv.OFF_SCREEN = True


def render_png(mesh, field: str, png_path: str, flat: bool) -> None:
    """Static view at the X-Y plane; slices mid-Z if the mesh has depth."""
    plotter = pv.Plotter(off_screen=True, window_size=(1000, 800))
    if flat:
        plotter.add_mesh(mesh, scalars=field, cmap="viridis", show_edges=False)
    else:
        center = mesh.center
        plane = mesh.slice(normal="z", origin=center)
        plotter.add_mesh(plane, scalars=field, cmap="viridis", show_edges=False)
    plotter.view_xy()
    plotter.add_scalar_bar(title=field)
    plotter.screenshot(png_path)
    plotter.close()


# Below this many triangles, decimation is pure quality loss: the ~1 MB
# embedded VTK.js viewer dwarfs the geometry's contribution to file size,
# so there's nothing to gain, and coarsening a mesh this small visibly
# blocks up any sharp gradient (e.g. behind a BFS step).
MIN_CELLS_FOR_DECIMATION = 5000


def prepare_surface(mesh, field: str, reduction: float):
    """Cell-data volume grid -> decimated triangulated surface carrying `field`.

    Falls back to the undecimated surface if decimation drops all cells,
    which vtkQuadricDecimation can do on degenerate/flat geometry.
    """
    surf = (mesh.cell_data_to_point_data()
                .extract_surface(algorithm="dataset_surface")
                .triangulate())

    if reduction <= 0.0:
        return surf
    if surf.n_cells < MIN_CELLS_FOR_DECIMATION:
        print(f"skipping decimation: surface has {surf.n_cells} triangles "
              f"(< {MIN_CELLS_FOR_DECIMATION}), not worth the quality loss")
        return surf
    try:
        decimated = surf.decimate(reduction)
    except Exception as exc:  # noqa: BLE001 - keep the surface, report why
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


def export_html(mesh, field: str, reduction: float, html_path: str, flat: bool) -> None:
    surf = prepare_surface(mesh, field, reduction)

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(surf, scalars=field, cmap="viridis", show_edges=False)
    plotter.add_scalar_bar(title=field)
    # Matches render_png: an isometric camera on an essentially flat mesh
    # just tilts a thin plane into a confusing sliver. Genuinely 3D fields
    # still get isometric so their shape is visible at all.
    plotter.view_xy() if flat else plotter.view_isometric()
    plotter.export_html(html_path)
    plotter.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vtk_path")
    ap.add_argument("outdir")
    ap.add_argument("--field", default="p")
    ap.add_argument("--decimate", type=float, default=0.75,
                    help="fraction of triangles to remove before HTML export")
    args = ap.parse_args()

    mesh = pv.read(args.vtk_path)
    field = args.field

    if field not in mesh.cell_data and field not in mesh.point_data:
        raise SystemExit(
            f"no '{field}' field in {args.vtk_path}; "
            f"cell_data={list(mesh.cell_data.keys())} "
            f"point_data={list(mesh.point_data.keys())}")

    zmin, zmax = mesh.bounds[4], mesh.bounds[5]
    flat = (zmax - zmin) < 1e-9

    png_path = f"{args.outdir}/{field}_xy.png"
    render_png(mesh, field, png_path, flat)
    print(f"wrote {png_path}")

    # Best-effort: never let a missing trame backend cost us the PNG.
    html_path = f"{args.outdir}/{field}_scene.html"
    try:
        export_html(mesh, field, args.decimate, html_path, flat)
        print(f"wrote {html_path}")
    except Exception as exc:  # noqa: BLE001
        print(f"WARN interactive HTML export failed: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
