#!/usr/bin/env python3
"""Render the pressure ("p") cell-data field at the X-Y plane from a
single legacy .vtk file written by XCALibre.jl, saving a PNG.

Usage: postprocess_vtk.py <input.vtk> <output.png>

Runs off-screen (no display needed) on a compute node. If the mesh is
already flat in Z (typical for XCALibre's 2D cases), it's viewed directly
from +Z; otherwise it's sliced at the mid-Z plane first.
"""
import sys

import pyvista as pv

pv.OFF_SCREEN = True


def main(vtk_path: str, png_path: str) -> None:
    mesh = pv.read(vtk_path)

    if "p" not in mesh.cell_data:
        raise SystemExit(f"no 'p' field in {vtk_path}; available: {list(mesh.cell_data.keys())}")

    zmin, zmax = mesh.bounds[4], mesh.bounds[5]
    flat = (zmax - zmin) < 1e-9

    plotter = pv.Plotter(off_screen=True, window_size=(1000, 800))
    if flat:
        plotter.add_mesh(mesh, scalars="p", cmap="viridis", show_edges=False)
        plotter.view_xy()
    else:
        mid_z = (zmin + zmax) / 2.0
        center = mesh.center
        plane = mesh.slice(normal="z", origin=(center[0], center[1], mid_z))
        plotter.add_mesh(plane, scalars="p", cmap="viridis", show_edges=False)
        plotter.view_xy()

    plotter.add_scalar_bar(title="p")
    plotter.screenshot(png_path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} <input.vtk> <output.png>")
    main(sys.argv[1], sys.argv[2])
