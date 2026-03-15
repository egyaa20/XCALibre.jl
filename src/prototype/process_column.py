#!/usr/bin/env python3
"""
plot_hydrostatic_column.py
==========================
Produces 3 vertical profile plots (at x=0.5, t=0.1s):
    1. p_rgh  vs  y
    2. p      vs  y
    3. |U|    vs  y

OpenFOAM: single combined file
    postProcessing_column/sample/0.1/verticalLine_p_p_rgh_U.xy
    columns:  y   p   p_rgh   Ux   Uy   Uz

XCALibre: time_1000.vtk (sampled at x=0.5)
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

try:
    import pyvista as pv
except ImportError:
    sys.exit("ERROR: pyvista not installed.  Run:  pip install pyvista")

# =============================================================================
# CONFIG
# =============================================================================
OF_SAMPLE_DIR = os.path.join("postProcessing_column", "sample", "0.1")
OF_SET_NAME   = "verticalLine"

VTK_FILE      = "time_1000.vtk"

XC_X        = 0.5
XC_Y_START  = 0.0
XC_Y_END    = 1.0
XC_N_POINTS = 200

OUT_DIR = "graphs"

C_OF = "#2196F3"
C_XC = "#FF6D00"
XC_HALO = [pe.Stroke(linewidth=7, foreground="white"), pe.Normal()]

# =============================================================================
# HELPERS
# =============================================================================
def setup_matplotlib():
    plt.rcParams.update({
        "figure.figsize": (7, 8), "font.size": 13,
        "axes.grid": True, "grid.alpha": 0.30,
        "legend.frameon": True, "legend.framealpha": 0.92,
        "font.family": "sans-serif",
    })

def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path}")

# =============================================================================
# OPENFOAM — combined file reader
# columns: y  p  p_rgh  Ux  Uy  Uz
# =============================================================================
def read_of_combined():
    path = os.path.join(OF_SAMPLE_DIR, f"{OF_SET_NAME}_p_p_rgh_U.xy")
    if not os.path.isfile(path):
        print(f"  WARNING: not found – {path}")
        if os.path.isdir(OF_SAMPLE_DIR):
            print(f"  Files present: {os.listdir(OF_SAMPLE_DIR)}")
        else:
            print(f"  Directory missing: {OF_SAMPLE_DIR}")
        return None
    data = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    data.append([float(v) for v in parts[:6]])
                except ValueError:
                    continue
    if not data:
        print(f"  WARNING: file empty – {path}")
        return None
    arr = np.array(data, dtype=float)
    print(f"  Loaded {path}  ({len(arr)} rows)")
    return {
        "y":     arr[:, 0],
        "p":     arr[:, 1],
        "p_rgh": arr[:, 2],
        "Umag":  np.sqrt(arr[:, 3]**2 + arr[:, 4]**2 + arr[:, 5]**2),
    }

# =============================================================================
# XCALIBRE — VTK reader
# =============================================================================
def read_mesh_any(filepath):
    m = pv.read(filepath)
    if isinstance(m, pv.MultiBlock):
        blocks = [b for b in m if b is not None]
        if not blocks:
            raise RuntimeError(f"{filepath}: empty MultiBlock")
        try:
            return pv.merge(blocks)
        except Exception:
            merged = blocks[0]
            for b in blocks[1:]:
                try:
                    merged = merged.merge(b)
                except Exception:
                    pass
            return merged
    return m

def ensure_point_field(mesh, field):
    if field in mesh.point_data:
        return mesh
    if field in mesh.cell_data:
        return mesh.cell_data_to_point_data()
    raise KeyError(f"'{field}' not found. Available: "
                   f"{list(mesh.point_data.keys()) + list(mesh.cell_data.keys())}")

def sample_xc_line(mesh, field):
    mesh = ensure_point_field(mesh, field)
    zmin, zmax = mesh.bounds[4], mesh.bounds[5]
    z_mid = 0.5 * (zmin + zmax)
    sampled = mesh.sample_over_line(
        (XC_X, XC_Y_START, z_mid),
        (XC_X, XC_Y_END,   z_mid),
        resolution=XC_N_POINTS - 1,
    )
    sampled = ensure_point_field(sampled, field)
    y    = np.asarray(sampled.points)[:, 1]
    vals = np.asarray(sampled.point_data[field], dtype=float)
    return y, vals

# =============================================================================
# MAIN
# =============================================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_matplotlib()

    # ── OpenFOAM ─────────────────────────────────────────────────────────── #
    print("\n── OpenFOAM ──────────────────────────────────")
    of = read_of_combined()
    of_y = of_p_rgh = of_p = of_Umag = None
    if of is not None:
        of_y, of_p_rgh, of_p, of_Umag = of["y"], of["p_rgh"], of["p"], of["Umag"]

    # ── XCALibre ─────────────────────────────────────────────────────────── #
    print("\n── XCALibre ──────────────────────────────────")
    xc_y_prgh = xc_prgh = xc_y_p = xc_p = xc_y_U = xc_Umag = None

    if not os.path.isfile(VTK_FILE):
        print(f"  WARNING: not found – {VTK_FILE}")
    else:
        mesh = read_mesh_any(VTK_FILE)
        print(f"  Loaded {VTK_FILE}")
        print(f"  Fields: {list(mesh.point_data.keys()) + list(mesh.cell_data.keys())}")

        for field, yvar, vvar in [("p_rgh", "xc_y_prgh", "xc_prgh"),
                                   ("p",     "xc_y_p",    "xc_p")]:
            try:
                y, v = sample_xc_line(mesh, field)
                locals()[yvar] = y
                locals()[vvar] = v.reshape(-1)
                print(f"  {field} sampled ({len(y)} pts)")
            except KeyError as e:
                print(f"  {e}")

        # re-sample individually to keep variables in scope
        try:
            xc_y_prgh, v = sample_xc_line(mesh, "p_rgh")
            xc_prgh = v.reshape(-1)
        except KeyError:
            pass
        try:
            xc_y_p, v = sample_xc_line(mesh, "p")
            xc_p = v.reshape(-1)
        except KeyError:
            pass
        try:
            xc_y_U, v = sample_xc_line(mesh, "U")
            xc_Umag = np.sqrt((v**2).sum(axis=1)) if v.ndim == 2 else np.abs(v)
            print(f"  U sampled ({len(xc_y_U)} pts)")
        except KeyError as e:
            print(f"  {e}")

    # ── Plot 1: p_rgh ────────────────────────────────────────────────────── #
    print("\n── Plotting ──────────────────────────────────")
    fig, ax = plt.subplots()
    ax.set_title("Dynamic Pressure  p_rgh  (t = 0.1 s)", fontsize=14)
    ax.set_xlabel("p_rgh  [Pa]", fontsize=12)
    ax.set_ylabel("y  [m]",      fontsize=12)
    if of_y is not None:
        ax.plot(of_p_rgh, of_y, color=C_OF, lw=3.5, ls="-", alpha=1.0, label="OpenFOAM", zorder=5)
    if xc_y_prgh is not None:
        ln, = ax.plot(xc_prgh, xc_y_prgh, color=C_XC, lw=5.0, ls="--", label="XCALibre", zorder=3)
        ln.set_path_effects(XC_HALO)
    ax.legend(fontsize=11)
    savefig(os.path.join(OUT_DIR, "profile_p_rgh.png"))

    # ── Plot 2: p ────────────────────────────────────────────────────────── #
    fig, ax = plt.subplots()
    ax.set_title("Total Pressure  p  (t = 0.1 s)", fontsize=14)
    ax.set_xlabel("p  [Pa]",  fontsize=12)
    ax.set_ylabel("y  [m]",   fontsize=12)
    if of_y is not None:
        ax.plot(of_p, of_y, color=C_OF, lw=3.5, ls="-", alpha=1.0, label="OpenFOAM", zorder=5)
    if xc_y_p is not None:
        ln, = ax.plot(xc_p, xc_y_p, color=C_XC, lw=5.0, ls="--", label="XCALibre", zorder=3)
        ln.set_path_effects(XC_HALO)
    ax.legend(fontsize=11)
    savefig(os.path.join(OUT_DIR, "profile_p.png"))

    # ── Plot 3: |U| ──────────────────────────────────────────────────────── #
    fig, ax = plt.subplots()
    ax.set_title("Velocity Magnitude  |U|  (t = 0.1 s)", fontsize=14)
    ax.set_xlabel("|U|  [m/s]", fontsize=12)
    ax.set_ylabel("y  [m]",     fontsize=12)
    if of_y is not None:
        ax.plot(of_Umag, of_y, color=C_OF, lw=3.5, ls="-", alpha=1.0, label="OpenFOAM", zorder=5)
    if xc_y_U is not None:
        ln, = ax.plot(xc_Umag, xc_y_U, color=C_XC, lw=5.0, ls="--", label="XCALibre", zorder=3)
        ln.set_path_effects(XC_HALO)
    ax.legend(fontsize=11)
    savefig(os.path.join(OUT_DIR, "profile_U_mag.png"))

    print(f"\n✅  Done. Plots saved to '{OUT_DIR}/'")


if __name__ == "__main__":
    main()