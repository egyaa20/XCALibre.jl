#!/usr/bin/env python3
import os
import re
import math
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

# =========================
# CONFIG
# =========================
VTK_PREFIX = "time_"
VTK_EXT = ".vtk"

REF_DIR = "ref_data"
OF_POST_DIR = "postProcessing"
OUT_DIR = "graphs"

# Physics / scaling (matches your gnuplot)
rho = 1000.0
g = 9.81
H = 0.3
t_scale = math.sqrt(g / H)          # t_nd = t * sqrt(g/H)
p_scale = 1.0 / (rho * g * H)       # p_nd = p / (rho*g*H)

# Time mapping for XCALibre (VTK): time_14000 -> 1.4 s
FINAL_STEP = 14000
FINAL_TIME_S = 1.4
dt = FINAL_TIME_S / FINAL_STEP      # 1e-4 s

# Fields in VTK
FIELD_ALPHA = "alpha"
FIELD_P = "p"

# Sampling definitions (copied from your system/sampling)
LINES = {
    "L1": {"start": (0.3,   0.0, 0.0), "end": (0.3,   0.6, 0.0), "n": 100},
    "L2": {"start": (1.114, 0.0, 0.0), "end": (1.114, 0.6, 0.0), "n": 100},
    "L3": {"start": (1.362, 0.0, 0.0), "end": (1.362, 0.6, 0.0), "n": 100},
    "ToeLine": {"start": (0.0, 1e-4, 0.0), "end": (1.61, 1e-4, 0.0), "n": 200},
}

PROBES = {
    "right":  [(1.6099, 0.003, 0.0), (1.6099, 0.03, 0.0), (1.6099, 0.08, 0.0)],
    "left":   [(1e-4,   0.003, 0.0), (1e-4,   0.03, 0.0), (1e-4,   0.08, 0.0)],
    "bottom": [(0.7,    1e-4,  0.0), (1.0,    1e-4, 0.0), (1.3,    1e-4, 0.0)],
}

# Palette (close to your gnuplot)
c_ref1 = "#1f77b4"  # Blue
c_ref2 = "#ff7f0e"  # Orange
c_ref3 = "#d62728"  # Red

# For fig8/9 experiments (P1/P2/P3 colors)
c_ref_p1 = "#d62728"  # Red
c_ref_p2 = "#9467bd"  # Purple
c_ref_p3 = "#7f7f7f"  # Gray

# OpenFOAM & XCALibre base
c_of_line = "#444444"
c_xc_line = "black"

# Use consistent probe colours in multi-probe plots (match the experiment probe colours)
probe_cols_lr = [c_ref_p1, c_ref_p2, c_ref_p3]

# =========================
# HELPERS (generic)
# =========================
def get_step_number(filename: str) -> int:
    m = re.search(rf"{re.escape(VTK_PREFIX)}(\d+)", filename)
    return int(m.group(1)) if m else -1

def list_vtk_files() -> list[str]:
    files = []
    for f in os.listdir("."):
        if f.startswith(VTK_PREFIX) and f.endswith(VTK_EXT):
            step = get_step_number(f)
            if step >= 0:
                files.append(f)
    files.sort(key=get_step_number)
    return files

def ensure_point_field(mesh: pv.DataSet, field: str) -> pv.DataSet:
    if field in mesh.point_data:
        return mesh
    if field in mesh.cell_data:
        return mesh.cell_data_to_point_data()
    raise KeyError(f"Field '{field}' not found in point_data or cell_data.")

def read_mesh_any(filepath: str) -> pv.DataSet:
    m = pv.read(filepath)
    if isinstance(m, pv.MultiBlock):
        blocks = [b for b in m if b is not None]
        if not blocks:
            raise RuntimeError(f"{filepath}: MultiBlock has no blocks.")
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

def height_from_alpha_profile(y: np.ndarray, a: np.ndarray, threshold=0.5) -> float:
    if len(y) < 2:
        return float("nan")

    if np.nanmax(a) < threshold:
        return 0.0
    if np.nanmin(a) >= threshold:
        return float(np.nanmax(y))

    a_prev, y_prev = a[0], y[0]
    for i in range(1, len(y)):
        a_cur, y_cur = a[i], y[i]
        if (a_prev >= threshold) and (a_cur < threshold):
            denom = (a_cur - a_prev)
            if abs(denom) < 1e-14:
                return float(y_prev)
            r = (threshold - a_prev) / denom
            return float(y_prev + r * (y_cur - y_prev))
        a_prev, y_prev = a_cur, y_cur

    mask = a >= threshold
    return float(np.nanmax(y[mask])) if np.any(mask) else 0.0

def toe_from_alpha_profile(x: np.ndarray, a: np.ndarray, threshold=0.5) -> float:
    mask = a >= threshold
    return float(np.nanmax(x[mask])) if np.any(mask) else float("nan")

def load_ref_xy(path: str) -> tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(xs, dtype=float), np.array(ys, dtype=float)

def setup_matplotlib():
    plt.rcParams.update({
        "figure.figsize": (10, 7),
        "font.size": 14,
        "axes.grid": True,
        "grid.alpha": 0.30,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
    })

def savefig(path: str):
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()

# =========================
# PLOT STYLE HELPERS
# =========================
XC_HALO = [pe.Stroke(linewidth=6, foreground="white"), pe.Normal()]

def plot_experiment(ax, x, y, color, label):
    # thin faint line + open circle markers (points always readable)
    ax.plot(x, y, color=color, lw=1.2, alpha=0.35, zorder=3, label="_nolegend_")
    ax.plot(
        x, y, linestyle="none", marker="o", ms=5,
        mfc="white", mec=color, mew=1.4,
        color=color, label=label, zorder=5
    )

def plot_openfoam(ax, x, y, label="OpenFOAM"):
    ax.plot(x, y, color=c_of_line, lw=2.2, ls="--", alpha=0.85, label=label, zorder=2)

def plot_xcalibre(ax, x, y, label="XCALibre"):
    (ln,) = ax.plot(x, y, color=c_xc_line, lw=3.6, ls="-", label=label, zorder=4)
    ln.set_path_effects(XC_HALO)

def plot_openfoam_probe(ax, x, y, color):
    ax.plot(x, y, color=color, lw=2.0, ls="--", alpha=0.80, label="_nolegend_", zorder=2)

def plot_xcalibre_probe(ax, x, y, color):
    (ln,) = ax.plot(x, y, color=color, lw=3.4, ls="-", label="_nolegend_", zorder=4)
    ln.set_path_effects(XC_HALO)

def add_solver_probe_legends(ax, exp_loc="upper right", style_loc="lower right"):
    # Legend 1: experiments (already labelled)
    handles_exp, labels_exp = ax.get_legend_handles_labels()
    leg1 = ax.legend(handles_exp, labels_exp, loc=exp_loc, title="Experiments")
    ax.add_artist(leg1)

    # Legend 2: solver styles + probe colours
    handles_style = [
        Line2D([0], [0], color="k", lw=2.2, ls="--", label="OpenFOAM"),
        Line2D([0], [0], color="k", lw=3.4, ls="-", label="XCALibre"),
        Line2D([0], [0], color=probe_cols_lr[0], lw=3.0, label="P1"),
        Line2D([0], [0], color=probe_cols_lr[1], lw=3.0, label="P2"),
        Line2D([0], [0], color=probe_cols_lr[2], lw=3.0, label="P3"),
    ]
    ax.legend(handles=handles_style, loc=style_loc, ncol=2, title="Models / probes")

# =========================
# XCALibre (VTK) extraction
# =========================
def sample_line_vtk(mesh: pv.DataSet, start, end, n: int, field: str):
    sampled = mesh.sample_over_line(start, end, resolution=max(n - 1, 1))
    sampled = ensure_point_field(sampled, field)
    pts = np.asarray(sampled.points)
    vals = np.asarray(sampled.point_data[field]).astype(float)
    return pts, vals

def probe_points_vtk(mesh: pv.DataSet, points, field: str) -> np.ndarray:
    mesh = ensure_point_field(mesh, field)

    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    L = max(xmax - xmin, ymax - ymin, zmax - zmin, 1.0)
    eps = 1e-9 * L  # tiny nudge

    pts = np.asarray(points, dtype=float).reshape(-1, 3)

    if (zmax - zmin) < 1e-12:
        pts[:, 2] = zmin
    else:
        pts[:, 2] = 0.5 * (zmin + zmax)

    pts[:, 0] = np.clip(pts[:, 0], xmin + eps, xmax - eps)
    pts[:, 1] = np.clip(pts[:, 1], ymin + eps, ymax - eps)
    if (zmax - zmin) >= 1e-12:
        pts[:, 2] = np.clip(pts[:, 2], zmin + eps, zmax - eps)

    cloud = pv.PolyData(pts)
    probed = cloud.sample(mesh)

    probed = ensure_point_field(probed, field)
    vals = np.asarray(probed.point_data[field], dtype=float).reshape(-1)

    mask = probed.point_data.get("vtkValidPointMask", None)
    if mask is not None:
        mask = np.asarray(mask).astype(int).reshape(-1)
        bad = np.where(mask == 0)[0]
        if bad.size:
            pdat = np.asarray(mesh.point_data[field], dtype=float).reshape(-1)
            for i in bad:
                idx = mesh.find_closest_point(pts[i])
                vals[i] = pdat[idx]
            print(f"⚠️ Invalid probes: {bad.size}/{len(vals)} (used nearest-point fallback)")

    return vals

def extract_xcalibre_timeseries(vtk_files: list[str]) -> pd.DataFrame:
    records = []
    for f in vtk_files:
        step = get_step_number(f)
        t = step * dt

        mesh = read_mesh_any(f)
        mesh = ensure_point_field(mesh, FIELD_ALPHA)
        mesh = ensure_point_field(mesh, FIELD_P)

        h_vals = {}
        for name in ["L1", "L2", "L3"]:
            start, end, n = LINES[name]["start"], LINES[name]["end"], LINES[name]["n"]
            pts, a = sample_line_vtk(mesh, start, end, n, FIELD_ALPHA)
            y = pts[:, 1]
            h_vals[name] = height_from_alpha_profile(y, a, threshold=0.5)

        start, end, n = LINES["ToeLine"]["start"], LINES["ToeLine"]["end"], LINES["ToeLine"]["n"]
        pts, a_toe = sample_line_vtk(mesh, start, end, n, FIELD_ALPHA)
        x = pts[:, 0]
        toe = toe_from_alpha_profile(x, a_toe, threshold=0.5)

        p_right = probe_points_vtk(mesh, PROBES["right"], FIELD_P)
        p_left = probe_points_vtk(mesh, PROBES["left"], FIELD_P)
        p_bottom = probe_points_vtk(mesh, PROBES["bottom"], FIELD_P)

        rec = {
            "step": step,
            "t": t,
            "t_nd": t * t_scale,
            "h_L1": h_vals["L1"],
            "h_L2": h_vals["L2"],
            "h_L3": h_vals["L3"],
            "toe_x": toe,
            "p_right_1": p_right[0], "p_right_2": p_right[1], "p_right_3": p_right[2],
            "p_left_1":  p_left[0],  "p_left_2":  p_left[1],  "p_left_3":  p_left[2],
            "p_bottom_1": p_bottom[0], "p_bottom_2": p_bottom[1], "p_bottom_3": p_bottom[2],
        }
        records.append(rec)
        print(f"✅ XCALibre {f}  step={step:>6d}  t={t:.6f}s")

    return pd.DataFrame(records).sort_values("step").reset_index(drop=True)

# =========================
# OpenFOAM extraction
# =========================
def _is_number_dir(name: str) -> bool:
    try:
        float(name)
        return True
    except Exception:
        return False

def _read_two_column_numeric(path: str) -> tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    if not os.path.isfile(path):
        return np.array([]), np.array([])
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(xs, dtype=float), np.array(ys, dtype=float)

def extract_openfoam_heights_and_toe(of_dir: str):
    base = os.path.join(of_dir, "sampleSets")
    if not os.path.isdir(base):
        return None

    time_dirs = [d for d in os.listdir(base) if _is_number_dir(d)]
    time_dirs.sort(key=lambda s: float(s))

    t_list = []
    hL1_list, hL2_list, hL3_list = [], [], []
    toe_list = []

    for td in time_dirs:
        t = float(td)

        hs = []
        for name in ["L1", "L2", "L3"]:
            fpath = os.path.join(base, td, f"{name}_alpha.water.xy")
            y, a = _read_two_column_numeric(fpath)
            if y.size == 0:
                hs.append(np.nan)
                continue
            hs.append(height_from_alpha_profile(y, a, threshold=0.5))

        ftoe = os.path.join(base, td, "ToeLine_alpha.water.xy")
        x, a_toe = _read_two_column_numeric(ftoe)
        toe = toe_from_alpha_profile(x, a_toe, threshold=0.5) if x.size else np.nan

        t_list.append(t)
        hL1_list.append(hs[0])
        hL2_list.append(hs[1])
        hL3_list.append(hs[2])
        toe_list.append(toe)

    t_arr = np.array(t_list, dtype=float)
    return {
        "t": t_arr,
        "t_nd": t_arr * t_scale,
        "h_L1": np.array(hL1_list, dtype=float),
        "h_L2": np.array(hL2_list, dtype=float),
        "h_L3": np.array(hL3_list, dtype=float),
        "toe_x": np.array(toe_list, dtype=float),
    }

def extract_openfoam_probe_file(path: str):
    if not os.path.isfile(path):
        return None

    t_list = []
    p1_list, p2_list, p3_list = [], [], []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                t_list.append(float(parts[0]))
                p1_list.append(float(parts[1]))
                p2_list.append(float(parts[2]))
                p3_list.append(float(parts[3]))
            except ValueError:
                continue

    if not t_list:
        return None

    t = np.array(t_list, dtype=float)
    p1 = np.array(p1_list, dtype=float)
    p2 = np.array(p2_list, dtype=float)
    p3 = np.array(p3_list, dtype=float)

    return {
        "t": t,
        "t_nd": t * t_scale,
        "p1": p1, "p2": p2, "p3": p3,
        "p1_nd": p1 * p_scale, "p2_nd": p2 * p_scale, "p3_nd": p3 * p_scale,
    }

def extract_openfoam_all(of_dir: str):
    if not os.path.isdir(of_dir):
        return None

    heights_toe = extract_openfoam_heights_and_toe(of_dir)
    pr = extract_openfoam_probe_file(os.path.join(of_dir, "rightWallP", "0", "p"))
    pl = extract_openfoam_probe_file(os.path.join(of_dir, "leftWallP", "0", "p"))
    pb = extract_openfoam_probe_file(os.path.join(of_dir, "bottomWallP", "0", "p"))

    if heights_toe is None and pr is None and pl is None and pb is None:
        return None

    return {"heights_toe": heights_toe, "right": pr, "left": pl, "bottom": pb}

# =========================
# PLOTTING (9 figs)
# =========================
def plot_all(xc: pd.DataFrame, of: dict | None):
    os.makedirs(OUT_DIR, exist_ok=True)
    setup_matplotlib()

    # XCALibre arrays
    t_xc = xc["t_nd"].to_numpy()
    h1_xc = (xc["h_L1"].to_numpy() / H)
    h2_xc = (xc["h_L2"].to_numpy() / H)
    h3_xc = (xc["h_L3"].to_numpy() / H)
    toe_xc = xc["toe_x"].to_numpy()

    pR_xc = [xc["p_right_1"].to_numpy() * p_scale,
             xc["p_right_2"].to_numpy() * p_scale,
             xc["p_right_3"].to_numpy() * p_scale]
    pL_xc = [xc["p_left_1"].to_numpy() * p_scale,
             xc["p_left_2"].to_numpy() * p_scale,
             xc["p_left_3"].to_numpy() * p_scale]
    pB_xc = [xc["p_bottom_1"].to_numpy() * p_scale,
             xc["p_bottom_2"].to_numpy() * p_scale,
             xc["p_bottom_3"].to_numpy() * p_scale]

    # OpenFOAM arrays (optional)
    of_ht = of["heights_toe"] if (of and of.get("heights_toe")) else None
    of_R = of["right"] if (of and of.get("right")) else None
    of_L = of["left"] if (of and of.get("left")) else None
    of_B = of["bottom"] if (of and of.get("bottom")) else None

    # --- FIG 1: Height L1 ---
    fig, ax = plt.subplots()
    ax.set_title("Water Height at L1 (h/H)")
    ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
    ax.set_ylabel("Dimensionless Height, h/H")
    ax.set_xlim(0, 8); ax.set_ylim(0, 1)

    x, y = load_ref_xy(os.path.join(REF_DIR, "fig1_lobovsky.txt"))
    plot_experiment(ax, x, y, c_ref1, "Lobovsky")
    x, y = load_ref_xy(os.path.join(REF_DIR, "fig1_zhang.txt"))
    plot_experiment(ax, x, y, c_ref2, "Zhang")

    if of_ht is not None:
        plot_openfoam(ax, of_ht["t_nd"], of_ht["h_L1"] / H, label="OpenFOAM")
    plot_xcalibre(ax, t_xc, h1_xc, label="XCALibre")

    ax.legend(loc="upper right")
    savefig(os.path.join(OUT_DIR, "Height_L1.png"))

    # --- FIG 2: Height L2 ---
    fig, ax = plt.subplots()
    ax.set_title("Water Height at L2 (h/H)")
    ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
    ax.set_ylabel("Dimensionless Height, h/H")
    ax.set_xlim(0, 8); ax.set_ylim(0, 1)

    x, y = load_ref_xy(os.path.join(REF_DIR, "fig2_lobovsky.txt"))
    plot_experiment(ax, x, y, c_ref1, "Lobovsky")
    x, y = load_ref_xy(os.path.join(REF_DIR, "fig2_zhang.txt"))
    plot_experiment(ax, x, y, c_ref2, "Zhang")
    x, y = load_ref_xy(os.path.join(REF_DIR, "fig2_nguyen.txt"))
    plot_experiment(ax, x, y, c_ref3, "Nguyen")

    if of_ht is not None:
        plot_openfoam(ax, of_ht["t_nd"], of_ht["h_L2"] / H, label="OpenFOAM")
    plot_xcalibre(ax, t_xc, h2_xc, label="XCALibre")

    ax.legend(loc="upper left")
    savefig(os.path.join(OUT_DIR, "Height_L2.png"))

    # --- FIG 3: Height L3 ---
    fig, ax = plt.subplots()
    ax.set_title("Water Height at L3 (h/H)")
    ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
    ax.set_ylabel("Dimensionless Height, h/H")
    ax.set_xlim(0, 8); ax.set_ylim(0, 1)

    x, y = load_ref_xy(os.path.join(REF_DIR, "fig3_lobovsky.txt"))
    plot_experiment(ax, x, y, c_ref1, "Lobovsky")
    x, y = load_ref_xy(os.path.join(REF_DIR, "fig3_zhang.txt"))
    plot_experiment(ax, x, y, c_ref2, "Zhang")
    x, y = load_ref_xy(os.path.join(REF_DIR, "fig3_nguyen.txt"))
    plot_experiment(ax, x, y, c_ref3, "Nguyen")

    if of_ht is not None:
        plot_openfoam(ax, of_ht["t_nd"], of_ht["h_L3"] / H, label="OpenFOAM")
    plot_xcalibre(ax, t_xc, h3_xc, label="XCALibre")

    ax.legend(loc="upper left")
    savefig(os.path.join(OUT_DIR, "Height_L3.png"))

    # --- FIG 4-6: Right wall pressures ---
    def plot_right_probe(fig_name, title, ylim, ref_files, probe_idx):
        fig, ax = plt.subplots()
        ax.set_title(title)
        ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
        ax.set_ylabel("Dimensionless Pressure, P*")
        ax.set_xlim(0, 7); ax.set_ylim(*ylim)

        for (rf, col, lab) in ref_files:
            x, y = load_ref_xy(os.path.join(REF_DIR, rf))
            plot_experiment(ax, x, y, col, lab)

        if of_R is not None:
            series = [of_R["p1_nd"], of_R["p2_nd"], of_R["p3_nd"]][probe_idx]
            plot_openfoam(ax, of_R["t_nd"], series, label="OpenFOAM")

        plot_xcalibre(ax, t_xc, pR_xc[probe_idx], label="XCALibre")

        ax.legend(loc="upper right")
        savefig(os.path.join(OUT_DIR, fig_name))

    plot_right_probe(
        "Pressure_Right_Probe1.png",
        "Right Wall Pressure (Probe 1)",
        (0, 5),
        [("fig4_lobovsky.txt", c_ref1, "Lobovsky"),
         ("fig4_zhang.txt",   c_ref2, "Zhang"),
         ("fig4_nguyen.txt",  c_ref3, "Nguyen")],
        probe_idx=0
    )
    plot_right_probe(
        "Pressure_Right_Probe2.png",
        "Right Wall Pressure (Probe 2)",
        (0, 2.4),
        [("fig5_lobovsky.txt", c_ref1, "Lobovsky"),
         ("fig5_zhang.txt",   c_ref2, "Zhang"),
         ("fig5_nguyen.txt",  c_ref3, "Nguyen")],
        probe_idx=1
    )
    plot_right_probe(
        "Pressure_Right_Probe3.png",
        "Right Wall Pressure (Probe 3)",
        (0, 1.0),
        [("fig6_lobovsky.txt", c_ref1, "Lobovsky"),
         ("fig6_zhang.txt",   c_ref2, "Zhang"),
         ("fig6_nguyen.txt",  c_ref3, "Nguyen")],
        probe_idx=2
    )

    # --- FIG 7: Toe ---
    fig, ax = plt.subplots()
    ax.set_title("Water Front Position (X_{foot})")
    ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
    ax.set_ylabel("Position (m)")
    ax.set_xlim(0, 3); ax.set_ylim(0.4, 1.8)

    x, y = load_ref_xy(os.path.join(REF_DIR, "fig7_garoosi.txt"))
    plot_experiment(ax, x, y, c_ref1, "Garoosi")

    if of_ht is not None:
        plot_openfoam(ax, of_ht["t_nd"], of_ht["toe_x"], label="OpenFOAM")
    plot_xcalibre(ax, t_xc, toe_xc, label="XCALibre")

    ax.legend(loc="upper left")
    savefig(os.path.join(OUT_DIR, "Toe_Position.png"))

    # --- FIG 8: Left wall (all probes) ---
    fig, ax = plt.subplots()
    ax.set_title("Left Wall Pressure History")
    ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
    ax.set_ylabel("Dimensionless Pressure, P*")
    ax.set_xlim(0, 8); ax.set_ylim(0, 1)

    # Experiments (keep labels, they populate legend 1)
    for rf, col, lab in [("fig8_hl1.txt", c_ref_p1, "Garoosi P1"),
                         ("fig8_hl2.txt", c_ref_p2, "Garoosi P2"),
                         ("fig8_hl3.txt", c_ref_p3, "Garoosi P3")]:
        x, y = load_ref_xy(os.path.join(REF_DIR, rf))
        plot_experiment(ax, x, y, col, lab)

    # Models: show all curves, but no legend spam
    for i, col in enumerate(probe_cols_lr):
        if of_L is not None:
            series = [of_L["p1_nd"], of_L["p2_nd"], of_L["p3_nd"]][i]
            plot_openfoam_probe(ax, of_L["t_nd"], series, color=col)
        plot_xcalibre_probe(ax, t_xc, pL_xc[i], color=col)

    add_solver_probe_legends(ax, exp_loc="upper right", style_loc="lower left")
    savefig(os.path.join(OUT_DIR, "Pressure_Left_All.png"))

    # --- FIG 9: Bottom wall (all probes) ---
    fig, ax = plt.subplots()
    ax.set_title("Bottom Wall Pressure History")
    ax.set_xlabel("Dimensionless Time, t(g/H)^{0.5}")
    ax.set_ylabel("Dimensionless Pressure, P*")
    ax.set_xlim(0, 8); ax.set_ylim(0, 1.4)

    for rf, col, lab in [("fig9_l1.txt", c_ref_p1, "Garoosi P1"),
                         ("fig9_l2.txt", c_ref_p2, "Garoosi P2"),
                         ("fig9_l3.txt", c_ref_p3, "Garoosi P3")]:
        x, y = load_ref_xy(os.path.join(REF_DIR, rf))
        plot_experiment(ax, x, y, col, lab)

    for i, col in enumerate(probe_cols_lr):
        if of_B is not None:
            series = [of_B["p1_nd"], of_B["p2_nd"], of_B["p3_nd"]][i]
            plot_openfoam_probe(ax, of_B["t_nd"], series, color=col)
        plot_xcalibre_probe(ax, t_xc, pB_xc[i], color=col)

    add_solver_probe_legends(ax, exp_loc="upper right", style_loc="upper left")
    savefig(os.path.join(OUT_DIR, "Pressure_Bottom_All.png"))

# =========================
# MAIN
# =========================
def main():
    if not os.path.isdir(REF_DIR):
        raise SystemExit(f"Missing '{REF_DIR}/' directory.")

    vtk_files = list_vtk_files()
    if not vtk_files:
        raise SystemExit("No VTK files found matching time_*.vtk in current directory.")

    print(f"🚀 Found {len(vtk_files)} XCALibre VTK files.")
    print(f"🕒 XCALibre dt = {dt:.8f} s (time_{FINAL_STEP} -> {FINAL_TIME_S}s)")

    xc = extract_xcalibre_timeseries(vtk_files)
    xc.to_csv("processed_xcalibre_timeseries.csv", index=False)
    print("🧾 Wrote processed_xcalibre_timeseries.csv")

    of = None
    if os.path.isdir(OF_POST_DIR):
        of = extract_openfoam_all(OF_POST_DIR)
        if of:
            if of.get("heights_toe"):
                pd.DataFrame(of["heights_toe"]).to_csv("processed_openfoam_heights_toe.csv", index=False)
            for key in ["right", "left", "bottom"]:
                if of.get(key):
                    pd.DataFrame({
                        "t": of[key]["t"],
                        "t_nd": of[key]["t_nd"],
                        "p1": of[key]["p1"], "p2": of[key]["p2"], "p3": of[key]["p3"],
                        "p1_nd": of[key]["p1_nd"], "p2_nd": of[key]["p2_nd"], "p3_nd": of[key]["p3_nd"],
                    }).to_csv(f"processed_openfoam_{key}.csv", index=False)
            print("🧾 Wrote processed_openfoam_*.csv")
        else:
            print("⚠️ postProcessing/ found, but nothing usable inside. Plotting without OpenFOAM.")
    else:
        print("ℹ️ No postProcessing/ folder. Plotting Experiments + XCALibre only.")

    plot_all(xc, of)
    print(f"📈 Done. All graphs saved to '{OUT_DIR}/'.")

if __name__ == "__main__":
    main()