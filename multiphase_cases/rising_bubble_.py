#!/usr/bin/env python3
"""
plot_bubble_front.py
====================
Tracks the leading edge of a rising bubble (alpha = 0 phase) from XCALibre
VTK output and plots its position against non-dimensional time, optionally
compared against reference datasets.

  T* = t * sqrt(g / R)
  Y  = bubble_front_coord  [m]  (raw position, not non-dimensionalised)

Usage
-----
  1. Edit the CONFIG section below.
  2. Run from the directory that contains your time_*.vtk files:
         python plot_bubble_front.py
"""

import os
import re
import math
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# =============================================================================
# USER CONFIG
# =============================================================================

# ---- VTK file naming --------------------------------------------------------
VTK_PREFIX = "time_"          # filename prefix
VTK_EXT    = ".vtk"           # filename extension

# ---- Time mapping -----------------------------------------------------------
# Map the final VTK step to a physical time, so the script infers dt for every
# intermediate step.
FINAL_STEP   = 4800          # integer step number of the last VTK file
FINAL_TIME_S = 0.96            # physical time (seconds) that FINAL_STEP represents

# ---- Physics / non-dimensionalisation ---------------------------------------
g_mag   = 9.80                # gravitational acceleration magnitude  [m/s²]
R       = 0.25                # reference length (e.g. bubble radius) [m]
T_SCALE = math.sqrt(g_mag / R)   # multiply t[s]  → T* (dimensionless)

# ---- Alpha field settings ---------------------------------------------------
FIELD_ALPHA      = "alpha"    # name of the volume-fraction field in the VTK
BUBBLE_THRESHOLD = 0.5        # alpha < BUBBLE_THRESHOLD  →  bubble phase

# ---- Bubble front tracking --------------------------------------------------
# The front is found by taking ALL mesh points where alpha < BUBBLE_THRESHOLD
# and returning the maximum value of FRONT_COORD_IDX among them.
# No scan line needed — the whole mesh is queried directly.
#
# FRONT_COORD_IDX: 0 = X, 1 = Y, 2 = Z
#   → vertically rising bubble   → 1 (Y)
#   → horizontally moving bubble → 0 (X)
FRONT_COORD_IDX = 1

# ---- Reference datasets -----------------------------------------------------
# Each entry is a dict with:
#   "label"  – legend label
#   "color"  – line/marker color
#   "data"   – list of (T*, Y*) tuples  ← already non-dimensional
#              OR set to None and provide "file" pointing to a two-column
#              CSV/TXT (comma- or whitespace-separated, no header) that
#              contains non-dimensional values directly.
#
# Set REF_DATASETS = [] to disable all reference overlays.

REF_DATASETS = [
    {
        "label": "Garoosi et al.",
        "color": "#1f77b4",         # blue
        "file":  None,              # set to e.g. "ref_data/garoosi.txt" to load from file
        "data": [
            (0.031250117346637234, 0.6227272727272726),
            (0.13281255168839956,  0.6227272727272726),
            (0.24218736635521815,  0.6254545731977982),
            (0.3671875377185617,   0.6327272588556461),
            (0.5078127677552238,   0.6436363913796164),
            (0.6250002607703048,   0.6527273004705255),
            (0.7265623970888615,   0.6599999861283735),
            (0.8437498901039425,   0.6709091186523436),
            (0.9921875004656612,   0.6845454822887073),
            (1.1328127305023232,   0.7),
            (1.2734373644925734,   0.7163636641068891),
            (1.4296873551793485,   0.7309091047807172),
            (1.5703125852160105,   0.7463636224920098),
            (1.7656250735744787,   0.7690908952192826),
            (1.9921876792795847,   0.7918181679465554),
            (2.2578121866099727,   0.8209090839732777),
            (2.4531246749684414,   0.8436363567005504),
            (2.6406247830018534,   0.8618181748823686),
            (2.851562628056847,    0.8827272935347124),
            (3.0234373793937337,   0.9027272657914596),
            (3.2343746284023145,   0.9281818389892578),
            (3.4609372341074205,   0.9545454545454546),
            (3.6874998398125256,   0.9790909160267223),
            (3.914062445517633,    1.0045454545454546),
            (4.179687548894432,    1.0400000138716265),
            (4.453125032596288,    1.0690909298983486),
            (4.7343748966231995,   1.1009090943769975),
            (4.968749882653364,    1.1281818216497248),
            (5.2187496293336375,   1.15363637750799),
            (5.484374732710437,    1.1827272761951795),
            (5.742186859715768,    1.2136363809758968),
            (5.953125300817173,    1.2363636537031695),
        ],
    },
    {
        "label": "Zhang et al.",
        "color": "#ff7f0e",         # orange
        "file":  None,
        "data": [
            (0.49999979138375533,  0.6336363358931105),
            (0.757812514435499,    0.659090909090909),
            (1.015625237487242,    0.6854545593261718),
            (1.2265624864958233,   0.7127272865988991),
            (1.499999970197679,    0.7454545454545454),
            (1.7499997168779537,   0.7718181956898081),
            (1.9921876792795847,   0.8027272657914595),
            (2.203124928288166,    0.8254545385187322),
            (2.4453122946433847,   0.8509091117165305),
            (2.6250000223517396,   0.8709090839732777),
            (2.828124891035265,    0.8954545454545454),
            (3.1250001117587014,   0.9290909160267222),
            (3.3671874781139204,   0.9636363636363636),
            (3.6874998398125256,   0.998181811246005),
            (3.914062445517633,    1.0254545385187324),
            (4.2031246898696,      1.057272720336914),
            (4.437499675899764,    1.083636370572177),
            (4.742187276948256,    1.1181818355213513),
            (5.0,                  1.1472727342085405),
            (5.304687005002082,    1.1781818216497248),
            (5.562500324100236,    1.2054545489224522),
            (5.835937807802092,    1.2418181766163219),
        ],
    },
    # ── Add more datasets here following the same pattern ──────────────────
    # {
    #     "label": "My experiment",
    #     "color": "#d62728",
    #     "file":  "ref_data/my_exp.txt",   # two-column file: T*  Y*
    #     "data":  None,
    # },
]

# ---- Output -----------------------------------------------------------------
OUT_DIR     = "graphs"
GRAPH_FNAME = "Bubble_Front_Position.png"
CSV_FNAME   = "processed_bubble_front.csv"

PLOT_XLIM = (0.0, 6.0)    # (T*_min, T*_max)
PLOT_YLIM = (0.6, 1.3)    # (Y_min,  Y_max)  [m]

# =============================================================================
# HELPERS
# =============================================================================

_dt = FINAL_TIME_S / FINAL_STEP    # seconds per step


def _step_number(filename: str) -> int:
    m = re.search(rf"{re.escape(VTK_PREFIX)}(\d+)", filename)
    return int(m.group(1)) if m else -1


def _list_vtk_files() -> list[str]:
    files = [
        f for f in os.listdir(".")
        if f.startswith(VTK_PREFIX) and f.endswith(VTK_EXT) and _step_number(f) >= 0
    ]
    files.sort(key=_step_number)
    return files


def _read_mesh(filepath: str) -> pv.DataSet:
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


def _to_point_data(mesh: pv.DataSet, field: str) -> pv.DataSet:
    if field in mesh.point_data:
        return mesh
    if field in mesh.cell_data:
        return mesh.cell_data_to_point_data()
    raise KeyError(f"Field '{field}' not found in point_data or cell_data.")


def _bubble_front(mesh: pv.DataSet) -> float:
    """
    Return the maximum coordinate (FRONT_COORD_IDX) across ALL mesh points
    where alpha < BUBBLE_THRESHOLD.  Returns NaN if no bubble phase found.
    """
    coords = np.asarray(mesh.points)[:, FRONT_COORD_IDX]
    alpha  = np.asarray(mesh.point_data[FIELD_ALPHA], dtype=float)

    mask = alpha < BUBBLE_THRESHOLD
    if not np.any(mask):
        return float("nan")

    return float(np.nanmax(coords[mask]))


def _load_ref_data(ds: dict) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (T*, Y*) arrays for a reference dataset dict.
    Prefers inline `data` list; falls back to `file` if data is None.
    """
    if ds.get("data") is not None:
        arr = np.array(ds["data"], dtype=float)
        return arr[:, 0], arr[:, 1]

    path = ds.get("file")
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(
            f"Reference dataset '{ds['label']}': no inline data and "
            f"file not found: {path!r}"
        )

    xs, ys = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in re.split(r"[,\s]+", line)]
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(xs, dtype=float), np.array(ys, dtype=float)


# =============================================================================
# EXTRACTION LOOP
# =============================================================================

def extract_timeseries(vtk_files: list[str]) -> pd.DataFrame:
    records = []
    for f in vtk_files:
        step = _step_number(f)
        t    = step * _dt

        mesh = _read_mesh(f)
        mesh = _to_point_data(mesh, FIELD_ALPHA)

        front = _bubble_front(mesh)

        records.append({
            "step":  step,
            "t":     t,
            "t_nd":  t * T_SCALE,
            "front": front,
        })
        print(f"✅  {f}   step={step:>7d}   t={t:.5f}s   "
              f"front={'nan' if math.isnan(front) else f'{front:.5f} m'}")

    return pd.DataFrame(records).sort_values("step").reset_index(drop=True)


# =============================================================================
# PLOT STYLE HELPERS
# =============================================================================

_HALO = [pe.Stroke(linewidth=6, foreground="white"), pe.Normal()]


def _plot_ref(ax, x, y, color, label):
    """Thin faint connecting line + open circle markers."""
    ax.plot(x, y, color=color, lw=1.2, alpha=0.35, zorder=3, label="_nolegend_")
    ax.plot(
        x, y, linestyle="none", marker="o", ms=5,
        mfc="white", mec=color, mew=1.4,
        color=color, label=label, zorder=5,
    )


def _plot_xcalibre(ax, x, y, label="XCALibre"):
    (ln,) = ax.plot(x, y, color="black", lw=3.2, ls="-", label=label, zorder=4)
    ln.set_path_effects(_HALO)


def _setup_matplotlib():
    plt.rcParams.update({
        "figure.figsize":    (10, 7),
        "font.size":         14,
        "axes.grid":         True,
        "grid.alpha":        0.30,
        "legend.frameon":    True,
        "legend.framealpha": 0.92,
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Arial", "DejaVu Sans"],
    })


# =============================================================================
# RMS ERROR
# =============================================================================
"""
Methodology — RMS error vs Garoosi et al.
==========================================
The Garoosi reference points are irregularly spaced in T*, while XCALibre
produces values at every simulation timestep.  To compare them on a common
basis the following procedure is used:

  1. Interpolation onto reference abscissae
     For each Garoosi point (T*_i, Y_ref_i) that falls within a sector, the
     XCALibre front position is linearly interpolated to T*_i using
     numpy.interp.  This avoids re-sampling the (sparser) reference data and
     means the RMS is evaluated exactly where experimental observations exist.

  2. Point-wise relative error
     e_i = ( Y_XCALibre(T*_i) - Y_Garoosi_i ) / Y_Garoosi_i  ×  100  [%]

  3. Sector RMS
     RMS = sqrt( mean( e_i^2 ) )   [%]   over all Garoosi points in the sector.

  4. Three sectors reported
     T* ∈ [0, 2],  [2, 4],  [4, 6]
     If a sector contains no Garoosi points the RMS is reported as NaN.
"""

RMS_SECTORS = [(0, 2), (2, 4), (4, 6)]
RMS_REF_LABEL = "Garoosi et al."       # must match a label in REF_DATASETS


def compute_rms(df: pd.DataFrame) -> dict:
    """
    Interpolate XCALibre onto Garoosi T* points and compute sector RMS.
    Returns a dict  { (t_lo, t_hi): rms_value }.
    """
    # Find Garoosi dataset
    garoosi = next((d for d in REF_DATASETS if d["label"] == RMS_REF_LABEL), None)
    if garoosi is None:
        print(f"⚠️  RMS: dataset '{RMS_REF_LABEL}' not found in REF_DATASETS — skipping.")
        return {}

    t_ref, y_ref = _load_ref_data(garoosi)

    t_xc = df["t_nd"].to_numpy()
    y_xc = df["front"].to_numpy()

    # Remove NaN rows from XCALibre (bubble not yet present)
    valid = ~np.isnan(y_xc)
    t_xc, y_xc = t_xc[valid], y_xc[valid]

    results = {}
    for (t_lo, t_hi) in RMS_SECTORS:
        mask = (t_ref >= t_lo) & (t_ref <= t_hi)
        if not np.any(mask):
            results[(t_lo, t_hi)] = float("nan")
            continue

        t_pts = t_ref[mask]
        y_pts = y_ref[mask]

        # Only interpolate within XCALibre's covered range
        in_range = (t_pts >= t_xc.min()) & (t_pts <= t_xc.max())
        if not np.any(in_range):
            results[(t_lo, t_hi)] = float("nan")
            continue

        t_pts  = t_pts[in_range]
        y_pts  = y_pts[in_range]

        y_interp       = np.interp(t_pts, t_xc, y_xc)
        errors_pct     = (y_interp - y_pts) / y_pts * 100.0   # % of reference
        results[(t_lo, t_hi)] = float(np.sqrt(np.mean(errors_pct ** 2)))

    return results


def print_rms(rms: dict):
    print("\n" + "=" * 52)
    print(f"  RMS error vs {RMS_REF_LABEL}")
    print(f"  (linear interp of XCALibre onto ref. abscissae)")
    print("=" * 52)
    print(f"  {'Sector (T*)':<20}  {'RMS [%]':>12}  {'N pts':>6}")
    print("-" * 52)

    garoosi = next((d for d in REF_DATASETS if d["label"] == RMS_REF_LABEL), None)
    t_ref, _ = _load_ref_data(garoosi) if garoosi else (np.array([]), np.array([]))

    for (t_lo, t_hi), rms_val in rms.items():
        n = int(np.sum((t_ref >= t_lo) & (t_ref <= t_hi)))
        val_str = f"{rms_val:.4f} %" if not math.isnan(rms_val) else "      NaN"
        print(f"  T* ∈ [{t_lo}, {t_hi}]        {val_str:>12}  {n:>6}")

    print("=" * 52 + "\n")


# =============================================================================
# PLOT
# =============================================================================

def plot_bubble_front(df: pd.DataFrame, rms: dict):
    os.makedirs(OUT_DIR, exist_ok=True)
    _setup_matplotlib()

    coord_labels = {0: "X", 1: "Y", 2: "Z"}
    coord_name   = coord_labels.get(FRONT_COORD_IDX, str(FRONT_COORD_IDX))

    fig, ax = plt.subplots()
    ax.set_title("Bubble Front Position")
    ax.set_xlabel(r"Dimensionless Time,  $T^* = t\,\sqrt{g/R}$")
    ax.set_ylabel(f"Bubble Front Position,  {coord_name} [m]")
    ax.set_xlim(*PLOT_XLIM)
    ax.set_ylim(*PLOT_YLIM)

    # -- Reference datasets ---------------------------------------------------
    for ds in REF_DATASETS:
        try:
            x_ref, y_ref = _load_ref_data(ds)
            _plot_ref(ax, x_ref, y_ref, color=ds["color"], label=ds["label"])
        except FileNotFoundError as e:
            print(f"⚠️  Skipping reference '{ds['label']}': {e}")

    # -- XCALibre -------------------------------------------------------------
    _plot_xcalibre(ax, df["t_nd"].to_numpy(), df["front"].to_numpy())

    # -- RMS annotations: shaded sectors + text ------------------------------
    sector_colors = ["#4daf4a", "#984ea3", "#e41a1c"]   # green / purple / red
    y_lo, y_hi = PLOT_YLIM
    for i, ((t_lo, t_hi), rms_val) in enumerate(rms.items()):
        col = sector_colors[i % len(sector_colors)]
        ax.axvspan(t_lo, t_hi, alpha=0.06, color=col, zorder=0)
        ax.axvline(t_lo, color=col, lw=0.8, ls=":", alpha=0.5, zorder=1)
        ax.axvline(t_hi, color=col, lw=0.8, ls=":", alpha=0.5, zorder=1)
        label = "nan" if math.isnan(rms_val) else f"{rms_val:.2f} %"
        ax.text(
            (t_lo + t_hi) / 2, y_lo + 0.97 * (y_hi - y_lo),
            f"RMS = {label}",
            ha="center", va="top", fontsize=10, color=col,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=col, alpha=0.85),
        )

    ax.legend(loc="upper left")
    plt.tight_layout()

    out_path = os.path.join(OUT_DIR, GRAPH_FNAME)
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"📈  Graph saved → {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    vtk_files = _list_vtk_files()
    if not vtk_files:
        raise SystemExit(
            f"No VTK files matching '{VTK_PREFIX}*{VTK_EXT}' found in the "
            "current directory."
        )

    print(f"🚀  Found {len(vtk_files)} XCALibre VTK file(s).")
    print(f"🕒  dt = {_dt:.8f} s  (step {FINAL_STEP} → {FINAL_TIME_S} s)")
    print(f"📐  T_SCALE = sqrt({g_mag}/{R}) = {T_SCALE:.6f}")
    print(f"📐  Tracking coord index {FRONT_COORD_IDX} across full mesh")
    print(f"📚  Reference datasets: "
          f"{len(REF_DATASETS)} "
          f"({', '.join(d['label'] for d in REF_DATASETS) or 'none'})\n")

    df = extract_timeseries(vtk_files)
    df.to_csv(CSV_FNAME, index=False)
    print(f"\n🧾  Data saved → {CSV_FNAME}")

    rms = compute_rms(df)
    print_rms(rms)

    plot_bubble_front(df, rms)
    print("✅  Done.")


if __name__ == "__main__":
    main()