#!/usr/bin/env python3
# =============================================================================
# Transient diagnostics — is this run steadily converging, oscillating, or
# diverging?
# -----------------------------------------------------------------------------
# Walks the heated-run VTU series (+ times.csv) and builds time histories of
# T, p and |U| at a few fixed probe points (near-wall downstream — past the
# pseudo-critical line — and the outlet), plus global field metrics per
# snapshot (max|U|, max T, and the spatial std of T as a checkerboard/speckle
# indicator). Each history is classified:
#
#   STEADY        -> asymptotes to a constant. A steady solver is valid & faster.
#   LIMIT CYCLE   -> settles into a regular oscillation. Genuinely unsteady
#                    (physical instability); the transient is the answer.
#   DIVERGENCE    -> grows without bound / speckle worsens. Numerical, not
#                    physics (the ψ / pseudo-critical problem).
#   DEVELOPING    -> still in the initial transient; run longer to decide.
#
# Usage:
#     python test/tatsumoto/steady_state_analysis.py [path/to/case.toml]
#     python test/tatsumoto/steady_state_analysis.py --dir path/to/<case>_results
# =============================================================================

import os
import re
import sys
import argparse
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

try:
    import tomllib                      # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib             # fallback: pip install tomli

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CAMPAIGN   = os.path.dirname(SCRIPT_DIR)

# --- Resolve the results directory -------------------------------------------
ap = argparse.ArgumentParser(description="Steady-state / limit-cycle / divergence diagnostics.")
ap.add_argument("case", nargs="?", default=os.path.join(CAMPAIGN, "configs", "supercritical.toml"),
                help="path to case.toml (default: sibling case.toml)")
ap.add_argument("--dir", default=None, help="results dir override (skips case.toml lookup)")
ap.add_argument("--stride", type=int, default=1, help="use every Nth snapshot")
args = ap.parse_args()

if args.dir:
    RESULTS_DIR = os.path.abspath(args.dir)
    CASE_NAME = os.path.basename(RESULTS_DIR).replace("_results", "")
    D_HYD = None
else:
    with open(args.case, "rb") as f:
        CFG = tomllib.load(f)
    CASE_NAME = CFG["case"]["name"]
    D_HYD = float(CFG["flow"]["D_hyd"])
    RUN_DIR = CFG["case"].get("run_dir") or os.path.join(CAMPAIGN, "runs", CASE_NAME)
    RESULTS_DIR = os.path.join(RUN_DIR, "vtk")

TIMES_CSV = os.path.join(RESULTS_DIR, "times.csv")
OUT_PNG   = os.path.join(os.path.dirname(RESULTS_DIR), "steady_state_analysis.png")


# --- IO helpers --------------------------------------------------------------
def load_times(path):
    """Return {iteration: time} from the solver's times.csv sidecar."""
    if not os.path.isfile(path):
        sys.exit(f"times.csv not found: {path}\nRun heated.jl first.")
    it2t = {}
    with open(path) as f:
        next(f)  # header: iteration,time
        for line in f:
            line = line.strip()
            if not line:
                continue
            a, b = line.split(",")[:2]
            it2t[int(a)] = float(b)
    return it2t


def cell_field(grid, name):
    """Fetch a field as a cell-centred array (interpolating point->cell if needed)."""
    if name in grid.cell_data:
        return np.asarray(grid.cell_data[name])
    if name in grid.point_data:
        return np.asarray(grid.point_to_cell_data(name).cell_data[name])
    return None


def snapshot_centres(grid):
    return np.asarray(grid.cell_centers().points)


# --- Time-series classifier --------------------------------------------------
def _autocorr_peak(x):
    """Largest normalised autocorrelation over lags [2, n/2]; ~periodic if high."""
    x = x - x.mean()
    n = len(x)
    if n < 8 or np.allclose(x, 0):
        return 0.0
    ac = np.correlate(x, x, mode="full")[n - 1:]
    ac = ac / (ac[0] + 1e-30)
    lo = 2
    hi = max(lo + 1, n // 2)
    return float(np.max(ac[lo:hi])) if hi > lo else 0.0


def classify(t, x, *, steady_tol=0.02, osc_tol=0.05):
    """Classify a scalar time history. Returns (label, metrics dict)."""
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    n = len(x)
    if n < 6:
        return "DEVELOPING (too few snapshots)", {}
    if not np.all(np.isfinite(x)):
        return "DIVERGENCE (non-finite values)", {"finite": False}

    late = x[n // 2:]
    tl = t[n // 2:]
    mean_late = float(np.mean(late))
    denom = abs(mean_late) + 1e-30

    # late-window linear trend
    A = np.polyfit(tl, late, 1)
    slope = float(A[0])
    drift = slope * (tl[-1] - tl[0]) / denom          # fractional change across late window
    detr = late - np.polyval(A, tl)
    osc = (detr.max() - detr.min()) / denom            # detrended peak-to-peak (fractional)

    # whole-series magnitude growth (early fifth vs late fifth)
    k = max(3, n // 5)
    early_mag = float(np.mean(np.abs(x[:k])))
    late_mag = float(np.mean(np.abs(x[-k:])))
    growth = (late_mag - early_mag) / (abs(early_mag) + 1e-30)

    # is the oscillation amplitude itself growing? (divergence even w/o mean drift)
    half = len(detr) // 2
    amp_early = np.std(detr[:half]) if half >= 2 else 0.0
    amp_late = np.std(detr[half:]) if half >= 2 else 0.0
    amp_growth = (amp_late - amp_early) / (amp_early + 1e-30)

    period = _autocorr_peak(detr)

    m = dict(drift=drift, osc=osc, growth=growth, amp_growth=amp_growth,
             period=period, mean_late=mean_late, slope=slope)

    # --- decision tree -------------------------------------------------------
    if amp_growth > 1.0 and osc > osc_tol:
        return "DIVERGENCE (growing oscillation)", m
    if abs(drift) > 0.15 and growth > 0.3:
        return "DIVERGENCE (unbounded growth)", m
    if abs(drift) < steady_tol and osc < steady_tol:
        return "STEADY (asymptotes to constant)", m
    if osc > osc_tol and abs(drift) < 0.05 and period > 0.3:
        return "LIMIT CYCLE (regular oscillation)", m
    if osc > osc_tol and abs(drift) < 0.05:
        return "OSCILLATING (irregular, no net drift)", m
    return "DEVELOPING (still in transient)", m


# --- Probe geometry ----------------------------------------------------------
def build_probes(grid):
    """Fixed-coordinate probes; flow axis = z. Returns list of (name, point)."""
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    span_z = zmax - zmin
    # quarter-pipe radius from bounds if D_hyd unknown
    R = (D_HYD / 2.0) if D_HYD else max(xmax - xmin, ymax - ymin)
    rw = 0.92 * R / np.sqrt(2.0)        # near-wall offset along x=y diagonal
    rc = 0.05 * R / np.sqrt(2.0)        # near-axis (core)
    z_out = zmin + 0.95 * span_z
    z_mid = zmin + 0.50 * span_z
    return [
        ("wall_outlet", np.array([xmin + rw, ymin + rw, z_out])),   # near-wall, past T_c'
        ("wall_mid",    np.array([xmin + rw, ymin + rw, z_mid])),
        ("outlet_core", np.array([xmin + rc, ymin + rc, z_out])),
        ("core_mid",    np.array([xmin + rc, ymin + rc, z_mid])),
    ]


def nearest_cells(centres, probes):
    idx = {}
    for name, pt in probes:
        idx[name] = int(np.argmin(np.sum((centres - pt) ** 2, axis=1)))
    return idx


# --- Main --------------------------------------------------------------------
def main():
    it2t = load_times(TIMES_CSV)
    pat = re.compile(r"^time_(\d+)\.vtu$")
    files = []
    for fn in os.listdir(RESULTS_DIR):
        mobj = pat.match(fn)
        if not mobj:
            continue
        it = int(mobj.group(1))
        if it in it2t:
            files.append((it2t[it], os.path.join(RESULTS_DIR, fn)))
    files.sort()
    files = files[:: max(1, args.stride)]
    if len(files) < 4:
        sys.exit(f"Need >=4 snapshots to judge; found {len(files)} in {RESULTS_DIR}")

    print(f"Analysing {len(files)} snapshots from {RESULTS_DIR}")
    times = np.array([t for t, _ in files])

    # fixed topology -> probe cells from the first snapshot
    g0 = pv.read(files[0][1])
    centres = snapshot_centres(g0)
    R = (D_HYD / 2.0) if D_HYD else None
    if R is not None:
        rr = np.sqrt(centres[:, 0] ** 2 + centres[:, 1] ** 2)
        nearwall_mask = rr >= 0.85 * R
    else:
        nearwall_mask = np.ones(len(centres), bool)
    probes = build_probes(g0)
    pcell = nearest_cells(centres, probes)
    print("Probe cells (name -> nearest cell centre):")
    for name, _ in probes:
        print(f"  {name:12s} @ {np.round(centres[pcell[name]], 5)}")

    # time histories
    series = {name: dict(T=[], p=[], U=[]) for name, _ in probes}
    g_maxU, g_maxT, g_stdT_wall = [], [], []

    for _, path in files:
        g = pv.read(path)
        T = cell_field(g, "T")
        U = cell_field(g, "U")
        p = cell_field(g, "p_rgh")
        if p is None:
            p = cell_field(g, "p")
        Umag = np.linalg.norm(U, axis=1) if U is not None else None

        for name, _ in probes:
            ci = pcell[name]
            series[name]["T"].append(float(T[ci]) if T is not None else np.nan)
            series[name]["p"].append(float(p[ci]) if p is not None else np.nan)
            series[name]["U"].append(float(Umag[ci]) if Umag is not None else np.nan)

        g_maxU.append(float(np.nanmax(Umag)) if Umag is not None else np.nan)
        g_maxT.append(float(np.nanmax(T)) if T is not None else np.nan)
        g_stdT_wall.append(float(np.nanstd(T[nearwall_mask])) if T is not None else np.nan)

    # --- verdicts ------------------------------------------------------------
    print("\n" + "=" * 70)
    print("PER-PROBE VERDICTS")
    print("=" * 70)
    overall = []
    for name, _ in probes:
        for fld in ("T", "U", "p"):
            label, m = classify(times, series[name][fld])
            overall.append(label.split()[0])
            extra = (f"  drift={m.get('drift', 0):+.3f} osc={m.get('osc', 0):.3f} "
                     f"period={m.get('period', 0):.2f}") if m else ""
            print(f"  {name:12s} {fld:2s} : {label}{extra}")

    glabel, gm = classify(times, g_stdT_wall)
    print(f"\n  near-wall std(T) (speckle) : {glabel}")
    ulabel, um = classify(times, g_maxU)
    print(f"  global max|U|              : {ulabel}")

    # overall summary
    print("\n" + "=" * 70)
    if "DIVERGENCE" in overall or "DIVERGENCE" in glabel or "DIVERGENCE" in ulabel:
        verdict = ("DIVERGENCE detected -> numerical, not physics. The run is not "
                   "trustworthy past the blow-up; tighten the ψ/pressure coupling or dt.")
    elif all(v == "STEADY" for v in overall):
        verdict = ("STEADY everywhere -> a steady-state solver is valid and would be "
                   "much faster. The transient is just a path to the fixed point.")
    elif "LIMIT" in overall:
        verdict = ("LIMIT CYCLE present -> genuinely unsteady (physical oscillation). "
                   "A steady solver is inappropriate; report the transient/oscillation.")
    else:
        verdict = ("MIXED / still DEVELOPING -> run longer, then re-check. No clean "
                   "asymptote or limit cycle yet.")
    print("OVERALL:", verdict)
    print("=" * 70)

    # --- plot ----------------------------------------------------------------
    plt.rcParams.update({"figure.dpi": 120, "font.size": 9})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    fig.suptitle(f"Transient diagnostics — {CASE_NAME}", fontweight="bold")

    axT, axU, axP, axG = axes.flat
    for name, _ in probes:
        axT.plot(times, series[name]["T"], label=name, lw=1.2)
        axU.plot(times, series[name]["U"], label=name, lw=1.2)
        axP.plot(times, series[name]["p"], label=name, lw=1.2)
    axT.set(title="Probe T(t) [K]", xlabel="t [s]", ylabel="T")
    axU.set(title="Probe |U|(t) [m/s]", xlabel="t [s]", ylabel="|U|")
    axP.set(title="Probe p_rgh(t) [Pa]", xlabel="t [s]", ylabel="p_rgh")
    for ax in (axT, axU, axP):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)

    axG.plot(times, g_maxT, label="max T", color="crimson")
    axG2 = axG.twinx()
    axG2.plot(times, g_stdT_wall, label="near-wall std(T)", color="navy", ls="--")
    axG2.plot(times, g_maxU, label="max|U|", color="green", ls=":")
    axG.set(title="Global metrics", xlabel="t [s]", ylabel="max T [K]")
    axG2.set_ylabel("std(T) [K] / max|U| [m/s]")
    axG.grid(alpha=0.3)
    lines = axG.get_lines() + axG2.get_lines()
    axG.legend(lines, [ln.get_label() for ln in lines], fontsize=7, loc="upper left")

    fig.text(0.5, 0.005, verdict, ha="center", fontsize=8, wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(OUT_PNG)
    print(f"\nFigure -> {OUT_PNG}")


if __name__ == "__main__":
    main()
