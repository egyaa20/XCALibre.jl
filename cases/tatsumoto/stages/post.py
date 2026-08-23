#!/usr/bin/env python3
# =============================================================================
# Stage 4 — POSTPROCESS / VALIDATION
# -----------------------------------------------------------------------------
# Reads the heated-run VTU series (+ times.csv) from <case>_results/, extracts
# the peak heated-wall temperature per snapshot, and plots the boiling/heat-
# transfer curve  q vs ΔT_L  against the digitised Tatsumoto (2010) data
# selected by operating pressure.
#
# All parameters come from case.toml (same file the Julia stages read), so
# Q0/tau/pressure/T_in are never duplicated.
#
# Usage:
#     python test/tatsumoto/post.py [path/to/case.toml]
# =============================================================================

import os
import re
import sys
import math
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

try:
    import tomllib                      # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib             # fallback: pip install tomli

# --- Load configuration ------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CAMPAIGN   = os.path.dirname(SCRIPT_DIR)
CASE_FILE  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(CAMPAIGN, "configs", "supercritical.toml")
with open(CASE_FILE, "rb") as f:
    CFG = tomllib.load(f)

CASE_NAME     = CFG["case"]["name"]
Q0            = float(CFG["heating"]["Q0"])
TAU           = float(CFG["heating"]["tau"])
PRESSURE      = float(CFG["thermo"]["pressure"])
T_IN          = float(CFG["flow"]["T_in"])
D_HYD         = float(CFG["flow"]["D_hyd"])
DISCARD_UNTIL = float(CFG["run"]["discard_until"])
POST_STRIDE   = int(CFG["run"].get("post_stride", 1))   # 1 = use every snapshot

RUN_DIR     = CFG["case"].get("run_dir") or os.path.join(CAMPAIGN, "runs", CASE_NAME)
RESULTS_DIR = os.path.join(RUN_DIR, "vtk")
TIMES_CSV   = os.path.join(RESULTS_DIR, "times.csv")
# summary/ (not RUN_DIR directly): hpc/poller.sh publishes everything under
# <rundir>/summary/ to the pipeline-status branch after each job -- writing
# here means validation.png shows up there automatically, no SSH needed to
# go looking for it on scratch.
SUMMARY_DIR = os.path.join(RUN_DIR, "summary")
os.makedirs(SUMMARY_DIR, exist_ok=True)
OUT_PNG     = os.path.join(SUMMARY_DIR, "validation.png")

R_WALL       = D_HYD / 2.0          # pipe radius
NEARWALL_FAC = 0.9                  # cells with r >= 0.9R count as near-wall

# =============================================================================
# Reference data — Tatsumoto (2010), digitised  (ΔT_L [K], q [W/m^2])
# =============================================================================
TATSUMOTO_3p5MPa = np.array(sorted([   # Fig. 5, P=3.5 MPa, Tin=78 K (supercritical)
    (3.913743637129656,   6277.011260085071),
    (4.543673900348009,   7331.092117354604),
    (5.427097986443437,   8562.181810139109),
    (6.390802916802306,   9619.358517243929),
    (7.262917064184462,  10531.05775100815),
    (8.022704224879188,  12299.508678079948),
    (8.675042144138516,  12952.677228103266),
    (9.65090183828954,   14551.957697468608),
    (10.736536462273268, 15931.156597142464),
    (11.859705360152178, 17441.07258962322),
    (13.574252802742134, 19342.69412176486),
    (14.994288575060652, 21451.651781482473),
    (16.68100336211878,  24731.952383696815),
    (19.22872890408915,  27075.98623972674),
    (21.39177793549197,  31622.78050286181),
    (24.65899021591013,  35527.27795575324),
    (27.04573863833884,  39400.86063777456),
    (30.95570111641362,  43135.17768071612),
    (35.18003806393228,  49092.072054095224),
    (40.55316098689205,  49731.236145771974),
    (47.75429892144337,  52372.25198801486),
    (57.03914171761611,  61166.94687094465),
    (66.69206044648602,  67836.08757096462),
    (75.7931217983432,   70520.38600012632),
    (84.31912679381838,  74265.3891117446),
    (95.14701147639234,  84521.31252096554),
    (108.13108893651477, 94957.24649657779),
    (118.59711790505723, 108070.68007401533),
    (131.93816516258795, 121414.24888908741),
    (152.08936506635288, 134652.2265567523),
    (166.81012406588786, 165615.58204306275),
    (182.9557081702406,  193426.89349300967),
    (200.66402647631259, 217309.4587016483),
]))

TATSUMOTO_1MPa = np.array(sorted([    # Fig. 4(a), P=1.0 MPa, Tin=78 K, Tsat=103.7 K
    ( 1.22835650807441,    5458.7004669872495),
    ( 1.3994886115958445,  6282.346922128502),
    ( 1.6105400862839656,  7230.276810407023),
    ( 1.8721085412031135,  8321.227517279649),
    ( 2.154434002892006,   9680.89541297105),
    ( 2.4179210564387965, 11021.822563360989),
    ( 2.768635315628035,  12822.739826516095),
    ( 3.2021852858959585, 14598.870451500688),
    ( 3.611886374579019,  16442.292787632265),
    ( 3.933428146875538,  18319.378504970748),
    ( 4.348547434070129,  19974.19966691078),
    ( 4.807476857402031,  21544.343938842758),
    ( 5.422563643671501,  22988.06392709914),
    ( 6.085743775565414,  25890.80291970923),
    ( 6.694344189386675,  28846.549846170037),
    ( 7.475465834277752,  33559.97672025372),
    ( 8.264398524210804,  38623.75389056604),
    ( 9.36866024798557,   41659.913884112546),
    (10.514448429960979,  46415.89471643547),
    (11.979297477812988,  53419.46775118258),
    (13.310150206467748,  61479.79161563335),
    (14.937983045231162,  65599.62099944215),
    (16.267804240168132,  70756.31668041172),
    (19.004963495305677,  84116.79134323457),
]))

# Subcooled FLOW BOILING, P=0.5 MPa, Tin=77.84 K, ΔTsub=16.2 K (Fig. 3b).
# NOT sorted — this is a measured trajectory with a DNB fold near ΔT_L~17-19 K
# and a post-DNB decline past the q_DNB peak (~401 kW/m^2 at ΔT_L~40 K). Sorting
# would destroy the fold/decline, so the points are kept in measurement order.
TATSUMOTO_0p5MPa = np.array([   # (ΔT_L [K], q [W/m^2])
    ( 2.72,  27826), ( 2.93,  29573), ( 3.15,  31050), ( 3.48,  34229),
    ( 3.84,  36379), ( 4.26,  40594), ( 4.79,  45299), ( 5.28,  50548),
    ( 5.83,  55048), ( 6.47,  62180), ( 7.27,  63714), ( 7.92,  70236),
    ( 8.84,  77426), ( 9.76,  84319), (11.24,  96411), (12.78, 111588),
    (14.02, 123012), (15.66, 138950), (17.49, 147677), (18.71, 162795),
    (19.65, 181660), (18.59, 205196), (17.70, 186141), (17.17, 207711),
    (16.96, 231782), (18.25, 249359), (20.26, 265021), (22.08, 278256),
    (23.77, 299358), (26.06, 322060), (29.11, 342288), (32.91, 363786),
    (36.09, 372759), (40.06, 401028), (45.29, 396172), (48.75, 368246),
    (54.45, 359381), (60.81, 350730), (69.61, 342288), (80.66, 346483),
    (91.20, 338143),
])

def reference_curve(pressure_pa):
    """Pick the Tatsumoto dataset + label closest to the operating pressure."""
    if pressure_pa >= 3.0e6:
        return TATSUMOTO_3p5MPa, "Tatsumoto et al., 2010  (3.5 MPa)", 126.8, r"$T_c'=126.8$ K"
    elif pressure_pa <= 0.75e6:   # 0.5 MPa subcooled boiling case
        return TATSUMOTO_0p5MPa, "Tatsumoto et al., 2010  (0.5 MPa)", 94.02, r"$T_{sat}=94.02$ K"
    else:
        return TATSUMOTO_1MPa, "Tatsumoto et al., 2010  (1.0 MPa)", 103.7, r"$T_{sat}=103.7$ K"

# =============================================================================
# VTU EXTRACTION
# =============================================================================
def load_times(path):
    """Return {iteration: time} from the solver's times.csv sidecar."""
    mapping = {}
    if not os.path.isfile(path):
        sys.exit(f"times.csv not found: {path}\nRun heated.jl first.")
    with open(path) as f:
        next(f)  # header: iteration,time
        for line in f:
            it, t = line.split(",")
            mapping[int(it)] = float(t)
    return mapping


# The near-wall mask is computed ONCE — the mesh geometry is identical across
# every snapshot, so there is no need to recompute cell centres per file.
_MASK = None        # boolean array selecting near-wall cells/points


def get_temperature(grid):
    """T field, whether XCALibre wrote it as cell or point data."""
    if "T" in grid.cell_data:
        return np.asarray(grid.cell_data["T"]), grid.cell_centers().points
    if "T" in grid.point_data:
        return np.asarray(grid.point_data["T"]), grid.points
    raise KeyError("No 'T' field in VTU")


def build_mask(first_vtu):
    """Read ONE snapshot (pyvista) to build the near-wall mask. Reused for all."""
    global _MASK
    _, pts = get_temperature(pv.read(first_vtu))
    r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)    # tube axis along z
    mask = r >= NEARWALL_FAC * R_WALL
    _MASK = mask if np.any(mask) else np.ones(len(r), dtype=bool)


def read_T_fast(vtu_path):
    """Stream-parse ONLY the `T` DataArray, skipping points/connectivity and the
    other 9 fields. Far cheaper than a full pyvista parse when we just need T."""
    vals = []
    capturing = False
    with open(vtu_path) as f:
        for line in f:
            if not capturing:
                if "DataArray" in line and 'Name="T"' in line:
                    capturing = True            # values start on the next line
                continue
            if "</DataArray>" in line:
                break
            vals.extend(line.split())
    return np.asarray(vals, dtype=float)


def wall_peak_temperature(vtu_path):
    """Peak near-wall temperature = heated-wall surface temp (fast path)."""
    T = read_T_fast(vtu_path)
    if _MASK is not None and len(T) == len(_MASK):
        return float(np.nanmax(T[_MASK]))
    # Streaming count mismatch (unexpected layout) -> robust pyvista fallback.
    Tarr, pts = get_temperature(pv.read(vtu_path))
    r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    nw = r >= NEARWALL_FAC * R_WALL
    if not np.any(nw):
        nw = np.ones_like(Tarr, dtype=bool)
    return float(np.nanmax(Tarr[nw]))


def collect_points():
    it2time = load_times(TIMES_CSV)
    pat = re.compile(r"^time_(\d+)\.vtu$")

    # 1) gather the kept snapshots (filtered by discard time), sorted by time
    kept = []
    for fname in os.listdir(RESULTS_DIR):
        m = pat.match(fname)
        if not m:
            continue
        it = int(m.group(1))
        if it not in it2time or it2time[it] < DISCARD_UNTIL:
            continue
        kept.append((it2time[it], fname))
    kept.sort()

    # 2) optional stride for quick previews (POST_STRIDE = 1 keeps all)
    if POST_STRIDE > 1:
        kept = kept[::POST_STRIDE]
    if not kept:
        return []

    # 3) build the near-wall mask ONCE, then fast-read T for every snapshot
    build_mask(os.path.join(RESULTS_DIR, kept[0][1]))
    points = []
    for t, fname in kept:
        T_w = wall_peak_temperature(os.path.join(RESULTS_DIR, fname))
        points.append((t, T_w - T_IN, Q0 * math.exp(t / TAU), T_w))
    return points

# =============================================================================
# PLOT  (style mirrors the heat-flux poster figure)
# =============================================================================
def setup_matplotlib():
    plt.rcParams.update({
        "figure.figsize": (7.4, 5.4), "figure.dpi": 120,
        "font.size": 16, "axes.labelsize": 18,
        "xtick.labelsize": 14, "ytick.labelsize": 14,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif", "serif"],
        "mathtext.fontset": "cm",
        "axes.grid": True, "grid.alpha": 0.30, "grid.linewidth": 0.5,
        "axes.linewidth": 1.0,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 6, "ytick.major.size": 6,
        "xtick.minor.size": 3, "ytick.minor.size": 3,
        "legend.frameon": False, "legend.fontsize": 14,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    })


def main():
    points = collect_points()
    if not points:
        sys.exit("No VTU points found past discard_until — nothing to plot.")
    print(f"Collected {len(points)} points (t >= {DISCARD_UNTIL} s)")
    for t, dT, q, T_w in points:
        print(f"  t={t:7.3f}s  T_w={T_w:8.3f}K  dT={dT:8.3f}K  q={q:.3e} W/m^2")

    dT = np.array([p[1] for p in points])
    q  = np.array([p[2] for p in points])

    ref, ref_label, t_crit, crit_label = reference_curve(PRESSURE)

    setup_matplotlib()
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$\Delta T_L = T_w - T_\mathrm{in}$  [K]")
    ax.set_ylabel(r"Heat flux  $q_w$  [W/m$^2$]")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1.0e-1, 1.0e3); ax.set_ylim(1.0e2, 1.0e6)

    # Reference: neutral black, open squares + thin line.
    ax.plot(ref[:, 0], ref[:, 1], color="0.15", lw=1.6, ls="-",
            marker="s", ms=6, mfc="white", mec="0.15", mew=1.2,
            label=ref_label, zorder=4)
    # Pseudo-critical / saturation marker.
    ax.axvline(t_crit - T_IN, color="0.55", lw=0.9, ls="--", alpha=0.8,
               label=crit_label)
    # XCALibre: saturated navy, filled circles.
    ax.plot(dT, q, color="#1f4e89", lw=2.0, ls="-",
            marker="o", ms=5.0, mfc="#1f4e89", mec="#1f4e89",
            label="XCALibre", zorder=5)

    ax.legend(loc="upper left", borderpad=0.5)
    for spine in ax.spines.values():
        spine.set_color("0.2")
    ax.tick_params(colors="0.2")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved -> {OUT_PNG}")


if __name__ == "__main__":
    main()
