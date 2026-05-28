#!/usr/bin/env python3
"""
plot_verification.py
====================
Generates the two figures that present the discretisation-error verification
results for the α-transport scheme:

    Fig. 1  alpha_profile_quad40.png
            α(x) at t = 100 s on quad40 for the three configurations
            (upwind / vanleer / compressed), with the analytical erfc
            reference for the upwind case overlaid.

    Fig. 2  convergence_alpha.png
            log-log L¹ and RMS profile error vs Δx on quad20–quad160 for
            all three configurations. Reference slopes from the modified-
            equation analysis are overlaid.

The α profiles are read from the JSON dump produced by
`stage1_2_dump_profiles.jl` (run that first). The convergence numbers are
hardcoded from the stage-3 output and live in CONFIG below — update them
if you re-run stage 3 with new meshes.

Style follows VOF_cases/rising_bubble.py.
"""

import json
import math
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# =============================================================================
# CONFIG — update if stage 3 is rerun with different meshes / numbers
# =============================================================================

OUT_DIR    = "graphs"
PROFILE_JS = "alpha_profiles_quad40.json"     # produced by Julia helper script

# Modified-equation parameters (matches stage 1)
U_FREESTREAM = 1.0
DX_QUAD40    = 25.0
T_END        = 100.0
COURANT      = U_FREESTREAM * 1.0 / DX_QUAD40                # = 0.04
D_NUM        = 0.5 * U_FREESTREAM * DX_QUAD40 * (1.0 - COURANT)
SIGMA_PRED   = math.sqrt(2.0 * D_NUM * T_END)                # 48.99 m
X_INTERFACE_END = 600.0

# Stage-3 results (Δx, L¹, RMS) for each configuration, geometric mesh series
# quad20 → quad40 → quad80 → quad160.
STAGE3 = {
    "upwind": {
        "dx":  [50.0, 25.0, 12.5, 6.25],
        "L1":  [53.0365, 38.2826, 27.3528, 19.4424],
        "RMS": [1.2324e-1, 1.0521e-1, 8.9203e-2, 7.5332e-2],
        "color": "#1f77b4",
        "label": "Pure upwind",
        "marker": "o",
    },
    "vanleer": {
        "dx":  [50.0, 25.0, 12.5, 6.25],
        "L1":  [19.4184, 10.1708, 5.1454, 2.5746],
        "RMS": [6.0040e-2, 4.4287e-2, 3.1648e-2, 2.2393e-2],
        "color": "#ff7f0e",
        "label": "van-Leer HO",
        "marker": "s",
    },
    "compressed": {
        "dx":  [50.0, 25.0, 12.5, 6.25],
        "L1":  [5.4495, 2.7247, 1.3624, 0.6812],
        "RMS": [1.7233e-2, 1.2185e-2, 8.6163e-3, 6.0927e-3],
        "color": "#2ca02c",
        "label": "Compressed",
        "marker": "^",
    },
}

# Predicted slopes (modified-equation σ-scaling for a step IC)
SLOPE_PRED = {
    "upwind":     {"L1": 0.5, "RMS": 0.25},
    "vanleer":    {"L1": 1.0, "RMS": 0.5},
    "compressed": {"L1": 1.0, "RMS": 0.5},
}


# =============================================================================
# STYLE
# =============================================================================

def _setup_matplotlib():
    plt.rcParams.update({
        "figure.figsize":      (7.2, 5.2),
        "font.size":           14,
        "axes.titlesize":      19,
        "axes.labelsize":      18,
        "xtick.labelsize":     15,
        "ytick.labelsize":     15,
        "font.family":         "serif",
        "font.serif":          ["Computer Modern Roman", "DejaVu Serif", "serif"],
        "mathtext.fontset":    "cm",
        "axes.grid":           True,
        "grid.alpha":          0.25,
        "grid.linewidth":      0.6,
        "axes.linewidth":      0.8,
        "xtick.direction":     "in",
        "ytick.direction":     "in",
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.top":           True,
        "ytick.right":         True,
        "legend.frameon":      True,
        "legend.framealpha":   0.9,
        "legend.edgecolor":    "0.7",
        "legend.fontsize":     12,
    })


# =============================================================================
# FIG 1 — α-profile overlay on quad40
# =============================================================================

def plot_alpha_profiles():
    if not os.path.isfile(PROFILE_JS):
        print(f"⚠️  {PROFILE_JS} not found — run stage1_2_dump_profiles.jl first.")
        return

    with open(PROFILE_JS, "r") as f:
        profiles = json.load(f)

    fig, ax = plt.subplots()
    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$\alpha$")
    ax.set_xlim(0.0, 1000.0)
    ax.set_ylim(-0.05, 1.05)

    # Analytical step IC — drawn first, behind everything, but visible.
    ax.axvline(X_INTERFACE_END, color="black", lw=1.5, ls=(0, (5, 3)), alpha=0.9,
               label=fr"analytical step $x={X_INTERFACE_END:.0f}$ m", zorder=3)

    # Analytical erfc reference for the upwind modified-equation prediction.
    # Use black so it stands clear of the blue upwind data series.
    x_fine = np.linspace(0.0, 1000.0, 1000)
    alpha_erfc = 0.5 * erfc((x_fine - X_INTERFACE_END) / (SIGMA_PRED * math.sqrt(2.0)))
    ax.plot(x_fine, alpha_erfc, color="black", lw=2.2, ls=":",
            label=fr"erfc, $\sigma_\mathrm{{pred}}={SIGMA_PRED:.1f}$ m", zorder=5)

    # Overlay simulation profiles
    plot_order = [("upwind", "o"), ("vanleer", "s"), ("compressed", "^")]
    for key, marker in plot_order:
        if key not in profiles:
            continue
        cfg = STAGE3[key]
        xs = np.asarray(profiles[key]["x"], dtype=float)
        αs = np.asarray(profiles[key]["alpha"], dtype=float)
        ax.plot(xs, αs, color=cfg["color"], lw=1.6, marker=marker, ms=4.5,
                mfc="white", mec=cfg["color"], label=cfg["label"], zorder=4)

    ax.legend(loc="lower left")
    plt.tight_layout()

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "alpha_profile_quad40.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"📈  Saved → {out_path}")


# =============================================================================
# FIG 2 — log-log convergence (one figure per metric)
# =============================================================================

def _fit_slope(dx, err):
    logx = np.log(dx)
    logy = np.log(err)
    p = np.polyfit(logx, logy, 1)
    return p[0], p[1]   # slope, intercept


def _plot_one_metric(metric, ylabel, fname):
    fig, ax = plt.subplots()

    for key in ("upwind", "vanleer", "compressed"):
        cfg  = STAGE3[key]
        dx   = np.asarray(cfg["dx"])
        err  = np.asarray(cfg[metric])
        slope, intercept = _fit_slope(dx, err)
        slope_pred = SLOPE_PRED[key][metric]

        ax.loglog(dx, err, color=cfg["color"], lw=1.6,
                  marker=cfg["marker"], ms=7, mfc="white", mec=cfg["color"],
                  label=fr"{cfg['label']}: $p={slope:.2f}$ "
                        fr"(pred {slope_pred:g})",
                  zorder=4)

    # Reference slope guides — anchor at top-right of data range.
    # upwind has slope ½, vanleer/compressed share slope 1, so plot one
    # guide per distinct predicted slope.
    dx_ref = np.array([6.0, 60.0])
    for key, ls in [("upwind", "--"), ("vanleer", ":")]:
        cfg = STAGE3[key]
        anchor_x = max(cfg["dx"])
        anchor_y = max(cfg[metric])
        slope_pred = SLOPE_PRED[key][metric]
        y_ref = anchor_y * (dx_ref / anchor_x) ** slope_pred
        ax.loglog(dx_ref, y_ref, color=cfg["color"], lw=0.9, ls=ls, alpha=0.4)

    ax.set_xlabel(r"$\Delta x$ [m]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(5.0, 60.0)
    ax.legend(loc="lower right", fontsize=11)

    plt.tight_layout()

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, fname)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"📈  Saved → {out_path}")


def plot_convergence():
    _plot_one_metric("L1",  r"$L^1$ error", "convergence_L1.png")
    _plot_one_metric("RMS", r"RMS error",   "convergence_RMS.png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    _setup_matplotlib()
    plot_alpha_profiles()
    plot_convergence()
    print("✅  Done.")


if __name__ == "__main__":
    main()
