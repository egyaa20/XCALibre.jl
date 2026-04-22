# ============================================================
# Stage 4 — 2D Rayleigh-Taylor instability (variable-density mixing)
# ============================================================
#
# Purpose
# -------
# Force the mixture solver to behave like a variable-density single-fluid
# flow and verify that buoyancy-driven instability grows at the classical
# inviscid rate and produces the canonical mushroom structure.
#
# Trick
# -----
# Use a tiny "particle" diameter (d = 1e-6 m) so the Stokes/SN drift velocity
# is O(1e-9) m/s — effectively zero slip. The mixture then reduces to a
# buoyant miscible fluid and the α transport tests the advective path + MULES
# boundedness in a stiff, lateral-instability regime.
#
# Setup
# -----
# * 2D rectangle 0.25 m (W) × 1.0 m (H), mesh `RTI_coarse.unv` (or medium).
# * Dense phase on top, light on bottom, separated by a smoothed tanh interface
#   with a single-mode cosine perturbation at y = H/2:
#       y_iface(x) = H/2 + A·cos(2π x/W),   A = 0.01·W
#       α_d(x,y)   = 0.5·(1 + tanh((y − y_iface(x))/(2Δy)))
#   Because the solver tracks α_c (continuous = light), we initialise
#       α_c = 1 − α_d.
# * All walls no-slip for U, zero-gradient for α, zero-gradient for p_rgh.
# * Laminar. Hydrostatic p_rgh anchored with pref=0.
# * Fixed dt (adaptive blows up in closed box; see Stages 2–3).
#
# Hardcoded constants to change in Solvers_2_MULTIPHASE.jl
# --------------------------------------------------------
#     diameter = 1.0e-6     # <-- make slip negligible for RTI
#     standard_manninen = true
#
# Phases
# ------
#     Phase(rho=1000.0, mu=1.0e-3)   # light (tracked as α_c)
#     Phase(rho=2000.0, mu=1.0e-3)   # heavy (dispersed)
#
# Analytical targets
# ------------------
# Atwood number  A_t = (ρ_d − ρ_c)/(ρ_d + ρ_c) = 1000/3000 = 0.333
# Wavenumber     k   = 2π/W = 25.13 rad/m
# Inviscid linear growth rate
#     σ = √(A_t · g · k) ≈ 9.07 1/s
# Viscous correction knocks this down ~5–10 % → accept σ ∈ [8.0, 9.5].
#
# Goncharov terminal bubble velocity
#     U_bubble = √(A_t g / ((1 + A_t) k)) ≈ 0.313 m/s → accept 0.25–0.35.
#
# Acceptance
# ----------
# * Linear growth rate within 15 % of 9.07 over t ∈ [0.05, 0.25] s
# * Clear mushroom by t = 1.0 s
# * Mean α_c conserved to machine precision (MULES)
# * α ∈ [0, 1] (MULES-bounded)
#
# Symmetry test (diagnostic, not automated):
#   Re-run with the perturbation centered (A·cos(2π x/W) gives a bubble at x=0
#   by periodicity; to test centred-mushroom, replace with A·cos(π x/W) so the
#   single mushroom sits at x = W/2). It must stay centred. If it walks to a
#   wall, there is a discretisation asymmetry — not a turbulence issue.
# ============================================================
ENV["XCALIBRE_MMP_DIAG"] = "0"
ENV["XCALIBRE_MMP_DIAG_EVERY"] = "1000"

using XCALibre
using Test
using LinearAlgebra
using StaticArrays
using Printf
using Statistics: mean

# ------------------------------------------------------------
# Analytical references
# ------------------------------------------------------------
function rt_growth_rate(A_t, g, k; nu=0.0)
    σ_inv  = sqrt(A_t * g * k)
    # Chandrasekhar viscous RT correction (approximate): σ² + 2νk²σ − A_t g k = 0
    # Solve: σ = −νk² + √(ν²k⁴ + A_t g k)
    σ_visc = -nu*k^2 + sqrt(nu^2 * k^4 + A_t*g*k)
    (inviscid=σ_inv, viscous=σ_visc)
end
goncharov_bubble(A_t, g, k) = sqrt(A_t * g / ((1 + A_t) * k))

# ------------------------------------------------------------
# Simulation setup
# ------------------------------------------------------------
grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
mesh_file = joinpath(grids_dir, "unit_test_stage4.unv")   # swap to RTI_medium.unv for finer
mesh = UNV2D_mesh(mesh_file, scale=1.0)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

# --- Physical parameters ---
ρ_cont = 1000.0   # light (tracked as α_c)
ρ_disp = 2000.0   # heavy
μ_cont = 1.0e-3
μ_disp = 1.0e-3
gvec   = [0.0, -9.81, 0.0]

# --- Domain / perturbation geometry ---
W           = 0.25
H           = 1.0
A_pert      = 0.01 * W          # perturbation amplitude (2.5 mm)
y_mid       = H / 2.0
k_wave      = 2π / W
A_t         = (ρ_disp - ρ_cont) / (ρ_disp + ρ_cont)
σ_theory    = rt_growth_rate(A_t, 9.81, k_wave; nu=μ_cont/ρ_cont)
U_bubble_t  = goncharov_bubble(A_t, 9.81, k_wave)

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(rho=ρ_cont, mu=μ_cont),
            Phase(rho=ρ_disp, mu=μ_disp),
        ),
        gravity = Gravity(gvec),
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev,
)

noSlip = [0.0, 0.0, 0.0]

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:inlet,  noSlip),
            Wall(:outlet, noSlip),
            Wall(:left,   noSlip),
            Wall(:right,  noSlip),
        ],
        p_rgh = [
            Wall(:inlet),
            Wall(:outlet),
            Wall(:left),
            Wall(:right),
        ],
        alpha = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:left),
            Zerogradient(:right),
        ],
    ),
)

schemes = (
    U     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    p     = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
)

solvers = (
    U     = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
    p_rgh = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-10,
                        itmax=2000),
    alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
)

# t_end = 1.5 s at dt = 2e-4 → 7500 iter. Tighten dt if U spikes.
runtime = Runtime(
    iterations=7500, time_step=2.0e-4, write_interval=250)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime,
    hardware=hardware, boundaries=BCs)

GC.gc()

# ------------------------------------------------------------
# Spatial initialisation: tanh interface + cosine perturbation
# α_c (tracked) = 1 − α_d
# ------------------------------------------------------------
cells_host = model.domain.cells
ncells     = length(cells_host)
alpha_host = zeros(Float64, ncells)

# Rough Δy estimate from the domain height / expected cell count (~256 rows).
# Use cached y-extent to set interface thickness ~2 cells.
yc_all = [cells_host[i].centre[2] for i in 1:ncells]
Δy_est = (maximum(yc_all) - minimum(yc_all)) / 256

for i in 1:ncells
    xc = cells_host[i].centre[1]
    yc = cells_host[i].centre[2]
    y_iface = y_mid + A_pert * cos(2π * xc / W)
    α_d = 0.5 * (1.0 + tanh((yc - y_iface) / (2.0 * Δy_est)))
    alpha_host[i] = 1.0 - α_d        # tracked = continuous = light
end

initialise!(model.momentum.U, noSlip)
copyto!(model.fluid.alpha.values, alpha_host)

initial_mean_α = sum(alpha_host .* [cells_host[i].volume for i in 1:ncells]) /
                 sum(cells_host[i].volume for i in 1:ncells)

@info "Running Stage 4 (Rayleigh-Taylor)" A_t k_wave σ_theory U_bubble_t initial_mean_α

@time residuals = run!(model, config, pref=0.0)

# ------------------------------------------------------------
# Post-processing
# ------------------------------------------------------------
α_vals = Array(model.fluid.alpha.values)
U      = model.momentum.U
Uy_vals = Array(U.y.values)

mean_α = mean(α_vals)
min_α  = minimum(α_vals)
max_α  = maximum(α_vals)

# Bin-based interface tracking: at x = 0 (bubble, tracked α_c rising from
# below) and x = W/2 (spike, heavy sinking). We report the y-location of the
# α_c = 0.5 contour in a vertical slice.
x_all = [cells_host[i].centre[1] for i in 1:ncells]
# Late-time-safe tip trackers: by t_end the mushrooms have rolled up so a
# vertical slice has many α=0.5 crossings. Use extremal-y heuristics instead.
#   bubble tip @ x=0    = highest y with α_c ≥ 0.5 (top of rising light plume)
#   spike tip  @ x=W/2  = lowest  y with α_c ≤ 0.5 (bottom of sinking dense finger)
function tip_y_extremal(xc_target, half_width, light_side::Bool)
    # light_side=true  → pick MAX y among cells with α_c ≥ 0.5 (top of light plume → bubble)
    # light_side=false → pick MIN y among cells with α_c ≤ 0.5 (bottom of dense finger → spike)
    y_best = NaN
    for i in 1:ncells
        abs(x_all[i] - xc_target) < half_width || continue
        α = α_vals[i]
        hit = light_side ? (α >= 0.5) : (α <= 0.5)
        hit || continue
        y_best = isnan(y_best) ? yc_all[i] :
                 (light_side ? max(y_best, yc_all[i]) : min(y_best, yc_all[i]))
    end
    y_best
end

y_bubble = tip_y_extremal(0.0, W/16, true)    # bubble tip (light rises at x=0)
y_spike  = tip_y_extremal(W/2, W/16, false)   # spike tip  (dense sinks at x=W/2)
h_bubble = isnan(y_bubble) ? NaN : (y_bubble - y_mid)
h_spike  = isnan(y_spike)  ? NaN : (y_mid - y_spike)

@testset "Stage 4 — Rayleigh-Taylor instability" begin
    @test isapprox(mean_α, initial_mean_α; atol=5e-5)
    @test min_α > -1e-4
    @test max_α <  1 + 1e-4
    # Qualitative: by t_end mushrooms are well developed → amplitudes grew
    # far beyond A_pert = 2.5 mm.
    @test !isnan(h_bubble) && h_bubble > 5 * A_pert
    @test !isnan(h_spike)  && h_spike  > 5 * A_pert

    @printf("\n=== Stage 4 results ===\n")
    @printf("A_t, k, g          = %.3f, %.3f rad/m, 9.81 m/s²\n", A_t, k_wave)
    @printf("σ_inviscid         = %.3f 1/s   (expected growth rate)\n", σ_theory.inviscid)
    @printf("σ_viscous_corr     = %.3f 1/s\n", σ_theory.viscous)
    @printf("U_bubble (Goncharov) = %.3f m/s\n", U_bubble_t)
    @printf("mean(α_c)          = %.6e  (init %.6e, drift %.2e)\n",
            mean_α, initial_mean_α, mean_α - initial_mean_α)
    @printf("min/max α          = %.4e / %.4e\n", min_α, max_α)
    @printf("bubble tip y       = %.4e m   (h = %.4e)\n", y_bubble, h_bubble)
    @printf("spike tip  y       = %.4e m   (h = %.4e)\n", y_spike,  h_spike)
    @printf("max |U|            = %.4e m/s\n", maximum(sqrt.(Array(U.x.values).^2 .+ Uy_vals.^2)))
    @printf("=======================\n")

    @info "Stage 4 summary" A_t σ_theory U_bubble_t mean_α h_bubble h_spike
end

# ------------------------------------------------------------
# Notes
# ------------------------------------------------------------
# 1) Linear growth rate: to measure it quantitatively, re-run with
#    write_interval = 50 (every 0.01 s), extract the α = 0.5 contour at
#    x = 0 or x = W/2 each snapshot, and fit exp(σt) to |y − H/2| over
#    t ∈ [0.05, 0.25] s.  Not automated here — requires on-disk snapshots.
# 2) Symmetry diagnostic: swap `cos(2π x/W)` → `cos(π (x − W/2)/ (W/2))`
#    (single central mushroom). A centred mushroom that walks to a wall
#    indicates a discretisation asymmetry in the mixture solver.
# 3) Grid convergence: switch mesh_file to RTI_medium.unv and confirm the
#    growth rate and bubble speed converge (not a drift).
# 4) Turbulence check: once laminar passes, swap RANS{Laminar}() →
#    RANS{KOmegaSST}() or similar to see how turbulence modelling affects
#    the mixing layer — this is the Becker-relevant case.
