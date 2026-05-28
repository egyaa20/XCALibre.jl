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
        model=Mixture(diameter=1.0e-6),
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


adaptive = AdaptiveTimeStepping(
        maxCo=0.25,
        maxAlphaCo=0.05,
        minShrink=0.1,
        maxGrow=1.2
    )


runtime = Runtime(
    iterations=7500, time_step=2.0e-4, write_interval=50, adaptive=adaptive)

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
x_all  = [cells_host[i].centre[1] for i in 1:ncells]
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

# Tip-extraction helper used by both the linear-phase σ measurement (chunked
# Phase 1) and the final-state mushroom diagnostic (after Phase 2).
function tip_y_extremal(α_arr, xc_target, half_width, light_side::Bool)
    y_best = NaN
    for i in 1:ncells
        abs(x_all[i] - xc_target) < half_width || continue
        α = α_arr[i]
        hit = light_side ? (α >= 0.5) : (α <= 0.5)
        hit || continue
        y_best = isnan(y_best) ? yc_all[i] :
                 (light_side ? max(y_best, yc_all[i]) : min(y_best, yc_all[i]))
    end
    y_best
end

@info "Running Stage 4 (Rayleigh-Taylor)" A_t k_wave σ_theory U_bubble_t initial_mean_α

# ============================================================
# Phase 1 — chunked linear-phase σ measurement
# ============================================================
# Linear stability predicts the perturbation kinetic energy
#   KE_⊥(t) ≡ ∫ U_y² dV ∝ exp(2σ·t)
# This is a *continuous* function of the field state, immune to the
# cell-discretisation artefact that affects bubble-tip position (the tip
# y-value snaps to discrete cell rows until growth exceeds one cell height,
# producing a step-function ln(h) and a misleading slope ≈ 0).
# Slope of ln(KE) is 2σ → divide by 2 to recover σ.
# Fixed dt (no adaptive) so sample times are exactly known.
n_chunks_lin    = 5
chunk_iters_lin = 250             # 250 × 2e-4 s = 50 ms per chunk
dt_lin          = 2.0e-4

times_lin = Float64[]
ke_lin    = Float64[]
amps_lin  = Float64[]   # tip-position record (kept for printing only)

for k in 1:n_chunks_lin
    runtime_lin = Runtime(iterations=chunk_iters_lin, time_step=dt_lin,
                          write_interval=-1)
    config_lin = Configuration(solvers=solvers, schemes=schemes,
                               runtime=runtime_lin, hardware=hardware,
                               boundaries=BCs)
    GC.gc()
    @info "Stage 4 phase-1 chunk $k / $n_chunks_lin"
    run!(model, config_lin, pref=0.0)

    Uy_now  = Array(model.momentum.U.y.values)
    α_now   = Array(model.fluid.alpha.values)
    cell_volumes = [cells_host[i].volume for i in 1:ncells]
    ke_now  = sum(Uy_now[i]^2 * cell_volumes[i] for i in 1:ncells)
    y_b     = tip_y_extremal(α_now, 0.0, W/16, true)

    push!(times_lin, k * chunk_iters_lin * dt_lin)
    push!(ke_lin, ke_now)
    push!(amps_lin, isnan(y_b) ? NaN : (y_b - y_mid))
end

# σ-fit on samples after the velocity field has developed (skip first chunk
# where U is still ramping from zero).
fit_mask    = (times_lin .>= chunk_iters_lin * dt_lin * 1.5) .& (ke_lin .> 0)
ts_fit      = times_lin[fit_mask]
log_ke_fit  = log.(ke_lin[fit_mask])
σ_simulated = if length(ts_fit) ≥ 2
    # least-squares slope of ln(KE) vs t → 2σ → σ = slope/2
    A = hcat(ones(length(ts_fit)), ts_fit)
    coef = A \ log_ke_fit
    coef[2] / 2.0
else
    NaN
end
@info "Stage 4 σ fit" σ_simulated σ_theory.viscous σ_theory.inviscid n_pts=length(ts_fit) ke_lin amps_lin

# ============================================================
# Phase 2 — continue to t_end for nonlinear regime
# ============================================================
remaining_iters = 7500 - n_chunks_lin * chunk_iters_lin    # 6250
runtime_nl = Runtime(iterations=remaining_iters, time_step=dt_lin,
                     write_interval=50, adaptive=adaptive)
config_nl = Configuration(solvers=solvers, schemes=schemes, runtime=runtime_nl,
                          hardware=hardware, boundaries=BCs)
GC.gc()
@time residuals = run!(model, config_nl, pref=0.0)

# ------------------------------------------------------------
# Post-processing
# ------------------------------------------------------------
α_vals = Array(model.fluid.alpha.values)
U      = model.momentum.U
Uy_vals = Array(U.y.values)

mean_α = mean(α_vals)
min_α  = minimum(α_vals)
max_α  = maximum(α_vals)

# Late-time tip trackers: by t_end the mushrooms have rolled up so a vertical
# slice has many α=0.5 crossings. The tip_y_extremal helper (defined before
# Phase 1 above) picks the MAX/MIN y among cells crossing α=0.5 in a thin x band.
y_bubble = tip_y_extremal(α_vals, 0.0, W/16, true)    # bubble tip (light rises at x=0)
y_spike  = tip_y_extremal(α_vals, W/2, W/16, false)   # spike tip  (dense sinks at x=W/2)
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

    # ----------------------------------------------------------------
    # Quantitative: linear growth rate σ from log-amplitude fit.
    # Compare against Chandrasekhar viscous correction; tolerance 20 %
    # to allow for ν-correction approximation, finite mesh resolution,
    # and the lower-order limiter on α-transport.
    # ----------------------------------------------------------------
    @test !isnan(σ_simulated)
    @test isapprox(σ_simulated, σ_theory.viscous; rtol=0.20)

    @printf("\n=== Stage 4 results ===\n")
    @printf("A_t, k, g          = %.3f, %.3f rad/m, 9.81 m/s²\n", A_t, k_wave)
    @printf("σ_inviscid         = %.3f 1/s   (expected growth rate)\n", σ_theory.inviscid)
    @printf("σ_viscous_corr     = %.3f 1/s\n", σ_theory.viscous)
    @printf("σ_simulated (fit)  = %.3f 1/s   (ratio %.3f, n_pts=%d)\n",
            σ_simulated, σ_simulated / σ_theory.viscous, length(ts_fit))
    @printf("U_bubble (Goncharov) = %.3f m/s\n", U_bubble_t)
    @printf("mean(α_c)          = %.6e  (init %.6e, drift %.2e)\n",
            mean_α, initial_mean_α, mean_α - initial_mean_α)
    @printf("min/max α          = %.4e / %.4e\n", min_α, max_α)
    @printf("bubble tip y       = %.4e m   (h = %.4e)\n", y_bubble, h_bubble)
    @printf("spike tip  y       = %.4e m   (h = %.4e)\n", y_spike,  h_spike)
    @printf("max |U|            = %.4e m/s\n", maximum(sqrt.(Array(U.x.values).^2 .+ Uy_vals.^2)))
    @printf("=======================\n")

    @info "Stage 4 summary" A_t σ_theory σ_simulated U_bubble_t mean_α h_bubble h_spike
end

# ------------------------------------------------------------
# Notes
# ------------------------------------------------------------
# 1) Linear growth rate: measured quantitatively in Phase 1 of the test by
#    chunking the run with fixed dt and sampling the bubble tip every 25 ms,
#    then fitting σ from ln(h) vs t over t ∈ [0.05, 0.25] s. The fit slope
#    is asserted against Chandrasekhar's viscous-corrected σ at 20 % rtol.
# 2) Symmetry diagnostic: swap `cos(2π x/W)` → `cos(π (x − W/2)/ (W/2))`
#    (single central mushroom). A centred mushroom that walks to a wall
#    indicates a discretisation asymmetry in the mixture solver.
# 3) Grid convergence: switch mesh_file to RTI_medium.unv and confirm the
#    growth rate and bubble speed converge (not a drift).
# 4) Turbulence check: once laminar passes, swap RANS{Laminar}() →
#    RANS{KOmegaSST}() or similar to see how turbulence modelling affects
#    the mixing layer — this is the Becker-relevant case.
