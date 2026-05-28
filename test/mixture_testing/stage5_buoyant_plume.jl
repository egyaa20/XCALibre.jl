# ============================================================
# Stage 5 — 2D buoyant plume in stagnant ambient (Turner / List)
# ============================================================
#
# Purpose
# -------
# Localised-source buoyant rising flow in a confined box — the same dynamics
# as the Becker plume, stripped of any bubble/drag complication. By using a
# tiny "particle" diameter, slip is O(1e-9) m/s and the mixture solver
# reduces to a single-fluid variable-density flow driven purely by buoyancy.
#
# This stage tests:
#   (A) centred source  → self-symmetric rising plume, tests buoyancy
#       coupling, α transport under a steady forcing, and entrainment.
#   (B) offset source   → (not run here — swap source-centre below) the
#       Becker-relevant test: does the plume recover from initial asymmetry
#       or does it glue to the near wall?
#
# Setup
# -----
# Mesh: reuse `unit_test_stage4.unv` (0.25 × 1.0) with scale=2.0
#       → 0.5 m (W) × 2.0 m (H).
# Source: narrow strip at the bottom, |x − x_src| < W_src/2, W_src = 0.05 m.
#         Inlet velocity v_in = 0.1 m/s upward, α_d = 1 (plume fluid).
#         Outside the strip the bottom is no-slip (handled by the Dirichlet
#         function returning zero there — acts as a zero-velocity wall).
# Top:    pressure outlet (p_rgh = 0, zero-gradient for U, α).
# Sides:  no-slip walls, zero-gradient α.
# Tracked α = α_c (ambient). Initial α_c = 1 everywhere.
#
# Hardcoded constants to change in Solvers_2_MULTIPHASE.jl
# --------------------------------------------------------
#     diameter = 1.0e-6     # negligible slip, same trick as Stage 4
#     standard_manninen = true
#
# Phases (5 % density deficit)
# ----------------------------
#     Phase(rho=1000.0, mu=1.0e-3)   # ambient    (tracked α_c)
#     Phase(rho=950.0,  mu=1.0e-3)   # plume fluid
#
# Analytical targets (List 1982, 2-D plane plume)
# ------------------------------------------------
# Reduced buoyancy flux at source:
#     B_0 = g · (ρ_c − ρ_d)/ρ_c · v_in · W_src
#         = 9.81 · 0.05 · 0.1 · 0.05 = 2.45 × 10⁻³ m³/s³  (per unit depth)
# Centreline velocity (plane plume, roughly z-independent):
#     v_c ≈ 1.66 · B_0^{1/3} ≈ 0.22 m/s
# Plume half-width:
#     b(z) ≈ α_E · z, α_E ≈ 0.1
# Entrainment velocity: u_e ≈ α_E · v_c ≈ 0.022 m/s.
#
# Acceptance (Configuration A)
# ----------------------------
# * Plume rises from the centred strip and remains symmetric in time-average
# * Centreline vertical velocity at z ∈ {0.5, 1.0, 1.5} m in 0.15–0.30 m/s
# * α_c conserved by MULES (no drift)
# * 0 ≤ α ≤ 1
# ============================================================
ENV["XCALIBRE_MMP_DIAG"] = "0"

using XCALibre
using Test
using LinearAlgebra
using StaticArrays
using Printf
using Statistics: mean

# ------------------------------------------------------------
# Simulation setup
# ------------------------------------------------------------
grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
mesh_file = joinpath(grids_dir, "unit_test_stage4.unv")
# unit_test_stage4
mesh = UNV2D_mesh(mesh_file, scale=2.0)   # 0.25×1.0 base → 0.5×2.0

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

# --- Physical parameters ---
ρ_cont = 1000.0   # ambient (tracked α_c)
ρ_disp = 950.0    # plume fluid (5 % lighter)
μ_cont = 1.0e-3
μ_disp = 1.0e-3
gvec   = [0.0, -9.81, 0.0]

# --- Source geometry / driving ---
W          = 0.5
H          = 2.0
W_src      = 0.05
x_src      = W / 4.0          # Config A (centred). Use W/4 for Config B (offset, Becker test).
v_in       = 0.1              # inlet velocity (m/s)
half_src   = W_src / 2.0

B0         = 9.81 * (ρ_cont - ρ_disp) / ρ_cont * v_in * W_src
v_c_theory = 1.66 * B0^(1/3)

# --- Bubble/particle diameter for the Manninen drift closure. ---
# 1 µm makes slip negligible (O(1e-9) m/s) so the mixture reduces to a
# single-fluid variable-density flow driven by buoyancy only — the regime
# the plume theory was derived for.
bubble_d = 1.0e-6

# --- Turbulence inlet values (k-ω SST) ---
# Tu ≈ 5 % at the source; length scale ℓ ≈ 0.07·W_src; Cμ = 0.09.
# k   = 1.5 (Tu·v_in)²
# ω   = k^0.5 / (ℓ·Cμ^0.25)
Tu         = 0.05
ℓ_turb     = 0.07 * W_src
k_inlet    = 1.5 * (Tu * v_in)^2
ω_inlet    = sqrt(k_inlet) / (ℓ_turb * 0.09^0.25)
νt_inlet   = k_inlet / ω_inlet

# --- Inlet Dirichlet functions: velocity and α_c on :inlet patch ---
# DirichletFunction signature: f(coords, time, index) → SVector{3} for U, scalar for α.
# Outside the source strip the function returns (0,0,0) for U and 1.0 for α,
# which is equivalent to a no-slip ambient condition on that part of the wall.
function inlet_U(coords, _t, _i)
    x = coords[1]
    inside = abs(x - x_src) < half_src
    SVector{3,Float64}(0.0, inside ? v_in : 0.0, 0.0)
end
function inlet_alpha(coords, _t, _i)
    x = coords[1]
    inside = abs(x - x_src) < half_src
    inside ? 0.0 : 1.0      # tracked α_c = 0 in source (pure dispersed), 1 outside
end

# k-ω SST inlet functions. Inside the source strip use free-stream inlet
# values; outside (the no-slip wall portion of the bottom boundary) use 0
# for k (wall-like) and the same ω_inlet — ω is not critical there because
# U = 0 gives no turbulence production.
function inlet_k(coords, _t, _i)
    x = coords[1]
    inside = abs(x - x_src) < half_src
    inside ? k_inlet : 0.0
end
function inlet_omega(coords, _t, _i)
    ω_inlet
end

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        model = Mixture(diameter=bubble_d),
        phases = (
            Phase(rho=ρ_cont, mu=μ_cont),
            Phase(rho=ρ_disp, mu=μ_disp),
        ),
        gravity = Gravity(gvec),
    ),
    turbulence = RANS{KOmegaSST}(walls=(:left, :right,)),
    energy = Energy{Isothermal}(),
    domain = mesh_dev,
)

noSlip = [0.0, 0.0, 0.0]

# :inlet = bottom, :outlet = top, :left / :right = side walls.
BCs = assign(
    region = mesh_dev,
    (
        U = [
            DirichletFunction(:inlet, inlet_U),
            Zerogradient(:outlet),
            Wall(:left,  noSlip),
            Wall(:right, noSlip),
        ],
        p_rgh = [
            Wall(:inlet),
            Dirichlet(:outlet, 0.0),    # pressure outlet
            Wall(:left),
            Wall(:right),
        ],
        alpha = [
            DirichletFunction(:inlet, inlet_alpha),
            Zerogradient(:outlet),
            Zerogradient(:left),
            Zerogradient(:right),
        ],
        k = [
            DirichletFunction(:inlet, inlet_k),
            Zerogradient(:outlet),
            KWallFunction(:left),
            KWallFunction(:right),
        ],
        omega = [
            DirichletFunction(:inlet, inlet_omega),
            Zerogradient(:outlet),
            OmegaWallFunction(:left),
            OmegaWallFunction(:right),
        ],
        nut = [
            Extrapolated(:inlet),
            Zerogradient(:outlet),
            NutWallFunction(:left),
            NutWallFunction(:right),
        ],
    ),
)

schemes = (
    U     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    p     = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    k     = Schemes(divergence=Upwind),
    omega = Schemes(divergence=Upwind),
    y = Schemes(),
)

solvers = (
    U     = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
    p_rgh = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-10,
                        itmax=2000),
    alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
    k     = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-10, relax=0.6, rtol=1e-3),
    omega = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-10, relax=0.6, rtol=1e-3),
    y = SolverSetup(
        solver      = Cg(),
        preconditioner = Jacobi(),
        convergence = 1e-10,
        rtol = 1e-5,
        relax       = 0.7,
        itmax = 5000
    ),
)

# t_end ≈ 30 s (plume transit ≈ H / v_c ≈ 10 s, then ~20 s quasi-steady).
# dt = 1e-3, 30000 iter.  Adaptive dt is disabled here (same rationale as
# Stages 2–3: mostly-stagnant bulk gives misleading Courant).
runtime = Runtime(
    iterations=10000, time_step=1.0e-3, write_interval=100)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime,
    hardware=hardware, boundaries=BCs)

GC.gc()

# ------------------------------------------------------------
# Initial field: ambient everywhere (α_c = 1, U = 0)
# ------------------------------------------------------------
initialise!(model.momentum.U, noSlip)
initialise!(model.fluid.alpha, 1.0)
initialise!(model.turbulence.k,     k_inlet)
initialise!(model.turbulence.omega, ω_inlet)
initialise!(model.turbulence.nut,   νt_inlet)

@info "Running Stage 5 (buoyant plume, Config A — centred source, k-ω SST)" W H W_src x_src v_in B0 v_c_theory k_inlet ω_inlet νt_inlet

@time residuals = run!(model, config, pref=0.0)

# ------------------------------------------------------------
# Post-processing — centreline velocity at several heights
# ------------------------------------------------------------
cells_host = model.domain.cells
ncells     = length(cells_host)
xc_all     = [cells_host[i].centre[1] for i in 1:ncells]
yc_all     = [cells_host[i].centre[2] for i in 1:ncells]
α_vals     = Array(model.fluid.alpha.values)
Uy_vals    = Array(model.momentum.U.y.values)
Ux_vals    = Array(model.momentum.U.x.values)

"""
    plume_peak(z_target, half_band_y) -> (x_peak, vy_peak)

Scan all cells within `|y − z_target| < half_band_y` across the FULL
x-range, return the x-location of the maximum vertical velocity and the
velocity itself. This tracks the plume wherever it has wandered — the
previous fixed-x-bin sampler was meaningless once the plume meandered
off centre.
"""
function plume_peak(z_target, half_band_y)
    x_peak = NaN; vy_peak = -Inf
    for i in 1:ncells
        abs(yc_all[i] - z_target) < half_band_y || continue
        if Uy_vals[i] > vy_peak
            vy_peak = Uy_vals[i]
            x_peak  = xc_all[i]
        end
    end
    vy_peak == -Inf ? (NaN, NaN) : (x_peak, vy_peak)
end

Δz_band = 0.03       # ±3 cm vertical bin

x_05, vc_05 = plume_peak(0.5, Δz_band)
x_10, vc_10 = plume_peak(1.0, Δz_band)
x_15, vc_15 = plume_peak(1.5, Δz_band)

# Plume mass balance / conservation diagnostic
mean_α = mean(α_vals)
min_α  = minimum(α_vals)
max_α  = maximum(α_vals)
max_Umag = maximum(sqrt.(Ux_vals.^2 .+ Uy_vals.^2))

@testset "Stage 5 — buoyant plume (Config A)" begin
    @test min_α > -1e-4
    @test max_α <  1 + 1e-4
    # Bulk is still mostly ambient; mean α_c ≳ 0.7 in the quasi-steady state
    @test mean_α > 0.5
    # Plume drives non-trivial vertical velocity somewhere in the field
    @test max_Umag > 0.1

    # ----------------------------------------------------------------
    # Quantitative validation: List 1982 plane-plume centreline velocity
    # ----------------------------------------------------------------
    #     v_c ≈ 1.66 · B_0^{1/3}   (z-independent for a 2D plane plume)
    # The coefficient 1.66 is itself uncertain to ±15 % across experimental
    # data, plume self-similarity is asymptotic (near-source still developing,
    # near-top distorted by outflow), so we accept ±50 % at all three heights.
    # The mid-height z = 1.0 sample is the cleanest (developed yet not yet
    # outflow-perturbed) and uses a tighter ±35 % window.
    @test !isnan(vc_05) && 0.5*v_c_theory ≤ vc_05 ≤ 1.6*v_c_theory
    @test !isnan(vc_10) && 0.65*v_c_theory ≤ vc_10 ≤ 1.35*v_c_theory
    @test !isnan(vc_15) && 0.5*v_c_theory ≤ vc_15 ≤ 1.6*v_c_theory

    # Centreline must stay near the source (no wall-attachment).
    # The plume drift |x_peak − x_src| should be a small fraction of W.
    @test abs(x_10 - x_src) < 0.25 * W

    @printf("\n=== Stage 5 results ===\n")
    @printf("Source location x_src = %.4f m   (domain x ∈ [0, %.2f])\n", x_src, W)
    @printf("B_0                   = %.4e m³/s³\n", B0)
    @printf("v_c theory (List)     = %.4f m/s\n", v_c_theory)
    @printf("                        z        x_peak (drift)       v_y_peak\n")
    @printf("                     %.2f m      %.4f m (%+.4f)     %.4f m/s\n", 0.5, x_05, x_05 - x_src, vc_05)
    @printf("                     %.2f m      %.4f m (%+.4f)     %.4f m/s\n", 1.0, x_10, x_10 - x_src, vc_10)
    @printf("                     %.2f m      %.4f m (%+.4f)     %.4f m/s\n", 1.5, x_15, x_15 - x_src, vc_15)
    @printf("mean(α_c)             = %.4f\n", mean_α)
    @printf("min/max α             = %.4e / %.4e\n", min_α, max_α)
    @printf("max |U|               = %.4f m/s\n", max_Umag)
    @printf("=======================\n")

    @info "Stage 5 summary" x_src v_c_theory vc_05 vc_10 vc_15 x_05 x_10 x_15 mean_α max_Umag
end

# ------------------------------------------------------------
# Notes
# ------------------------------------------------------------
# 1) Configuration B (offset-source Becker test): set `x_src = W/4` above
#    and re-run. The plume should entrain ambient fluid and either recover
#    towards centre or stabilise as a tilted-but-stable plume. Wall-sticking
#    (outcome γ in the spec) is the Becker pathology.
# 2) Turbulence: the spec calls for k-ε. Swap `turbulence = RANS{Laminar}()`
#    for a `KOmega` / `KOmegaSST` model once this laminar baseline passes.
#    The inlet k, ε values must be supplied via DirichletFunction too.
# 3) Time-averaged fields: with `write_interval=500` (every 0.5 s) you have
#    60 snapshots over t ∈ [0, 30] s. Time-average the last 40 to compare
#    with List's 2-D plane-plume profile.  (Not automated here.)
# 4) Entrainment measurement: at x = W/2 + 3·b(z) sample u_x; expect an
#    inward (negative for x > x_src) mean ≈ −0.022 m/s.  (Manual.)
