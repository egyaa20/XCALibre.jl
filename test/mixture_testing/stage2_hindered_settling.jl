# ============================================================
# Stage 2 — Hindered (Richardson-Zaki) settling validation
# ============================================================
#
# Purpose
# -------
# Verify that the mixture solver reproduces the Richardson-Zaki hindered
# settling law at uniform volume fraction:
#
#     u_s(α_d) = u_t · (1 - α_d)^n ,   n ≈ 4.65 (low Re, Re_p < 0.2)
#
# Setup
# -----
# * Tall column ~ 0.05 m (width) × 1.0 m (height), grid `unit_test_stage2.unv`
# * Uniform initial α_d = 0.10 (→ α_continuous = 0.90)
# * Side walls (`:left`, `:right`) : Symmetry  (free-slip, no drift across)
# * Top    (`:outlet` or `:top`)   : Wall      (no-flux boundary)
# * Bottom (`:inlet`  or `:bottom`): Wall      (no-flux boundary)
# * Laminar turbulence (no dispersion)
# * Zero gravity anchor: grid origin near column centre (so |g·h| is bounded)
#
# In the absence of through-flow, the mixture U should stay ≈ 0 (closed box)
# while the slip field Ur drives α settling at a speed given by Richardson-Zaki.
# After a short transient we sample mean |Ur_y| over interior cells and compare
# to the analytical hindered settling velocity.
#
# Hardcoded constants to change in Solvers_2_MULTIPHASE.jl
# --------------------------------------------------------
# The test uses 200 µm glass beads in water (matching Stage 1). Edit
# `Solvers_2_MULTIPHASE.jl` near line ~180 so that:
#
#     diameter = 2.0e-4     # 200 µm  (was 3.0e-3 for 3 mm bubbles)
#     C_vm     = 0.5        # default
#     Sc_t     = 0.7        # irrelevant for laminar but keep consistent
#
# Phases in this script:
#     Phase(rho=1000.0, mu=1.0e-3)   # water  → primary (α tracked = water fraction)
#     Phase(rho=2500.0, mu=1.0e-3)   # glass beads (μ_d only affects mixture μ blend)
#
# Expected result
# ---------------
# u_t (solver's buoyancy form, α_d→0):
#     (ρ_d-ρ_c)/(ρ_c + C_vm·ρ_c) = 1500/1500 = 1.0
#     u_t_solver ≈ (ρ_d d²/18μ) * 1 * g / f_drag(Re=6.5≈1.58) ≈ 3.46e-2 m/s
# At α_d = 0.10, (1 - α_d)^4.65 ≈ 0.606
#     u_s_expected ≈ 0.606 * u_t_solver ≈ 2.10e-2 m/s (~21 mm/s)
#
# Acceptance: 5–10 % on the slip magnitude.
# Interpretation:
#   - Much larger than expected → buoyancy form issue (see Stage 1 note).
#   - Much smaller  than expected → something is damping Ur (e.g. turbulent
#     dispersion not disabled, or alpha clamping is eating gradient at walls).
# ============================================================
ENV["XCALIBRE_MMP_DIAG"] = "1"
ENV["XCALIBRE_MMP_DIAG_EVERY"] = "500"

using XCALibre
using Test
using LinearAlgebra
using StaticArrays
using Printf

# ------------------------------------------------------------
# Reference: Richardson-Zaki hindered-settling speed
# ------------------------------------------------------------
richardson_zaki(u_t, α_d; n=4.65) = u_t * (1 - α_d)^n

function sn_factor(Re)
    f = Re < 1000 ? 1 + 0.15*Re^0.687 : 0.0183*Re
    max(f, 1.0)
end

function sn_terminal_solver(ρd, ρc, d, μc, g, Cvm; tol=1e-12, itmax=500)
    τd = ρd * d^2 / (18*μc)
    buoy = (ρd - ρc) / (ρc + Cvm*ρc)
    u = g*(ρd-ρc)*d^2/(18*μc)
    for _ in 1:itmax
        Re = ρc*abs(u)*d/μc
        f  = sn_factor(Re)
        u_new = (τd / f) * buoy * g
        abs(u_new - u) < tol && return u_new
        u = u_new
    end
    u
end

# ------------------------------------------------------------
# Simulation setup
# ------------------------------------------------------------
grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
mesh_file = joinpath(grids_dir, "unit_test_stage2.unv")
mesh = UNV2D_mesh(mesh_file, scale=1.0)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

# --- Physical parameters (must match diameter set in MULTIPHASE.jl) ---
ρ_cont = 1000.0       # water
ρ_disp = 2500.0       # glass beads
μ_cont = 1.0e-3
μ_disp = 1.0e-3
d_part = 2.0e-4       # 200 µm
Cvm    = 0.5
gvec   = [0.0, -9.81, 0.0]
α_d0   = 0.10         # initial uniform dispersed-phase fraction
α_c0   = 1.0 - α_d0   # tracked (continuous) fraction

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(rho=ρ_cont, mu=μ_cont),   # tracked phase = continuous (water)
            Phase(rho=ρ_disp, mu=μ_disp)    # dispersed     = glass beads
        ),
        gravity = Gravity(gvec)
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
)

noSlip = [0.0, 0.0, 0.0]

# Hermetic column: no through-flow, slip sides, no-flux top/bottom.
BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:inlet,  noSlip),
            Wall(:outlet, noSlip),
            Symmetry(:left,  noSlip),
            Symmetry(:right, noSlip),
        ],
        p_rgh = [
            Wall(:inlet),
            Wall(:outlet),
            Zerogradient(:left),
            Zerogradient(:right),
        ],
        alpha = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:left),
            Zerogradient(:right),
        ],
    )
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
    # Closed box → pure-Neumann p_rgh has a null mode. Use Bicgstab (robust
    # to slight asymmetry) rather than Cg, which errors out when the matrix
    # is not strictly SPD.
    p_rgh = SolverSetup(solver=Bicgstab(),       preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-10,
                        itmax=2000),
    alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
)

# Fixed dt — adaptive time stepping blows up in a closed box because the
# mixture U stays ≈ 0, Courant reads ≈ 0, and dt grows by `maxGrow` per step.
# With u_t ≈ 35 mm/s and cell size ~ few mm, dt = 1 ms gives CFL ≲ 0.02.
runtime = Runtime(
    iterations=5000, time_step=1.0e-3, write_interval=500)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime,
    hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.U,  noSlip)
initialise!(model.fluid.alpha, α_c0)      # tracked field = continuous fraction

@info "Running Stage 2 (hindered settling)"
@time residuals = run!(model, config, pref=0.0)

# ------------------------------------------------------------
# Post-processing — compare sampled |Ur_y| to Richardson-Zaki
# ------------------------------------------------------------

u_t_solver = sn_terminal_solver(ρ_disp, ρ_cont, d_part, μ_cont, 9.81, Cvm)
u_rz       = richardson_zaki(u_t_solver, α_d0)

# Sample mixture velocity U (closed box → should stay near 0)
U = model.momentum.U
Uy = Array(U.y.values)
mean_U = sum(Uy) / length(Uy)

# Sample α field (should stratify but mean is conserved)
α_vals = Array(model.fluid.alpha.values)
mean_α = sum(α_vals) / length(α_vals)

# Look at central region only (avoid wall boundary layers if any)
ncells = length(α_vals)

@testset "Stage 2 — hindered settling diagnostics" begin
    # Volume conservation on α (tracked = continuous → should stay near α_c0)
    @test isapprox(mean_α, α_c0; atol=2e-3)

    # Closed box: mean mixture velocity should be tiny
    @test abs(mean_U) < 1e-3

    # Alpha bounded
    @test all(α_vals .>= -1e-6)
    @test all(α_vals .<= 1 + 1e-6)

    @printf("u_t (solver, Stage 2): %.4e m/s\n", u_t_solver)
    @printf("u_s Richardson-Zaki:   %.4e m/s\n", u_rz)
    @printf("mean(U_y) (should ≈0): %.4e m/s\n", mean_U)
    @printf("mean(α_c):             %.4e  (init %.4e)\n", mean_α, α_c0)

    @info "Stage 2 summary" u_t_solver u_rz mean_U mean_α α_c0
end

# NOTE: A quantitative check against u_rz requires sampling Ur from the
# solver at end-of-run, but Ur is a solver-internal working field (not stored
# on `model`). If you want to enforce a tight slip-velocity tolerance here,
# export Ur via postprocessing or expose it on the residuals NamedTuple.
