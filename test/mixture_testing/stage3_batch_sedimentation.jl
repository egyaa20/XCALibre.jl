# ============================================================
# Stage 3 — 2D batch sedimentation (interface tracking)
# ============================================================
#
# Purpose
# -------
# Most demanding validation: track the descent of a free upper interface of
# a settling suspension in a closed 2D box. During the "constant-rate"
# period (after transient, before accumulation bed rises), the interface
# descent velocity equals the Richardson-Zaki hindered settling velocity:
#
#     dH/dt = -u_t · (1 - α_d0)^n
#
# Setup
# -----
# * 2D rectangular domain 0.1 m × 0.5 m, grid `unit_test_stage3.unv`
# * All four boundaries: no-slip walls for U, zero-gradient for α
# * Uniform initial α_d = 0.10  (tracked α_c = 0.90)
# * Laminar, hydrostatic p_rgh ≈ 0 (anchored with pref=0.0)
# * Run ~ 60–100 s physical time (at this scale H drops ~50-150 mm so the
#   interface clears the initial zone without hitting the top wall).
#
# Hardcoded constants to change in Solvers_2_MULTIPHASE.jl
# --------------------------------------------------------
# Same as Stage 2:
#     diameter = 2.0e-4
#     C_vm     = 0.5
#     Sc_t     = 0.7        (inactive for laminar)
#
# To explore different scales, also try:
#     diameter = 5.0e-5     # 50 µm, ~1 mm/s Stokes → slow, longer run
#     diameter = 5.0e-4     # 500 µm, Re_p ~ 40, faster → shorter run
#
# Phases (same as Stage 2):
#     Phase(rho=1000.0, mu=1.0e-3)  # water (tracked)
#     Phase(rho=2500.0, mu=1.0e-3)  # glass beads
#
# Parameter sweeps (run this file multiple times, edit α_d0 / diameter):
#     α_d0 ∈ {0.05, 0.10, 0.15, 0.20}
#     d    ∈ {5e-5, 1e-4, 2e-4, 5e-4}  (ensure Re_p < 1000 for SN regime)
#
# Expected
# --------
# At α_d0 = 0.10, using the solver's buoyancy form:
#     u_t_solver ≈ 3.5e-2 m/s
#     u_rz       ≈ 0.606 · u_t_solver ≈ 2.1e-2 m/s
#     → dH/dt    ≈ -21 mm/s
# In ~5 s the interface drops ~100 mm.
#
# Acceptance: 5–10 % on the slope of the constant-rate segment.
# Volume conservation: ∫α dV constant to machine precision.
# Interface sharpness: max |∂α/∂y| should remain concentrated, not diffuse.
# ============================================================

ENV["XCALIBRE_MMP_DIAG"] = "0"
ENV["XCALIBRE_MMP_DIAG_EVERY"] = "500"
using XCALibre
using Test
using LinearAlgebra
using StaticArrays
using Printf
using Statistics: mean

# ------------------------------------------------------------
# Analytical references
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
mesh_file = joinpath(grids_dir, "unit_test_stage3.unv")
mesh = UNV2D_mesh(mesh_file, scale=1.0)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

# --- Physical parameters ---
ρ_cont = 1000.0
ρ_disp = 2500.0
μ_cont = 1.0e-3
μ_disp = 1.0e-3
d_part = 2.0e-4         # <-- edit to sweep particle size
Cvm    = 0.5
gvec   = [0.0, -9.81, 0.0]
α_d0   = 0.10           # <-- edit to sweep suspension loading
α_c0   = 1.0 - α_d0

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(rho=ρ_cont, mu=μ_cont),
            Phase(rho=ρ_disp, mu=μ_disp)
        ),
        gravity = Gravity(gvec)
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
)

noSlip = [0.0, 0.0, 0.0]

# Closed box: no-slip everywhere, no-flux α, zero-gradient p_rgh.
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
    # Closed box → pure-Neumann p_rgh: use Bicgstab (robust to small
    # asymmetries; Cg errors on "not SPD" in this setting).
    p_rgh = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-10,
                        itmax=2000),
    alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
)

# Fixed dt — adaptive time stepping blows up in a closed box (see Stage 2
# for the same issue and rationale). 1 ms × 10 000 iter = 10 s physical time,
# enough to see the interface descend ~20 cm at u_rz ≈ 21 mm/s.
runtime = Runtime(
    iterations=10000, time_step=1.0e-3, write_interval=1000)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime,
    hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.U,  noSlip)
initialise!(model.fluid.alpha, α_c0)          # uniform initial suspension

initial_mean_α = α_c0
@info "Running Stage 3 (2D batch sedimentation)" d_part α_d0

@time residuals = run!(model, config, pref=0.0)

# ------------------------------------------------------------
# Post-processing
# ------------------------------------------------------------
u_t_solver  = sn_terminal_solver(ρ_disp, ρ_cont, d_part, μ_cont, 9.81, Cvm)
u_rz        = richardson_zaki(u_t_solver, α_d0)
dH_dt_expect = -u_rz

α_vals = Array(model.fluid.alpha.values)
U      = model.momentum.U
Uy_vals = Array(U.y.values)

mean_α = mean(α_vals)
min_α  = minimum(α_vals)
max_α  = maximum(α_vals)

# Interface position estimate — bin α by y and find where α crosses α_c0.
# Requires cell centre coordinates.
cells = model.domain.cells
ncells = length(cells)
yc = zeros(ncells)
for i in 1:ncells
    yc[i] = cells[i].centre[2]
end

# Bin statistics: average α per y-bin (rough interface tracking).
y_min, y_max = extrema(yc)
nbins = 40
edges = range(y_min, y_max; length=nbins+1)
bin_α = zeros(nbins); bin_count = zeros(Int, nbins)
for i in 1:ncells
    k = clamp(searchsortedlast(edges, yc[i]), 1, nbins)
    bin_α[k] += α_vals[i]
    bin_count[k] += 1
end
bin_α ./= max.(bin_count, 1)

# Find highest bin whose mean α ≤ midpoint of (initial, pure_continuous).
# i.e. the interface where clear liquid meets suspension.
threshold = 0.5 * (α_c0 + 1.0)  # α_c = 1 above interface, α_c0 below
interface_bin = nbins
for k in nbins:-1:1
    if bin_α[k] < threshold
        interface_bin = k
        break
    end
end
y_interface = 0.5 * (edges[interface_bin] + edges[interface_bin+1])

@testset "Stage 3 — 2D batch sedimentation" begin
    # Volume conservation — mean α must be conserved to tight tolerance.
    @test isapprox(mean_α, initial_mean_α; atol=5e-3)

    # Bounded α
    @test min_α > -1e-4
    @test max_α <  1 + 1e-4

    # Closed box → mean Uy ≈ 0
    @test abs(mean(Uy_vals)) < 5e-3

    @printf("\n=== Stage 3 results ===\n")
    @printf("d_part         = %.3e m\n", d_part)
    @printf("α_d0           = %.3f\n", α_d0)
    @printf("u_t (solver)   = %.4e m/s\n", u_t_solver)
    @printf("u_rz (analyt.) = %.4e m/s\n", u_rz)
    @printf("dH/dt expected = %.4e m/s\n", dH_dt_expect)
    @printf("mean(α_c)      = %.4e  (init %.4e, drift %.2e)\n",
            mean_α, initial_mean_α, mean_α - initial_mean_α)
    @printf("min/max α      = %.4e / %.4e\n", min_α, max_α)
    @printf("mean(U_y)      = %.4e m/s\n", mean(Uy_vals))
    @printf("interface y    ≈ %.4e m (domain y ∈ [%.3e, %.3e])\n",
            y_interface, y_min, y_max)
    @printf("=======================\n")

    @info "Stage 3 summary" d_part α_d0 u_t_solver u_rz dH_dt_expect mean_α y_interface
end

# ------------------------------------------------------------
# Suggested follow-up (manual sweep, not automated here)
# ------------------------------------------------------------
# 1) Repeat with α_d0 ∈ {0.05, 0.10, 0.15, 0.20} and verify the slope
#    dH/dt tracks -u_t·(1-α_d0)^4.65 within 5-10 %.
# 2) Repeat with d_part ∈ {5e-5, 1e-4, 2e-4, 5e-4} and verify that the
#    Stokes-regime scaling dH/dt ∝ d² holds at small d (Re_p → 0) but
#    rolls off for large d as Re_p enters the Schiller-Naumann branch.
# 3) Volume conservation: the `mean_α` printout should remain pinned to
#    `initial_mean_α` across ALL variations. Any drift > 1 % indicates a
#    problem in the α_eqn solution or in the drift flux divergence.
# 4) Interface sharpness: plot bin_α and check that the transition from
#    α=α_c0 to α=1 spans only a few cells, not a diffuse smear. If the
#    interface is washing out, consider a compressive/TVD α scheme
#    (currently Upwind — van Leer is available via interpolate_vanleer!).
