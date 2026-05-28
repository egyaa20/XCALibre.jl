# ============================================================
# Stage 2 — Hindered (Richardson-Zaki) settling validation
# ============================================================
#
# Purpose
# -------
# Verify that the mixture solver reproduces the Richardson-Zaki hindered
# settling law at uniform volume fraction:
#
#     u_s(α_d) = u_t · (1 - α_d)^n ,   n ≈ 3.78 (Re_p ≈ 4 — Garside-Al-Dibouni)
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
# Expected result (standard Manninen buoyancy: (ρ_d - ρ_m)/ρ_d, α_d→0)
# --------------------------------------------------------------------
#     buoy   = (ρ_d - ρ_c)/ρ_d = 1500/2500 = 0.6
#     u_t    ≈ (ρ_d d²/18μ_c) * 0.6 * g / f_drag(Re≈4.6 → f≈1.43) ≈ 2.29e-2 m/s
# At α_d = 0.10, (1 - α_d)^3.78 ≈ 0.661
#     u_rz   ≈ 0.606 * u_t ≈ 1.39e-2 m/s (~14 mm/s)
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
richardson_zaki(u_t, α_d; n=3.78) = u_t * (1 - α_d)^n

function sn_factor(Re)
    f = Re < 1000 ? 1 + 0.15*Re^0.687 : 0.0183*Re
    max(f, 1.0)
end

"""
    sn_terminal_manninen(ρd, ρc, d, μc, g; tol=1e-12, itmax=500)

Schiller-Naumann fixed-point terminal velocity for the standard Manninen
buoyancy form `(ρ_d - ρ_m)/ρ_d`, evaluated at α_d → 0 so ρ_m = ρ_c.
Matches the solver kernel when invoked with `standard_manninen=true`.
"""
function sn_terminal_manninen(ρd, ρc, d, μc, g; tol=1e-12, itmax=500)
    ρm = ρc
    τd = ρd * d^2 / (18*μc)
    buoy = (ρd - ρm) / ρd
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
mesh_file = joinpath(grids_dir, "batch_medium.unv")
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
        model = Mixture(diameter = d_part),
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

u_t_solver = sn_terminal_manninen(ρ_disp, ρ_cont, d_part, μ_cont, 9.81)
u_rz       = richardson_zaki(u_t_solver, α_d0)

# Sample mixture velocity U (closed box → should stay near 0)
U = model.momentum.U
Uy = Array(U.y.values)
mean_U = sum(Uy) / length(Uy)

# Sample α field (should stratify but mean is conserved)
α_vals = Array(model.fluid.alpha.values)
mean_α = sum(α_vals) / length(α_vals)
ncells = length(α_vals)

# ------------------------------------------------------------
# Kernel diagnostic — probe compute_Ur! independently of the simulation
# ------------------------------------------------------------
# The full simulation result is sensitive to α-eqn discretisation, MULES,
# wall-zeroing, and aspect ratio. To isolate the kernel's behaviour we call
# compute_Ur! directly on a small mesh seeded with α_c = α_c0 and ρ_m derived
# from blend_properties. The fixed-point Ur should equal u_t_kernel*α_c^(n-1)
# with hindering, or u_t_kernel without — where u_t_kernel uses the actual
# Manninen buoyancy (ρ_d−ρ_m)/ρ_d at the seeded ρ_m, not the dilute limit.
let
    α_probe   = ScalarField(mesh_dev); initialise!(α_probe, α_c0)
    ρ_probe   = ScalarField(mesh_dev)
    XCALibre.Solvers.blend_properties!(ρ_probe, α_probe, ρ_cont, ρ_disp)
    Ur_probe  = VectorField(mesh_dev); initialise!(Ur_probe, [0.0, 0.0, 0.0])
    DUm_probe = VectorField(mesh_dev); initialise!(DUm_probe, [0.0, 0.0, 0.0])
    g_probe   = SVector{3,Float64}(0.0, -9.81, 0.0)
    τd_probe  = ρ_disp * d_part^2 / (18 * μ_cont)
    cfg_probe = config

    function picard_u_slip(buoy_eff, τd, ρc_, d_, μc_, g_; itmax=500, tol=1e-12)
        u = 0.02
        for _ in 1:itmax
            Re = ρc_ * abs(u) * d_ / μc_
            f  = max(1.0, 1 + 0.15 * Re^0.687)
            u_new = (τd / f) * buoy_eff * g_
            abs(u_new - u) < tol && return u_new
            u = u_new
        end
        return u
    end

    function probe_kernel(hindering_on::Bool)
        initialise!(Ur_probe, [0.0, 0.0, 0.0])
        for _ in 1:200
            XCALibre.Solvers.compute_Ur!(
                Ur_probe, α_probe, ρ_probe, g_probe, DUm_probe,
                ρ_cont, ρ_disp, μ_cont, d_part, τd_probe,
                true, hindering_on, 3.78, cfg_probe,
            )
        end
        uy_probe = Array(Ur_probe.y.values)
        idx_nz   = findfirst(x -> abs(x) > 1e-12, uy_probe)
        u_slip_kernel = idx_nz === nothing ? 0.0 : abs(uy_probe[idx_nz])

        ρ_sample  = Array(ρ_probe.values)[1]
        buoy_eff  = (ρ_disp - ρ_sample) / ρ_disp
        if hindering_on
            buoy_eff *= α_c0^2.78   # (n_RZ - 1)
        end
        # Kernel returns Ur = U_dm / α_c (slip velocity, not diffusion velocity),
        # so divide the analytical U_dm by α_c to compare like-for-like.
        u_slip_ref = picard_u_slip(buoy_eff, τd_probe, ρ_cont, d_part, μ_cont, 9.81) / α_c0

        @info "[probe] compute_Ur! at α_c=$α_c0" hindering=hindering_on ρ_m_seeded=ρ_sample u_slip_kernel u_slip_analytical=u_slip_ref ratio=u_slip_kernel/u_slip_ref
    end

    probe_kernel(true)
    probe_kernel(false)

    # Also sample the simulation's actual ρ field at end-of-run to verify
    # blend_properties has been keeping it consistent with α.
    ρ_sim = Array(model.fluid.rho.values)
    @info "[probe] simulation rho field" mean=sum(ρ_sim)/length(ρ_sim) min=minimum(ρ_sim) max=maximum(ρ_sim)
end

# ------------------------------------------------------------
# Quantitative R-Z check via mass-based interface descent
# ------------------------------------------------------------
# In a closed box with mixture U≈0 the dispersed-phase volumetric flux
# through any horizontal plane y* is:
#     J_d = α_d (1 - α_d) u_slip = α_d · V_interface
# where V_interface = (1-α_d)·u_slip is the Kynch shock speed of the upper
# (clear/suspension) interface. With the kernel hindered as u_slip = u_t·α_c^(n-1),
# V_interface = u_t·α_c^n = u_rz exactly.
#
# Integrating J_d through the mid-plane y_mid from t=0 to T_phys equals the
# mass that has left the upper half. Since the upper half started at α_d0:
#     f_clear = 1 - <α_d>_upper / α_d0
# and conservation gives  Δh_descent = f_clear · (y_max - y_mid),
# yielding V_int = Δh_descent / T_phys, *independent of upwind smearing of the
# interface profile*. Valid while the interface stays above y_mid (which
# requires u_rz · T_phys < (y_max - y_mid)/2 — true for both Stage 2 & 3).
cells       = model.domain.cells
yc          = [cells[i].centre[2] for i in 1:ncells]
cell_vols   = [cells[i].volume    for i in 1:ncells]
y_min, y_max = extrema(yc)
H_col       = y_max - y_min
y_mid       = 0.5 * (y_min + y_max)

upper_mask     = yc .> y_mid
V_upper        = sum(cell_vols[upper_mask])
α_d_vals       = 1.0 .- α_vals
mass_d_upper   = sum(α_d_vals[upper_mask] .* cell_vols[upper_mask])
α_d0_val       = 1.0 - α_c0
f_clear        = 1.0 - (mass_d_upper / V_upper) / α_d0_val
T_phys         = runtime.iterations * runtime.dt[1]
descent_dist   = f_clear * (y_max - y_mid)
descent_rate   = descent_dist / T_phys
descent_expect = u_rz

@testset "Stage 2 — hindered settling" begin
    # Volume conservation on α (tracked = continuous → should stay near α_c0)
    @test isapprox(mean_α, α_c0; atol=2e-3)

    # Closed box: mean mixture velocity should be tiny
    @test abs(mean_U) < 1e-3

    # Alpha bounded
    @test all(α_vals .>= -1e-6)
    @test all(α_vals .<= 1 + 1e-6)

    # Quantitative R-Z: mass-based interface descent matches u_rz to ~10%.
    # Tolerance covers numerical diffusion at the mid-plane (small) + small
    # start-up transient. With the (n-1) hindering form in the kernel this
    # should agree with R-Z directly.
    @test descent_dist > 0
    @test isapprox(descent_rate, descent_expect; rtol=0.10)

    @printf("u_t (solver, Stage 2):  %.4e m/s\n", u_t_solver)
    @printf("u_s Richardson-Zaki:    %.4e m/s\n", u_rz)
    @printf("mean(U_y) (should ≈0):  %.4e m/s\n", mean_U)
    @printf("mean(α_c):              %.4e  (init %.4e)\n", mean_α, α_c0)
    @printf("f_clear (upper half):   %.4e  (= mass-based clear-water fraction)\n", f_clear)
    @printf("descent (mass-based):   %.4e m/s  (expected %.4e, ratio %.3f, t=%.2fs)\n",
            descent_rate, descent_expect, descent_rate / descent_expect, T_phys)

    @info "Stage 2 summary" u_t_solver u_rz mean_U mean_α α_c0 f_clear descent_rate descent_expect
end
