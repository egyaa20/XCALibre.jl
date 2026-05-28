# Stage 3 — Mesh-convergence rate of the α-transport scheme.
#
# Sweeps five uniform Cartesian meshes (40²–120²) and runs each at three
# configurations to verify prediction (4) of the modified-equation
# analysis. For a discontinuous (step) initial condition, the L¹ and
# RMS convergence rates are set by how the smear-width σ scales with
# Δx, *not* by the formal order of the scheme. Specifically:
#
#       L¹  ∝ σ                  RMS ∝ √(σ / L_x)
#
# This is the classical Kuznetsov half-order degradation for
# discontinuous solutions of hyperbolic conservation laws (LeVeque,
# Finite Volume Methods, §8.5–8.6).
#
#   Config A — PURE UPWIND        (pure_upwind=true,  cAlpha=0)
#       σ ∝ √Δx       →   L¹ slope ≈ 0.5,    RMS slope ≈ 0.25
#       (modified equation: diffusivity D_num = ½|U|Δx(1−C); spread is
#        an erfc of width √(2 D_num t))
#
#   Config B — VAN-LEER HO        (pure_upwind=false, cAlpha=0)
#       σ ∝ Δx        →   L¹ slope ≈ 1.0,    RMS slope ≈ 0.5
#       (HO flux cancels the leading-order diffusion in the interface
#        band; smear width is a fixed number of cells)
#
#   Config C — COMPRESSED         (pure_upwind=false, cAlpha=1)
#       σ ∝ Δx        →   L¹ slope ≈ 1.0,    RMS slope ≈ 0.5
#       (compression sharpens further; absolute error must be lower
#        than Config B at every mesh)
#
# Two error norms are reported:
#   L¹  = Σ |α_num − α_exact| · Δx       (interface-position error)
#   RMS = √(Σ (α_num − α_exact)² / N)    (profile error)
#
# Cell Courant number is held fixed across meshes (C ≈ 0.04, dt = C·Δx/|U|)
# so the comparison is purely spatial.

using XCALibre
using LinearAlgebra
using Statistics
using Printf

# =============================================================================
# CONFIGURATION
# =============================================================================

# (filename, Δx, label) — domain is 1000 m, Δx = 1000 / N
# Geometric progression (Δx halves each refinement) — equal spacing in
# log Δx gives a well-conditioned slope regression.
MESH_SPECS = [
    ("quad20.unv",  50.0,  "20×20"),
    ("quad40.unv",  25.0,  "40×40"),
    ("quad80.unv",  12.5,  "80×80"),
    ("quad160.unv",  6.25, "160×160"),
]

# Three configurations per mesh: (label, pure_upwind, cAlpha)
CONFIGS = [
    ("upwind",     true,  0.0),
    ("vanleer",    false, 0.0),
    ("compressed", false, 1.0),
]

Lx           = 1000.0
U_FREESTREAM = 1.0
X_INTERFACE  = 500.0
SAMPLING_Y   = 500.0
t_END        = 100.0

COURANT_TARGET = 0.04   # held fixed across meshes

# =============================================================================
# RUN ONE CASE — returns measured errors and α-bound info
# =============================================================================

function run_case(grid::String, Δx::Float64, pure_upwind::Bool, c_alpha::Float64,
                  t_end::Float64)
    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh_file = joinpath(grids_dir, grid)
    mesh = UNV2D_mesh(mesh_file, scale=1.0)

    backend  = CPU(); workgroup = AutoTune(); activate_multithread(backend)
    hardware = Hardware(backend=backend, workgroup=workgroup)
    mesh_dev = adapt(backend, mesh)

    # Hold cell Courant fixed across meshes for fair spatial comparison
    dt = COURANT_TARGET * Δx / U_FREESTREAM
    iterations = max(round(Int, t_end / dt), 1)
    dt_actual  = t_end / iterations

    gravity = Gravity([0.0, 0.0, 0.0])
    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model  = VOF(sigma=0.0, cAlpha=c_alpha, cycles=1, smooth=0,
                         pure_upwind=pure_upwind),
            phases = (
                Phase(rho=1.0, mu=1.0e-6),
                Phase(rho=1.0, mu=1.0e-6),
            ),
            gravity = gravity,
        ),
        turbulence = RANS{Laminar}(),
        energy     = Energy{Isothermal}(),
        domain     = mesh_dev,
    )

    inlet_velocity = [U_FREESTREAM, 0.0, 0.0]
    BCs = assign(
        region = mesh_dev,
        (
            U = [
                Dirichlet(:inlet, inlet_velocity),
                Zerogradient(:outlet),
                Zerogradient(:top),
                Zerogradient(:bottom),
            ],
            p_rgh = [
                Zerogradient(:inlet),
                Dirichlet(:outlet, 0.0),
                Zerogradient(:top),
                Zerogradient(:bottom),
            ],
            alpha = [
                Dirichlet(:inlet, 1.0),
                Zerogradient(:outlet),
                Zerogradient(:top),
                Zerogradient(:bottom),
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
        U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1e-7),
        p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                            convergence=1e-7, relax=1.0, rtol=0.0, atol=1e-9, itmax=5000),
        alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                            convergence=1e-7, relax=1.0, rtol=0.0, atol=1e-7),
    )

    runtime = Runtime(iterations=iterations, time_step=dt_actual,
                      write_interval=iterations)
    config = Configuration(
        solvers=solvers, schemes=schemes, runtime=runtime,
        hardware=hardware, boundaries=BCs)

    GC.gc()

    initialise!(model.fluid.p_rgh, 0.0)
    initialise!(model.momentum.U,  inlet_velocity)
    initialise!(model.fluid.alpha, 0.0)

    setField_Expression!(
        mesh        = mesh,
        field       = model.fluid.alpha,
        condition   = (x, y, z) -> x < X_INTERFACE,
        value_true  = 1.0,
        value_false = 0.0,
    )

    run!(model, config)

    # Extract centre-row α profile (same logic as stages 1 & 2)
    cells_cpu = Array(mesh.cells)
    alpha_cpu = Array(model.fluid.alpha.values)
    cy_all = [cell.centre[2] for cell in cells_cpu]
    y_target = cy_all[argmin(abs.(cy_all .- SAMPLING_Y))]

    sample = Tuple{Float64,Float64}[]
    for (i, cell) in enumerate(cells_cpu)
        cx, cy, _ = cell.centre
        if abs(cy - y_target) < Δx/4
            push!(sample, (cx, alpha_cpu[i]))
        end
    end
    sort!(sample, by = first)

    x_interface_now = X_INTERFACE + U_FREESTREAM * t_end

    # Errors vs analytical sharp step at x = x_interface_now
    L1_err  = 0.0
    SS_err  = 0.0
    for (cx, α) in sample
        α_exact = cx < x_interface_now ? 1.0 : 0.0
        L1_err += abs(α - α_exact) * Δx
        SS_err += (α - α_exact)^2
    end
    RMS_err = sqrt(SS_err / max(length(sample), 1))

    α_min, α_max = extrema(alpha_cpu)

    return (Δx=Δx, dt=dt_actual, iterations=iterations,
            L1_err=L1_err, RMS_err=RMS_err,
            α_min=α_min, α_max=α_max,
            n_samples=length(sample))
end

# =============================================================================
# REGRESSION HELPER
# =============================================================================

function loglog_slope(Δx_vals, err_vals)
    logΔx = log.(Δx_vals)
    logE  = log.(err_vals)
    x̄ = mean(logΔx); ȳ = mean(logE)
    slope = sum((logΔx .- x̄) .* (logE .- ȳ)) / sum((logΔx .- x̄).^2)
    intercept = ȳ - slope * x̄
    return slope, intercept
end

# =============================================================================
# MAIN
# =============================================================================

println("=" ^ 70)
println("Stage 3 — Mesh-convergence rate of α-transport")
println("=" ^ 70)
@printf("Domain Lx = %g m,  U = %g m/s,  t_end = %g s,  C_target = %.3f\n\n",
        Lx, U_FREESTREAM, t_END, COURANT_TARGET)

# results[config_label] :: Vector{NamedTuple}
results = Dict{String, Vector{NamedTuple}}()
for (cfg_label, _, _) in CONFIGS
    results[cfg_label] = NamedTuple[]
end

for (grid, Δx, mesh_label) in MESH_SPECS
    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    if !isfile(joinpath(grids_dir, grid))
        @warn "Skipping $grid — file not found in $grids_dir"
        continue
    end

    for (cfg_label, pure_upwind, c_alpha) in CONFIGS
        @info "Running $mesh_label / $cfg_label  (Δx=$(round(Δx,digits=3)), pure_upwind=$pure_upwind, cAlpha=$c_alpha)"
        @time res = run_case(grid, Δx, pure_upwind, c_alpha, t_END)
        push!(results[cfg_label], (mesh_label=mesh_label, res...))
    end
end

# -----------------------------------------------------------------------------
# Per-configuration tables
# -----------------------------------------------------------------------------
println()
println("=" ^ 70)
println("Error tables")
println("=" ^ 70)
for (cfg_label, _, _) in CONFIGS
    println()
    println("── Config: $cfg_label ──")
    @printf("%-10s  %-8s  %-10s  %-12s  %-12s  %-14s\n",
            "Mesh", "Δx [m]", "iterations", "L¹ error", "RMS error", "α bounds")
    println("-"^78)
    for r in results[cfg_label]
        @printf("%-10s  %-8.3f  %-10d  %-12.4f  %-12.4e  [%-.4f, %.4f]\n",
                r.mesh_label, r.Δx, r.iterations, r.L1_err, r.RMS_err,
                r.α_min, r.α_max)
    end
end

# -----------------------------------------------------------------------------
# Convergence rates
# -----------------------------------------------------------------------------
println()
println("=" ^ 70)
println("Convergence rates  (slope of log(error) vs log(Δx))")
println("=" ^ 70)
@printf("%-12s  %-18s  %-18s\n", "Config", "L¹ slope (intc)", "RMS slope (intc)")
println("-"^60)

slopes = Dict{String, NamedTuple}()
for (cfg_label, _, _) in CONFIGS
    rs = results[cfg_label]
    if length(rs) < 2
        @printf("%-12s  insufficient meshes\n", cfg_label)
        continue
    end
    Δx_vals  = [r.Δx for r in rs]
    L1_vals  = [r.L1_err for r in rs]
    RMS_vals = [r.RMS_err for r in rs]
    pL1 , bL1  = loglog_slope(Δx_vals, L1_vals)
    pRMS, bRMS = loglog_slope(Δx_vals, RMS_vals)
    slopes[cfg_label] = (L1_slope=pL1, L1_intercept=bL1,
                         RMS_slope=pRMS, RMS_intercept=bRMS)
    @printf("%-12s  %.3f  (%.3f)     %.3f  (%.3f)\n",
            cfg_label, pL1, bL1, pRMS, bRMS)
end

# -----------------------------------------------------------------------------
# Pass criteria — predicted slopes from σ vs Δx scaling
#
# For a step IC: L¹ ∝ σ, RMS ∝ √σ. Predicted (L¹, RMS) slopes are:
#   upwind     (σ ∝ √Δx):  (0.5, 0.25)
#   vanleer    (σ ∝ Δx) :  (1.0, 0.5)
#   compressed (σ ∝ Δx) :  (1.0, 0.5)
# -----------------------------------------------------------------------------
println()
println("=" ^ 70)
println("Pass criteria")
println("=" ^ 70)

# (config, predicted_L1_slope, predicted_RMS_slope, tolerance)
PREDICTIONS = [
    ("upwind",     0.5, 0.25, 0.10),
    ("vanleer",    1.0, 0.5 , 0.10),
    ("compressed", 1.0, 0.5 , 0.10),
]

slope_oks = Bool[]
for (cfg, pL1, pRMS, tol) in PREDICTIONS
    s = get(slopes, cfg, nothing)
    if s === nothing
        push!(slope_oks, false)
        @printf("[CRITERION] %-12s slope check (insufficient data):  FAIL\n", cfg)
        continue
    end
    okL1  = abs(s.L1_slope  - pL1)  ≤ tol
    okRMS = abs(s.RMS_slope - pRMS) ≤ tol
    push!(slope_oks, okL1 && okRMS)
    @printf("[CRITERION] %-12s L¹ slope %.3f  vs %.2f±%.2f:  %s\n",
            cfg, s.L1_slope, pL1, tol, okL1 ? "PASS" : "FAIL")
    @printf("[CRITERION] %-12s RMS slope %.3f vs %.2f±%.2f:  %s\n",
            cfg, s.RMS_slope, pRMS, tol, okRMS ? "PASS" : "FAIL")
end

# Compressed should also have strictly lower L¹ than vanleer at every mesh
local errs_lower = true
if length(results["compressed"]) == length(results["vanleer"]) > 0
    for (rC, rB) in zip(results["compressed"], results["vanleer"])
        if rC.L1_err > rB.L1_err
            global errs_lower = false; break
        end
    end
else
    global errs_lower = false
end

# Boundedness across all runs
local bounds_ok = true
for (cfg_label, _, _) in CONFIGS
    for r in results[cfg_label]
        if r.α_min < -1e-3 || r.α_max > 1.0 + 1e-3
            global bounds_ok = false; break
        end
    end
end

println("[CRITERION] compressed L¹ < vanleer L¹ ∀ mesh:   ",
        errs_lower    ? "PASS" : "FAIL")
println("[CRITERION] α ∈ [-1e-3, 1+1e-3] across all runs: ",
        bounds_ok     ? "PASS" : "FAIL")

status = all(slope_oks) && errs_lower && bounds_ok ? "PASS" : "FAIL"
println()
println("Stage 3: $status")

println("""

Modified-equation predictions for a step IC (σ-scaling-driven):
  upwind     →  L¹ slope ≈ 0.5,  RMS slope ≈ 0.25   (σ ∝ √Δx)
  vanleer    →  L¹ slope ≈ 1.0,  RMS slope ≈ 0.5    (σ ∝ Δx)
  compressed →  L¹ slope ≈ 1.0,  RMS slope ≈ 0.5    (σ ∝ Δx, lower abs. error)
""")
