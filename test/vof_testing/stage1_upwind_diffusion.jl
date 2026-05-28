# Stage 1 — Numerical-diffusion test on a planar α step (single mesh).
#
# Verifies predictions (1) and (2) of the modified-equation analysis on
# `quad40.unv` (1000×1000 m, Δx = 25 m). A planar α step at x = 500 m is
# advected by U = (1, 0, 0) m/s for t = 100 s; the interface position at
# t_end is x = 600 m, with the smeared profile fitted to an erfc.
#
# Part A — PURE UPWIND          (pure_upwind=true,  cAlpha=0)
#   The diagnostic toggle in `Solvers_2_MULTIPHASE.jl` collapses
#   alphaf_HO onto alphaf_upwind, so phiAf = phiHf − phiLf = 0 and MULES
#   degenerates to bare 1st-order upwind. The modified equation predicts
#   numerical diffusivity D_num = ½|U|Δx(1−C) and an erfc spread
#       σ_upwind(t) = √(2 D_num t).
#
# Part B — VAN-LEER HO          (pure_upwind=false, cAlpha=0)
#   Default XCALibre VOF: MULES blends LO upwind with HO van-Leer; no
#   compression. The modified equation predicts σ ≪ σ_upwind because
#   the HO flux cancels the leading-order diffusion in the interface
#   band.
#
# Compression at cAlpha = 1 is verified separately in stage 2.
# Mesh convergence across configurations is verified in stage 3.

using XCALibre
using LinearAlgebra
using Statistics
using Printf
using SpecialFunctions

# =============================================================================
# CONFIGURATION
# =============================================================================

GRID         = "quad40.unv"        # 40×40 cells, 1000×1000 m → Δx = 25 m
Lx           = 1000.0
Δx           = 25.0

U_FREESTREAM = 1.0                 # advection velocity, m/s
X_INTERFACE  = 500.0               # initial interface x position, m
SAMPLING_Y   = 500.0               # row to sample for the 1D profile, m

t_END        = 100.0               # physical time, s
dt           = 1.0
ITERATIONS   = round(Int, t_END / dt)

COURANT      = U_FREESTREAM * dt / Δx
D_NUM_PRED   = 0.5 * U_FREESTREAM * Δx * (1.0 - COURANT)
SIGMA_PRED   = sqrt(2.0 * D_NUM_PRED * t_END)

# =============================================================================
# RUN ONE CASE — returns measured σ (erfc-fit width) and α-bound info
# =============================================================================

function run_case(pure_upwind::Bool)
    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh_file = joinpath(grids_dir, GRID)
    mesh = UNV2D_mesh(mesh_file, scale=1.0)

    backend  = CPU(); workgroup = AutoTune(); activate_multithread(backend)
    hardware = Hardware(backend=backend, workgroup=workgroup)
    mesh_dev = adapt(backend, mesh)

    gravity = Gravity([0.0, 0.0, 0.0])
    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model  = VOF(sigma=0.0, cAlpha=0.0, cycles=1, smooth=0,
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

    runtime = Runtime(iterations=ITERATIONS, time_step=dt, write_interval=ITERATIONS)
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

    label = pure_upwind ? "PURE UPWIND" : "MULES + van-Leer HO"
    @info "Stage 1: running $label..."
    @time run!(model, config)

    # Extract centre-row α profile.
    # Find the row whose cell-centre y is closest to SAMPLING_Y, then take
    # all cells within Δx/4 of that y. This is robust to meshes that don't
    # have a cell centre exactly at SAMPLING_Y.
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

    # Diagnostics
    α_global_min, α_global_max = extrema(alpha_cpu)
    α_global_mean = sum(alpha_cpu) / length(alpha_cpu)
    @info @sprintf("Global α stats: min=%.4f  max=%.4f  mean=%.4f", α_global_min, α_global_max, α_global_mean)
    @info "Sampled row (y_target=$y_target, requested $SAMPLING_Y):  $(length(sample)) cells"
    if !isempty(sample)
        @info @sprintf("  x range: [%.1f, %.1f]   α range: [%.4f, %.4f]",
            sample[1][1], sample[end][1],
            minimum(s[2] for s in sample), maximum(s[2] for s in sample))
        # Print first non-1 and first non-0 cells along the row
        for (j, (cx, α)) in enumerate(sample)
            if α < 0.999
                @info @sprintf("  first α<0.999 at j=%d, x=%.1f, α=%.4f", j, cx, α)
                break
            end
        end
        for (j, (cx, α)) in enumerate(sample)
            if α < 0.001
                @info @sprintf("  first α<0.001 at j=%d, x=%.1f, α=%.4f", j, cx, α)
                break
            end
        end
    end

    return (
        xs = [s[1] for s in sample],
        αs = [s[2] for s in sample],
    )
end

# erfc fit
erfc_profile(σ, x, x0) = 0.5 * erfc((x - x0) / (σ * sqrt(2.0)))

function fit_sigma(xs, αs, x0; σmax)
    σ_grid = collect(range(0.05*Δx, σmax, length=2000))
    errs = [sum((erfc_profile(σ, x, x0) - α)^2 for (x, α) in zip(xs, αs))
            for σ in σ_grid]
    return σ_grid[argmin(errs)]
end

# =============================================================================
# MAIN — run both cases and report
# =============================================================================

x_now = X_INTERFACE + U_FREESTREAM * t_END

println("=" ^ 70)
println("Stage 1 — Interface diffusion under MULES α-transport")
println("=" ^ 70)
@printf("Δx = %g m,  U = %g m/s,  dt = %g s,  t = %g s,  C = %.3f\n\n",
        Δx, U_FREESTREAM, dt, t_END, COURANT)

@printf("Pure-upwind reference (modified-equation prediction):\n")
@printf("  D_num = ½·|U|·Δx·(1−C) = %.3f m²/s\n", D_NUM_PRED)
@printf("  σ_pred = √(2·D_num·t)  = %.3f m  (= %.2f Δx)\n\n", SIGMA_PRED, SIGMA_PRED/Δx)

# Part A — pure upwind diagnostic
println("─" ^ 70)
println("Part A: PURE UPWIND (toggle ON, cAlpha=0)")
println("─" ^ 70)
res_A = run_case(true)
σ_A   = fit_sigma(res_A.xs, res_A.αs, x_now; σmax=2.0*SIGMA_PRED)
α_min_A, α_max_A = extrema(res_A.αs)

@printf("\n  σ measured:  %.3f m  (= %.2f Δx)\n", σ_A, σ_A/Δx)
@printf("  σ predicted: %.3f m  (= %.2f Δx)\n", SIGMA_PRED, SIGMA_PRED/Δx)
rel_err_A = abs(σ_A - SIGMA_PRED) / SIGMA_PRED * 100
@printf("  relative error: %.1f %%\n", rel_err_A)
@printf("  α bounds: [%.4f, %.4f]\n", α_min_A, α_max_A)
A_pass = (rel_err_A < 25.0) && (α_max_A < 1.001) && (α_min_A > -0.001)
println("  Part A: ", A_pass ? "PASS" : "FAIL",
        "   (verifies upwind modified-equation prediction)")

# Part B — production MULES + HO
println()
println("─" ^ 70)
println("Part B: MULES + van-Leer HO (toggle OFF, cAlpha=0)")
println("─" ^ 70)
res_B = run_case(false)
σ_B   = fit_sigma(res_B.xs, res_B.αs, x_now; σmax=2.0*SIGMA_PRED)
α_min_B, α_max_B = extrema(res_B.αs)

@printf("\n  σ measured:  %.3f m  (= %.2f Δx)\n", σ_B, σ_B/Δx)
@printf("  σ vs upwind: ratio = %.3f  (HO should give ratio ≪ 1)\n", σ_B/SIGMA_PRED)
@printf("  α bounds: [%.4f, %.4f]\n", α_min_B, α_max_B)
B_pass = (σ_B/SIGMA_PRED < 0.5) && (σ_B/Δx ≤ 4.0) && (α_max_B < 1.001) && (α_min_B > -0.001)
println("  Part B: ", B_pass ? "PASS" : "FAIL",
        "   (verifies HO contribution dominates upwind diffusion)")

# Comparison
println()
println("─" ^ 70)
println("Comparison")
println("─" ^ 70)
@printf("  σ pure-upwind / σ MULES+HO = %.3f  (HO is %.1f× sharper)\n",
        σ_A / max(σ_B, 1e-12), σ_A / max(σ_B, 1e-12))

println()
status = A_pass && B_pass ? "PASS" : "FAIL"
println("Stage 1: $status")

if !A_pass
    println("\nPart A failure interpretation:")
    println("  - σ ≪ predicted: pure-upwind toggle isn't fully active. Check that")
    println("    `alphaf_HO ← alphaf_upwind` actually copies values, and that no")
    println("    other operator path bypasses MULES.")
    println("  - σ ≫ predicted: solver is more diffusive than upwind theory implies.")
    println("    Likely a boundary-condition leak or unintended flux source.")
end
