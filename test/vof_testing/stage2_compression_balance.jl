# Stage 2 — Interface compression with cAlpha = 1 (single mesh).
#
# Verifies prediction (3) of the modified-equation analysis: with
# compression active (cAlpha = 1), the anti-diffusive flux confined to
# the interface band cancels the leading-order numerical diffusion, so
# the interface settles to a finite width of ~1–2 Δx independent of the
# advected distance.
#
# Same setup as stage 1 (`quad40.unv`, planar α step at x = 500 m,
# U = (1, 0, 0) m/s, t_end = 100 s, uniform fluid) but with cAlpha = 1.
# Interface width is measured as the number of cells in the partial-α
# band 0.05 < α < 0.95, and an erfc fit is also reported for direct
# comparison with stage 1.
#
# Convergence under mesh refinement is verified in stage 3.

using XCALibre
using LinearAlgebra
using Statistics
using Printf

# =============================================================================
# CONFIGURATION (identical to stage 1 except cAlpha)
# =============================================================================

GRID         = "quad40.unv"
Lx           = 1000.0
Δx           = 25.0

U_FREESTREAM = 1.0
X_INTERFACE  = 500.0
SAMPLING_Y   = 500.0

t_END        = 100.0
dt           = 1.0
ITERATIONS   = round(Int, t_END / dt)

COURANT      = U_FREESTREAM * dt / Δx

# Stage 1's predicted spread, for comparison
SIGMA_UPWIND = sqrt(2.0 * 0.5 * U_FREESTREAM * Δx * (1.0 - COURANT) * t_END)

println("=== Stage 2 — Compression balance (cAlpha = 1) ===")
println("Δx           = $Δx m")
println("U            = $U_FREESTREAM m/s")
println("dt           = $dt s")
println("C (Courant)  = $COURANT")
@printf("σ (upwind, cα=0): %.3f m  (= %.2f Δx)  — for comparison\n", SIGMA_UPWIND, SIGMA_UPWIND/Δx)
println("σ predicted (cα=1): ≈ 1.5–2 Δx (steady ~$(round(1.5*Δx, digits=1))–$(round(2*Δx, digits=1)) m)")

# =============================================================================
# CASE
# =============================================================================

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
        model  = VOF(sigma=0.0, cAlpha=1.0, cycles=1, smooth=0),     # ← key: cAlpha = 1
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

@info "Running stage 2: $ITERATIONS steps, t_end = $t_END s, cAlpha = 1.0"
@time residuals = run!(model, config)

# =============================================================================
# POST-PROCESS
# =============================================================================

cells_cpu = Array(mesh.cells)
alpha_cpu = Array(model.fluid.alpha.values)

# Pick the row whose cell-centre y is closest to SAMPLING_Y, then take all
# cells within Δx/4 of that y. Robust to meshes without a centre at
# SAMPLING_Y exactly.
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
xs    = [s[1] for s in sample]
αs    = [s[2] for s in sample]

x_interface_now = X_INTERFACE + U_FREESTREAM * t_END

# Width = number of cells in the partial-α band (0.05 < α < 0.95).
# This is the most direct measure of "interface thickness" for a sharpened front.
band = count(0.05 .< αs .< 0.95)
band_width_cells = max(band, 0)
band_width_m     = band_width_cells * Δx

# Also fit erfc width for direct comparison with stage 1 (will be tiny).
using SpecialFunctions
function erfc_profile(σ, x, x0)
    return 0.5 * erfc((x - x0) / (σ * sqrt(2.0)))
end

σ_grid = collect(range(0.1*Δx, 1.5*SIGMA_UPWIND, length=500))
errs = [sum((erfc_profile(σ, x, x_interface_now) - α)^2 for (x, α) in zip(xs, αs)) for σ in σ_grid]
σ_meas = σ_grid[argmin(errs)]

# Boundedness check — MULES should keep α ∈ [0, 1]
α_max = maximum(αs)
α_min = minimum(αs)

println("\n=== RESULTS ===")
@printf("Interface position (analytical): %.2f m\n", x_interface_now)
@printf("Partial-α band width (0.05<α<0.95): %d cells = %.1f m\n",
        band_width_cells, band_width_m)
@printf("σ measured (erfc fit, for comparison): %.3f m  (= %.2f Δx)\n",
        σ_meas, σ_meas/Δx)
@printf("σ vs upwind reference (stage 1):       %.3f m  → ratio %.3f\n",
        SIGMA_UPWIND, σ_meas/SIGMA_UPWIND)
@printf("α bounds:                       min = %.4f, max = %.4f\n", α_min, α_max)

# Pass criteria:
# (a) interface band ≤ 4 cells (compression should keep it ≲ 2)
# (b) σ ratio (cα=1)/(cα=0) ≪ 1 — compression must visibly sharpen
# (c) MULES keeps α ∈ [-0.001, 1.001] — boundedness preserved
band_ok    = band_width_cells <= 4
ratio_ok   = (σ_meas / SIGMA_UPWIND) < 0.5
bounds_ok  = (α_max < 1.001) && (α_min > -0.001)

println("\n[CRITERION] partial-α band ≤ 4 cells:    ",   band_ok ? "PASS" : "FAIL", "  ($band_width_cells cells)")
println("[CRITERION] σ ratio < 0.5 vs cα=0:       ", ratio_ok ? "PASS" : "FAIL", "  (ratio $(round(σ_meas/SIGMA_UPWIND, digits=3)))")
println("[CRITERION] α ∈ [-0.001, 1.001]:         ", bounds_ok ? "PASS" : "FAIL")

status = band_ok && ratio_ok && bounds_ok ? "PASS" : "FAIL"
println("\nStage 2: $status")
