# =============================================================================
# Tier 4 — Mixture (drift-flux) model validation
#
# 4.1+4.2 combined as one solver-level experiment (Class A — no internals):
#   A quiescent dilute suspension (c₀ = 5 % heavy particles, Re_p ≈ 0.12)
#   settles under gravity. Kinematic (Kynch) theory: the clear-fluid interface
#   descends from the top at the drift speed, which in the Stokes regime is
#       U_St = Δρ·g·d² / (18 μ_c)
#   (Schiller–Naumann correction ≈ 3 % at Re_p = 0.12 — inside tolerance.)
#   Assert: front displacement = U_St·t within 2 cells.
#
# 4.1b Δρ-scaling: repeat with Δρ halved → front displacement ratio ≈ 0.5.
#   Ratio tests cancel common factors, so this isolates the *linearity in Δρ*
#   of the drag law regardless of any shared calibration error.
#
# 4.3 Conservation + boundedness of the dispersed phase ride along.
# =============================================================================

mix_mesh = UNV2D_mesh(vv_grid("quad40.unv"), scale=0.001)   # 1 m box, Δ = 0.025 m
mix_mesh_dev = adapt(CPU(), mix_mesh)
mix_hw = Hardware(backend=CPU(), workgroup=AutoTune())

# Settle a uniform suspension; return (front_y, c_total, α_min, α_max).
function settle_run(; rho_d, d_p, mu_c, c0, dt, nsteps)
    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model = Mixture(diameter=d_p),
            phases = (
                Phase(rho=1000.0, mu=mu_c),      # continuous (α = 1)
                Phase(rho=rho_d,  mu=1.8e-5),    # dispersed heavy phase
            ),
            gravity = Gravity([0.0, -9.81, 0.0]),
        ),
        turbulence = RANS{Laminar}(),
        energy = Energy{Isothermal}(),
        domain = mix_mesh_dev,
    )

    BCs = assign(
        region = mix_mesh_dev,
        (
            U = [
                Wall(:inlet,  [0.0, 0.0, 0.0]),
                Wall(:outlet, [0.0, 0.0, 0.0]),
                Wall(:top,    [0.0, 0.0, 0.0]),
                Wall(:bottom, [0.0, 0.0, 0.0]),
            ],
            p_rgh = [
                Zerogradient(:inlet),
                Zerogradient(:outlet),
                Zerogradient(:bottom),
                Dirichlet(:top, 0.0),
            ],
            alpha = [
                Zerogradient(:inlet),
                Zerogradient(:outlet),
                Zerogradient(:bottom),
                Zerogradient(:top),
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
                        convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-9),
        p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                            convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-8, itmax=5000),
        alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                            convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-9),
    )
    runtime = Runtime(iterations=nsteps, time_step=dt, write_interval=-1)
    config = Configuration(solvers=solvers, schemes=schemes, runtime=runtime,
                           hardware=mix_hw, boundaries=BCs)

    initialise!(model.fluid.p_rgh, 0.0)
    initialise!(model.momentum.U, [0.0, 0.0, 0.0])
    initialise!(model.fluid.alpha, 1.0 - c0)     # uniform suspension
    GC.gc()
    run!(model, config)

    a = Array(model.fluid.alpha.values)
    c = 1.0 .- a
    ys, c_rows = column_means(mix_mesh, c; axis=2)
    front_y = front_position(ys, c_rows, c0 / 2)   # topmost row still holding ≥ c₀/2
    vols = [cl.volume for cl in mix_mesh.cells]
    return front_y, sum(c .* vols), minimum(a), maximum(a)
end

@testset "Tier 4 — mixture (drift-flux) validation" begin
    Δ = 0.025           # cell size
    μ_c, d_p, c0 = 0.5, 3.0e-3, 0.05
    # dt: the mixture path's explicit drift/slip coupling is unstable beyond
    # dt ≈ 0.05 s on this setup (probed: 0.01 stable, 0.05 → NaN). 0.008 gives
    # margin; this dt limit is itself useful documentation of the solver.
    dt = 0.008
    t_end = 8.0
    nsteps = round(Int, t_end / dt)

    # Case 1: Δρ = 2000 → U_St = Δρ g d²/(18 μ_c) ≈ 0.0196 m/s, Re_p ≈ 0.12
    U_st1 = (3000.0 - 1000.0) * 9.81 * d_p^2 / (18μ_c)
    f1, ctot1, amin1, amax1 = settle_run(rho_d=3000.0, d_p=d_p, mu_c=μ_c,
                                          c0=c0, dt=dt, nsteps=nsteps)
    vols = [cl.volume for cl in mix_mesh.cells]
    c_init = c0 * sum(vols)

    @testset "4.1/4.2 settling front vs Stokes" begin
        y_expected = 1.0 - U_st1 * t_end     # Schiller–Naumann corr. ~3 % ⊂ band
        @test vv_check("front y (expected $(round(y_expected, digits=3)))", f1;
                       between=(y_expected - 2Δ, y_expected + 2Δ))
    end

    @testset "4.3 conservation + boundedness" begin
        # History: these two contracts FAILED when first written (2026-07-13):
        # drift 1.0e-3 rel and α_max = 1.0032. Root cause: the drift face factor
        # α(1−α) was upwinded by the MIXTURE flux direction (solver noise in
        # near-quiescent settling), extracting dispersed phase from clear cells
        # and injecting current noise that also drove the ∇·mdotf leak. Fixed by
        # drift-donor upwinding (`drift_donor_alpha!`, donor = sign of Urdotf):
        # drift → 9.3e-11, α_max ≤ 1, spurious currents ↓ 4 decades.
        @test vv_check("dispersed-volume drift (rel)",
                       abs(ctot1 - c_init) / c_init; below=1e-10)
        @test amin1 > -1e-9
        @test amax1 < 1.0 + 1e-9
    end

    @testset "4.1b Δρ linearity" begin
        # Case 2: Δρ halved, time doubled → SAME expected displacement.
        # Ratio target 1.0 separates linear-in-Δρ (1.0) from quadratic (0.5)
        # and Δρ-independent (2.0), with half the front-quantisation noise of
        # a direct half-displacement comparison.
        f2, _, _, _ = settle_run(rho_d=2000.0, d_p=d_p, mu_c=μ_c,
                                  c0=c0, dt=dt, nsteps=2nsteps)
        disp_ratio = (1.0 - f2) / (1.0 - f1)
        @test vv_check("equal-displacement Δρ-linearity ratio", disp_ratio;
                       between=(0.75, 1.3))
    end
end
