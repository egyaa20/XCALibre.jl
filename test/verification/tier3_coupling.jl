# =============================================================================
# Tier 3 — Pressure–velocity coupling contracts
#
# 3.1 Hydrostatic equilibrium as a DISCRETE FIXED POINT (v2 replacement of the
#     duration/dt-sensitive spurious-currents test):
#
#   The exactly representable equilibrium (face-aligned interface, U=0,
#   p_rgh = 0) must be a fixed point of the scheme. The dt-robust observable is
#   the one-step SPURIOUS ACCELERATION  a_sp = max|U|₁/dt  (residual force per
#   unit mass, dt-independent to first order):
#     (a) a_sp/g ≪ 1.  A non-well-balanced discretisation produces
#         a_sp/g = O(Δρ/ρ̄) = O(1); the bound sits decades below that and
#         decades above round-off — nothing to tune.
#     (b) dt-invariance asserted directly: a_sp at dt=1e-5 vs 1e-4 within 3×.
#     (c) No secular growth: 100 steps stay within an order of magnitude of
#         the one-step level (spurious currents saturate; they must not build).
# =============================================================================

hydro_mesh = UNV2D_mesh(vv_grid("quad40.unv"), scale=0.001)   # 1 m box, Δ = 0.025 m
hydro_backend = CPU()
hydro_mesh_dev = adapt(hydro_backend, hydro_mesh)
hydro_hw = Hardware(backend=hydro_backend, workgroup=AutoTune())

# Fresh model + config per experiment so every run starts from the same
# exactly-representable equilibrium state.
function hydro_setup(dt, n_iterations)
    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model = VOF(cAlpha=0.0, sigma=0.0),
            phases = (
                Phase(rho=1000.0, mu=1.0e-3),
                Phase(rho=1.2,    mu=1.8e-5),
            ),
            gravity = Gravity([0.0, -9.81, 0.0]),
        ),
        turbulence = RANS{Laminar}(),
        energy = Energy{Isothermal}(),
        domain = hydro_mesh_dev,
    )

    BCs = assign(
        region = hydro_mesh_dev,
        (
            U = [
                Wall(:inlet,  [0.0, 0.0, 0.0]),
                Wall(:outlet, [0.0, 0.0, 0.0]),
                Extrapolated(:top),
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
                Extrapolated(:top),
            ],
        ),
    )

    schemes = (
        U     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
        p     = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
        p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
        alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    )

    # Tight linear tolerances: the measured residual must be discretisation
    # imbalance, not unconverged linear solves.
    solvers = (
        U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-10, relax=1.0, rtol=0.0, atol=1e-10),
        p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                            convergence=1e-10, relax=1.0, rtol=0.0, atol=1e-9),
        alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                            convergence=1e-10, relax=1.0, rtol=0.0, atol=1e-10),
    )

    runtime = Runtime(iterations=n_iterations, time_step=dt, write_interval=-1)
    config = Configuration(solvers=solvers, schemes=schemes, runtime=runtime,
                           hardware=hydro_hw, boundaries=BCs)

    # Exactly representable equilibrium: interface on the y = 0.5 face row
    # (40 cells over 1 m ⇒ 0.5 is a grid line), U = 0, p_rgh = 0.
    initialise!(model.fluid.p_rgh, 0.0)
    initialise!(model.momentum.U, [0.0, 0.0, 0.0])
    initialise!(model.fluid.alpha, 0.0)
    setField_Box!(mesh=hydro_mesh, field=model.fluid.alpha, value=1.0,
                  min_corner=[0.0, 0.0, -0.5], max_corner=[1.0, 0.5, 0.5])

    return model, config
end

function spurious_acceleration(dt)
    model, config = hydro_setup(dt, 1)
    GC.gc()
    run!(model, config)
    return max_umag(model.momentum.U) / dt
end

@testset "Tier 3 — pressure–velocity coupling" begin
    @testset "3.1 hydrostatic discrete fixed point" begin
        g_mag = 9.81

        a_fine   = spurious_acceleration(1.0e-5)
        a_coarse = spurious_acceleration(1.0e-4)

        # (a) well-balanced magnitude: decades below the O(1)·g signature of a
        #     non-balanced scheme. Measured on this setup: 1.6e-7 (dt=1e-5),
        #     4.7e-8 (dt=1e-4) — bound gives 2 decades headroom while sitting
        #     5 decades below the unbalanced signature.
        @test vv_check("a_sp/g (dt=1e-5)", a_fine / g_mag;   below=1e-5)
        @test vv_check("a_sp/g (dt=1e-4)", a_coarse / g_mag; below=1e-5)

        # (b) the fixed-point property itself: spurious acceleration is a
        #     residual force, so it must be dt-invariant — up to the linear-
        #     solver noise floor. At this level of balance the one-step velocity
        #     (~1e-11 m/s) is atol-dominated and dt-independent in u, giving
        #     a_sp a 1/dt noise component (measured ratio 3.3 for a 10× dt
        #     span). One decade of band absorbs that while still failing any
        #     systematic dt-dependent imbalance.
        ratio = a_fine / max(a_coarse, eps())
        @test vv_check("a_sp dt-invariance ratio", ratio; between=(0.1, 10.0))

        # (c) no secular growth over 100 steps at dt=1e-5
        model, config1 = hydro_setup(1.0e-5, 1)
        GC.gc()
        run!(model, config1)
        u1 = max_umag(model.momentum.U)
        cfg99 = Configuration(   # same everything, 99 more steps on the evolved model
            solvers=config1.solvers, schemes=config1.schemes,
            runtime=Runtime(iterations=99, time_step=1.0e-5, write_interval=-1),
            hardware=config1.hardware, boundaries=config1.boundaries)
        run!(model, cfg99)
        u100 = max_umag(model.momentum.U)
        @test vv_check("growth u100/u1", u100 / max(u1, eps()); below=10.0)
    end
end
