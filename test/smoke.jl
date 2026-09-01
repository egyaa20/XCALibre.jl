# Fast pre-submission gate run by validate.sh on the login node before any
# sbatch. Exercises the APIs the tatsumoto stages depend on: HelmholtzTable at
# the case operating pressure (3.5 MPa, crossing the pseudo-critical line),
# Energy{HelmholtzEnthalpy}, the HeatFlux wall BC, JLD2 checkpointing and a
# short multiphase solve with t_end. Budget: well under a minute after load.

using XCALibre
using JLD2
using Test

include(joinpath(@__DIR__, "..", "cases", "tatsumoto", "stages", "stage_common.jl"))

@testset "smoke" begin
    table = HelmholtzTable(fluid=N2(), p_ref=3.5e6, T_min=90.0, T_max=170.0, n_points=300)
    @test issorted(table.h)
    @test table.rho[1] > 500.0 && table.rho[end] < 150.0
    T0 = 150.0
    @test isapprox(table_temperature(table, table_enthalpy(table, T0)), T0; atol=2*table.dT)

    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh = UNV2D_mesh(joinpath(grids_dir, "quad40.unv"), scale=1.0e-4)
    backend = CPU(); activate_multithread(backend)
    mesh_dev = adapt(backend, mesh)
    hardware = Hardware(backend=backend, workgroup=AutoTune())
    noSlip = [0.0, 0.0, 0.0]

    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model = VOF(cAlpha=0.0, sigma=0.0),
            phases = (Phase(rho=table, mu=HelmholtzMu()), Phase(rho=800.0, mu=1.6e-4)),
            gravity = Gravity([0.0, 0.0, 0.0])),
        turbulence = RANS{Laminar}(),
        energy = Energy{HelmholtzEnthalpy}(),
        domain = mesh_dev)

    BCs = assign(region = mesh_dev, (
        U = [Wall(:inlet, noSlip), Wall(:outlet, noSlip),
             Extrapolated(:top), Wall(:bottom, noSlip)],
        p_rgh = [Zerogradient(:inlet), Zerogradient(:outlet),
                 Dirichlet(:top, 0.0), Zerogradient(:bottom)],
        alpha = [Zerogradient(:inlet), Zerogradient(:outlet),
                 Zerogradient(:top), Zerogradient(:bottom)],
        h = [Zerogradient(:inlet), Zerogradient(:outlet),
             Zerogradient(:top), HeatFluxFunction(:bottom, (c, t, i) -> 40.0)],
    ))
    schemes = (
        U = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
        p = Schemes(time=Euler, gradient=Gauss, laplacian=Linear),
        p_rgh = Schemes(time=Euler, gradient=Gauss, laplacian=Linear),
        alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
        h = Schemes(time=Euler, divergence=Upwind, laplacian=Linear))
    solvers = (
        U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
        p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-8),
        alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
        h = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-8))
    runtime = Runtime(time_step=50.0, write_interval=-1, t_end=500.0,
                      adaptive=AdaptiveTimeStepping(maxCo=0.9, maxAlphaCo=0.85,
                                                    minShrink=0.1, maxGrow=1.1))
    config = Configuration(solvers=solvers, schemes=schemes,
                           runtime=runtime, hardware=hardware, boundaries=BCs)

    initialise!(model.momentum.U, noSlip)
    initialise!(model.fluid.p_rgh, 0.0)
    initialise!(model.fluid.alpha, 1.0)
    initialise!(model.energy.T, T0)

    residuals = run_solver!(model, config; inner_loops=2)

    T = model.energy.T
    @test all(isfinite, T.values)
    @test all(isfinite, model.momentum.U.x.values)
    @test maximum(T.values) > T0 + 0.001   # wall heating arrived
    @test minimum(T.values) > T0 - 0.05

    ckpt = joinpath(mktempdir(), "smoke.jld2")
    save_checkpoint(ckpt; time=1.23,
        U = model.momentum.U, T = model.energy.T, h = model.energy.h)
    T_saved = copy(T.values)
    initialise!(model.energy.T, 100.0)
    t = load_checkpoint!(ckpt; U = model.momentum.U, T = model.energy.T, h = model.energy.h)
    @test t == 1.23
    @test T.values == T_saved
end
println("smoke test passed")
