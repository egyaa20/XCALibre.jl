using XCALibre
using Test

erfc_x = collect(0.0:0.25:3.0)
erfc_v = [1.0, 0.723674, 0.479500, 0.288845, 0.157299, 0.077100, 0.033895,
          0.013346, 0.004678, 0.001463, 0.000407, 0.000101, 0.0000221]
function erfc_lerp(x)
    x >= 3.0 && return 0.0
    i = min(floor(Int, x/0.25) + 1, length(erfc_x) - 1)
    w = (x - erfc_x[i]) / 0.25
    erfc_v[i]*(1 - w) + erfc_v[i+1]*w
end

function build_sc_model(mesh_dev, table; gvec=[0.0, 0.0, 0.0])
    Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model = VOF(cAlpha=0.0, sigma=0.0),
            phases = (
                Phase(rho=table, mu=HelmholtzMu()),
                Phase(rho=800.0, mu=1.6e-4),
            ),
            gravity = Gravity(gvec)
        ),
        turbulence = RANS{Laminar}(),
        energy = Energy{HelmholtzEnthalpy}(),
        domain = mesh_dev
    )
end

sc_schemes = (
    U     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    p     = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    h     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
)

sc_solvers = (
    U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                    convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
    p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                    convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-8),
    alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                    convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-6),
    h = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                    convergence=1e-7, relax=1.0, rtol=0.0, atol=1.0e-8),
)

@testset "Multiphase Helmholtz enthalpy energy" begin

    table = HelmholtzTable(fluid=N2(), p_ref=4.0e6, T_min=90.0, T_max=170.0, n_points=400)

    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh = UNV2D_mesh(joinpath(grids_dir, "quad40.unv"), scale=0.001)
    backend = CPU(); workgroup = AutoTune()
    activate_multithread(backend)
    mesh_dev = adapt(backend, mesh)
    hardware = Hardware(backend=backend, workgroup=workgroup)
    noSlip = [0.0, 0.0, 0.0]

    T0 = 150.0
    cell_volumes = [c.volume for c in mesh.cells]
    yc = [c.centre[2] for c in mesh.cells]

    @testset "A: thermal equilibrium is preserved" begin
        model = build_sc_model(mesh_dev, table)
        BCs = assign(region = mesh_dev, (
            U = [Wall(:inlet, noSlip), Wall(:outlet, noSlip),
                 Wall(:top, noSlip), Wall(:bottom, noSlip)],
            p_rgh = [Zerogradient(:inlet), Zerogradient(:outlet),
                     Dirichlet(:top, 0.0), Zerogradient(:bottom)],
            alpha = [Zerogradient(:inlet), Zerogradient(:outlet),
                     Zerogradient(:top), Zerogradient(:bottom)],
            h = [Zerogradient(:inlet), Zerogradient(:outlet),
                 Zerogradient(:top), Zerogradient(:bottom)],
        ))
        runtime = Runtime(iterations=50, time_step=5.0e-3, write_interval=-1)
        config = Configuration(solvers=sc_solvers, schemes=sc_schemes,
                               runtime=runtime, hardware=hardware, boundaries=BCs)
        initialise!(model.momentum.U, noSlip)
        initialise!(model.fluid.p_rgh, 0.0)
        initialise!(model.fluid.alpha, 1.0)
        initialise!(model.energy.T, T0)

        run!(model, config)

        U = model.momentum.U
        vel = sqrt.(U.x.values.^2 .+ U.y.values.^2 .+ U.z.values.^2)
        rho_exact = XCALibre.ModelPhysics._table_lerp(
            table.rho, T0, table.T_min, table.dT, length(table.rho))

        @test maximum(vel) < 1e-7
        @test all(isapprox.(model.energy.T.values, T0; atol=0.05))
        @test all(isapprox.(model.fluid.rho.values, rho_exact; rtol=1e-3))
    end

    @testset "B: transient conduction from heated wall matches erfc profile" begin
        T_wall = 152.0
        t_end = 1.0e4
        mesh_B = UNV2D_mesh(joinpath(grids_dir, "quad40.unv"), scale=1.0e-4)
        mesh_B_dev = adapt(backend, mesh_B)
        vol_B = [c.volume for c in mesh_B.cells]
        yc_B = [c.centre[2] for c in mesh_B.cells]
        model = build_sc_model(mesh_B_dev, table)
        BCs = assign(region = mesh_B_dev, (
            U = [Wall(:inlet, noSlip), Wall(:outlet, noSlip),
                 Extrapolated(:top), Wall(:bottom, noSlip)],
            p_rgh = [Zerogradient(:inlet), Zerogradient(:outlet),
                     Dirichlet(:top, 0.0), Zerogradient(:bottom)],
            alpha = [Zerogradient(:inlet), Zerogradient(:outlet),
                     Zerogradient(:top), Zerogradient(:bottom)],
            h = [Zerogradient(:inlet), Zerogradient(:outlet),
                 Zerogradient(:top), Dirichlet(:bottom, table_enthalpy(table, T_wall))],
        ))
        runtime = Runtime(iterations=100, time_step=t_end/100, write_interval=-1)
        config = Configuration(solvers=sc_solvers, schemes=sc_schemes,
                               runtime=runtime, hardware=hardware, boundaries=BCs)
        initialise!(model.momentum.U, noSlip)
        initialise!(model.fluid.p_rgh, 0.0)
        initialise!(model.fluid.alpha, 1.0)
        initialise!(model.energy.T, T0)

        rho0 = XCALibre.ModelPhysics._table_lerp(
            table.rho, T0, table.T_min, table.dT, length(table.rho))
        mass0 = rho0 * sum(vol_B)
        run!(model, config)
        mass1 = sum(model.fluid.rho.values .* vol_B)

        T = model.energy.T
        @test all(isfinite, T.values)
        @test minimum(T.values) > T0 - 0.05
        @test maximum(T.values) < T_wall + 0.05
        @test mass1 < mass0

        Tmid = 0.5*(T0 + T_wall)
        i_mid = argmin(abs.(collect(range(90.0, 170.0, length=400)) .- Tmid))
        rho_m = table.rho[i_mid]; cp_m = table.cp[i_mid]; k_m = table.k[i_mid]
        a = k_m / (rho_m * cp_m)

        for i in eachindex(yc_B)
            eta = yc_B[i] / (2*sqrt(a*t_end))
            theta_exact = erfc_lerp(eta)
            theta = (T.values[i] - T0) / (T_wall - T0)
            @test isapprox(theta, theta_exact; atol=0.12)
        end
    end

    @testset "C: hydrostatic well-balance with table-driven density" begin
        model = build_sc_model(mesh_dev, table; gvec=[0.0, -9.81, 0.0])
        BCs = assign(region = mesh_dev, (
            U = [Wall(:inlet, noSlip), Wall(:outlet, noSlip),
                 Wall(:top, noSlip), Wall(:bottom, noSlip)],
            p_rgh = [Zerogradient(:inlet), Zerogradient(:outlet),
                     Dirichlet(:top, 0.0), Zerogradient(:bottom)],
            alpha = [Zerogradient(:inlet), Zerogradient(:outlet),
                     Zerogradient(:top), Zerogradient(:bottom)],
            h = [Zerogradient(:inlet), Zerogradient(:outlet),
                 Zerogradient(:top), Zerogradient(:bottom)],
        ))
        runtime = Runtime(iterations=100, time_step=5.0e-3, write_interval=-1)
        config = Configuration(solvers=sc_solvers, schemes=sc_schemes,
                               runtime=runtime, hardware=hardware, boundaries=BCs)
        initialise!(model.momentum.U, noSlip)
        initialise!(model.fluid.p_rgh, 0.0)
        initialise!(model.fluid.alpha, 1.0)
        initialise!(model.energy.T, T0)

        run!(model, config)

        U = model.momentum.U
        vel = sqrt.(U.x.values.^2 .+ U.y.values.^2 .+ U.z.values.^2)
        rho_exact = XCALibre.ModelPhysics._table_lerp(
            table.rho, T0, table.T_min, table.dT, length(table.rho))
        p = model.momentum.p

        @test maximum(vel) < 1e-6
        for i in eachindex(yc)
            @test isapprox(p.values[i], -rho_exact*9.81*yc[i]; rtol=1e-2)
        end
    end
end
