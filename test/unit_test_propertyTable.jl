using XCALibre
using Test

@testset "HelmholtzTable property lookup" begin

    # --- gas-phase N2 at 1 bar: table vs direct EOS and NIST anchors ---------
    table = HelmholtzTable(fluid=N2(), p_ref=1.0e5, T_min=200.0, T_max=400.0, n_points=201)
    eos = HelmholtzEnergy(name=N2())

    @testset "table matches direct EOS" begin
        res = eos(300.0, 1.0e5)   # T > Tc: branch 1
        i = findfirst(==(300.0), range(200.0, 400.0, length=201))
        @test table.rho[i] ≈ res.rho[1] rtol=1e-12
        @test table.mu[i]  ≈ res.mu[1]  rtol=1e-12
        @test table.k[i]   ≈ res.k[1]   rtol=1e-12
        @test table.cp[i]  ≈ res.cp[1]*1e3 rtol=1e-12
    end

    @testset "NIST anchors (N2, 300 K, 1 bar)" begin
        i = findfirst(==(300.0), range(200.0, 400.0, length=201))
        @test table.rho[i] ≈ 1.1233   rtol=1e-2
        @test table.mu[i]  ≈ 17.89e-6 rtol=2e-2
        @test table.k[i]   ≈ 0.02573  rtol=0.15
        @test table.cp[i]  ≈ 1041.0   rtol=2e-2
    end

    # --- supercritical N2 at 4 MPa across the pseudo-critical line ----------
    sc = HelmholtzTable(fluid=N2(), p_ref=4.0e6, T_min=90.0, T_max=160.0, n_points=400)

    @testset "supercritical branch continuity" begin
        @test all(isfinite, sc.rho) && all(>(0.0), sc.rho)
        @test issorted(sc.rho, rev=true)          # density falls monotonically with T
        @test issorted(sc.h)                      # h strictly increasing (inversion valid)
        @test sc.rho[1] > 500.0                   # liquid-like at 90 K
        @test sc.rho[end] < 150.0                 # gas-like at 160 K
    end

    # --- field kernels on a mesh --------------------------------------------
    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh = UNV2D_mesh(joinpath(grids_dir, "quad40.unv"), scale=0.001)
    backend = CPU(); workgroup = AutoTune()
    activate_multithread(backend)
    mesh_dev = adapt(backend, mesh)
    config = (hardware=Hardware(backend=backend, workgroup=workgroup),)

    phase = XCALibre.ModelPhysics.build_phase(
        Phase(rho=sc, mu=HelmholtzMu()), mesh_dev)

    T = ScalarField(mesh_dev)
    h = ScalarField(mesh_dev)
    T2 = ScalarField(mesh_dev)

    @testset "update_phase_properties! kernel" begin
        initialise!(T, 130.0)
        update_phase_properties!(phase, T, config)
        res = eos(130.0, 4.0e6)
        @test all(isapprox.(phase.rho.values, res.rho[1]; rtol=1e-4))
        @test all(isapprox.(phase.mu.values,  res.mu[1];  rtol=1e-4))
        @test all(isapprox.(phase.k.values,   res.k[1];   rtol=1e-4))
        @test all(isapprox.(phase.cp.values,  res.cp[1]*1e3; rtol=1e-4))
    end

    @testset "constant-property phase is a no-op" begin
        const_phase = XCALibre.ModelPhysics.build_phase(
            Phase(rho=746.0, mu=1.6e-4), mesh_dev)
        @test update_phase_properties!(const_phase, T, config) === nothing
        @test const_phase.rho[1] == 746.0
    end

    @testset "h <-> T roundtrip" begin
        for Tval in (95.0, 118.7, 128.9, 155.0)   # includes near pseudo-critical
            initialise!(T, Tval)
            enthalpy_from_temperature!(h, T, sc, config)
            temperature_from_enthalpy!(T2, h, sc, config)
            @test all(isapprox.(T2.values, Tval; atol=2*sc.dT))
        end
    end

    @testset "out-of-range T clamps to table edges" begin
        initialise!(T, 80.0)   # below T_min
        update_phase_properties!(phase, T, config)
        @test all(isapprox.(phase.rho.values, sc.rho[1]; rtol=1e-12))
        initialise!(T, 200.0)  # above T_max
        update_phase_properties!(phase, T, config)
        @test all(isapprox.(phase.rho.values, sc.rho[end]; rtol=1e-12))
    end
end
