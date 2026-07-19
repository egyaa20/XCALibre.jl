# =============================================================================
# Tier 5 — Energy + variable-density CONTRACTS (Class A)
#
# 5.1 GLOBAL ENERGY BALANCE — the highest bug-catch-per-line contract:
#     a plug flow through a channel with constant wall heat flux q on one wall
#     must, at steady state, satisfy the control-volume statement
#         ṁ · cp · (T̄_out − T_in) = q · A_wall
#     exactly, independent of mesh resolution or profile shape. Catches:
#     cp unit errors (historically 1000×), a missing ρ in the convection flux
#     (~800×), wall-flux BC mis-scaling, and broken T transport — all as
#     order-unity balance violations against a 3 % assert.
#
# Later contracts (activate with the corresponding solver capability):
#     5.2 ṁ·Δh = q·A with the h-formulation + real table
#     5.3 constant-density reduction of the compressible pressure equation
#     5.4 heated-duct continuity ρ̄ŪA = const
#
# Capability-gated: skips (with a message) if the installed multiphase solver
# has no energy stepping.
# =============================================================================

@testset "Tier 5 — energy + variable-density contracts" begin
    if !isdefined(XCALibre.Solvers, :init_multiphase_energy)
        @info "Tier 5 skipped: current multiphase solver is isothermal (no init_multiphase_energy)."
        @test_skip "multiphase energy support absent in installed solver"
    else
        @testset "5.1 global energy balance (ṁ·cp·ΔT = q·A)" begin
            en_mesh = UNV2D_mesh(vv_grid("quad40.unv"), scale=0.001)  # 1 m box
            en_mesh_dev = adapt(CPU(), en_mesh)
            en_hw = Hardware(backend=CPU(), workgroup=AutoTune())

            Uin, T_in, q_wall = 0.1, 300.0, 100.0
            ρ, cp_v, k_v = 1.0, 1000.0, 0.5      # gas-like: ΔT = qA/(ṁcp) = 1 K
            wall_flux(coords, time, fID) = q_wall

            model = Physics(
                time = Transient(),
                fluid = Fluid{Multiphase}(
                    model = VOF(cAlpha=0.0, sigma=0.0),
                    phases = (
                        Phase(rho=ρ, mu=1.0e-5, cp=cp_v, k=k_v),
                        Phase(rho=ρ, mu=1.0e-5, cp=cp_v, k=k_v),  # identical
                    ),
                    gravity = Gravity([0.0, 0.0, 0.0]),
                ),
                turbulence = RANS{Laminar}(),
                energy = Energy{MultiphaseTemperature}(T_init=T_in, Pr_t=0.85),
                domain = en_mesh_dev,
            )

            BCs = assign(
                region = en_mesh_dev,
                (
                    U = [
                        Dirichlet(:inlet, [Uin, 0.0, 0.0]),
                        Zerogradient(:outlet),
                        Zerogradient(:top),      # slip: clean plug flow so the
                        Zerogradient(:bottom),   # outlet mean is flux-weighted
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
                    T = [
                        Dirichlet(:inlet, T_in),
                        Zerogradient(:outlet),
                        Neumann(:top, 0.0),                  # adiabatic
                        NeumannFunction(:bottom, wall_flux), # heated wall
                    ],
                ),
            )

            schemes = (
                U     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
                p     = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
                p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
                alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
                T     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
            )
            solvers = (
                U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                                convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-9),
                p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                                    convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-8, itmax=5000),
                alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                                    convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-9),
                T = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                                convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-10),
            )

            # 4 flow-through times (transit L/U = 10 s) to reach steady state.
            runtime = Runtime(iterations=400, time_step=0.1, write_interval=-1)
            config = Configuration(solvers=solvers, schemes=schemes, runtime=runtime,
                                   hardware=en_hw, boundaries=BCs)

            initialise!(model.fluid.p_rgh, 0.0)
            initialise!(model.momentum.U, [Uin, 0.0, 0.0])
            initialise!(model.fluid.alpha, 1.0)
            initialise!(model.energy.T, T_in)
            GC.gc()
            run!(model, config)

            Tvals = Array(model.energy.T.values)
            @test all(isfinite, Tvals)

            T_out = patch_mean(en_mesh, :outlet, Tvals)
            A_in  = patch_area(en_mesh, :inlet)
            A_hot = patch_area(en_mesh, :bottom)
            mdot  = ρ * Uin * A_in

            convected = mdot * cp_v * (T_out - T_in)
            applied   = q_wall * A_hot
            balance_err = abs(convected - applied) / applied

            @test vv_check("ΔT_out (expected $(round(applied/(mdot*cp_v), digits=3)) K)",
                           T_out - T_in; between=(0.5, 1.5))
            @test vv_check("energy balance |ṁcpΔT − qA|/qA", balance_err; below=0.03)
        end
    end
end
