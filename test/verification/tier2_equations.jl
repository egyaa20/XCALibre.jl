# =============================================================================
# Tier 2 — Single-equation verification
#
# 2.3 α-step advection contract. Both phases identical (ρ=1, μ=1e-5) so α is a
#     passive tracer in a trivially uniform plug flow — the α equation is
#     exercised in isolation without a frozen-flux solver hook. Contracts:
#       (i)   boundedness:   α ∈ [−ε, 1+ε] with ε at round-off (MULES promise);
#       (ii)  conservation:  Σ α V exact while the front is inside the domain;
#       (iii) kinematics:    the 0.5-crossing of the smeared front sits at
#                             x₀ + U·t within 2 cells (upwind smears the width,
#                             not the centre).
#
# 2.1/2.2 (energy MMS + transient conduction) require the energy equation —
# capability-gated like tier 5, added when energy stepping returns.
# =============================================================================

adv_mesh = UNV2D_mesh(vv_grid("quad40.unv"), scale=0.001)   # 1 m box, Δ = 0.025 m
adv_mesh_dev = adapt(CPU(), adv_mesh)
adv_hw = Hardware(backend=CPU(), workgroup=AutoTune())

@testset "Tier 2 — single-equation verification" begin
    @testset "2.3 α-step advection (boundedness/conservation/kinematics)" begin
        Uin, x0, dt, nsteps = 0.1, 0.25, 0.05, 90        # front: 0.25 → 0.70 m
        Δ = 0.025

        model = Physics(
            time = Transient(),
            fluid = Fluid{Multiphase}(
                model = VOF(cAlpha=0.0, sigma=0.0),
                phases = (
                    Phase(rho=1.0, mu=1.0e-5),
                    Phase(rho=1.0, mu=1.0e-5),   # identical ⇒ α is passive
                ),
                gravity = Gravity([0.0, 0.0, 0.0]),
            ),
            turbulence = RANS{Laminar}(),
            energy = Energy{Isothermal}(),
            domain = adv_mesh_dev,
        )

        BCs = assign(
            region = adv_mesh_dev,
            (
                U = [
                    Dirichlet(:inlet, [Uin, 0.0, 0.0]),
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
                    Dirichlet(:inlet, 0.0),
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
                            convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-9),
            p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                                convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-8),
            alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                                convergence=1e-9, relax=1.0, rtol=0.0, atol=1e-9),
        )
        runtime = Runtime(iterations=nsteps, time_step=dt, write_interval=-1)
        config = Configuration(solvers=solvers, schemes=schemes, runtime=runtime,
                               hardware=adv_hw, boundaries=BCs)

        initialise!(model.fluid.p_rgh, 0.0)
        initialise!(model.momentum.U, [Uin, 0.0, 0.0])
        initialise!(model.fluid.alpha, 0.0)
        setField_Box!(mesh=adv_mesh, field=model.fluid.alpha, value=1.0,
                      min_corner=[0.0, 0.0, -0.5], max_corner=[x0, 1.0, 0.5])

        vols = [c.volume for c in adv_mesh.cells]
        a0 = Array(model.fluid.alpha.values)
        total0 = sum(a0 .* vols)

        GC.gc()
        run!(model, config)

        a = Array(model.fluid.alpha.values)
        @test vv_check("α min", minimum(a); between=(-1e-9, 1.0))
        @test vv_check("α max", maximum(a); between=(0.0, 1.0 + 1e-9))
        @test vv_check("α conservation drift (rel)",
                       abs(sum(a .* vols) - total0) / total0; below=1e-10)

        xs, col = column_means(adv_mesh, a; axis=1)
        xf = front_position(xs, col, 0.5)
        x_expected = x0 + Uin * dt * nsteps
        @test vv_check("front x (expected $(x_expected))", xf;
                       between=(x_expected - 2Δ, x_expected + 2Δ))
    end
end
