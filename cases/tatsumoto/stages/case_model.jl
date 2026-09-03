# =============================================================================
# Shared model/BC/numerics setup for the tatsumoto supercritical case.
# Include after stage_common.jl. Solver: Solvers_5_Multiphase.jl with
# HelmholtzTable properties, Energy{HelmholtzEnthalpy} and HeatFlux wall BCs.
# =============================================================================

function build_backend(hw_cfg)
    if uppercase(hw_cfg["backend"]) == "CUDA"
        error("CUDA backend not wired for this case yet — set [hardware] backend = \"CPU\"")
    end
    backend = CPU()
    activate_multithread(backend)
    return backend, Hardware(backend = backend, workgroup = AutoTune())
end

"""
    build_case(CFG, mesh_dev; wall_heat_bc, write_interval, t_end)

Physics model, boundary conditions and solver configuration for the
supercritical N2 pipe. `wall_heat_bc` is the h-equation BC applied on
Wall_Heated (HeatFlux / HeatFluxFunction). Returns (model, config, htable).
"""
function build_case(CFG, mesh_dev, hardware; wall_heat_bc, write_interval, t_end)
    flow_cfg   = CFG["flow"];   thermo  = CFG["thermo"]
    energy_cfg = CFG["energy"]; run_cfg = CFG["run"]
    adapt_cfg  = run_cfg["adaptive"]

    U_IN  = Float64(flow_cfg["U_in"])
    T_IN  = Float64(flow_cfg["T_in"])
    P_OP  = Float64(thermo["pressure"])
    PR_T  = Float64(energy_cfg["Pr_t"])

    TI    = Float64(flow_cfg["turb_intensity"])
    D_HYD = Float64(flow_cfg["D_hyd"])
    ELL   = 0.07 * D_HYD
    K_INLET     = 1.5 * (TI * U_IN)^2
    OMEGA_INLET = sqrt(K_INLET) / (ELL * 0.09^0.25)
    NUT_INLET   = K_INLET / OMEGA_INLET

    inletVelocity  = [0.0, 0.0, U_IN]
    noSlipVelocity = [0.0, 0.0, 0.0]
    gravity        = Gravity([0.0, 0.0, -9.81])

    htable = HelmholtzTable(
        fluid    = N2(),
        p_ref    = P_OP,
        T_min    = Float64(thermo["T_table_lo"]),
        T_max    = Float64(thermo["T_table_hi"]),
        n_points = Int(thermo["n_points"]),
        h_ref_T  = T_IN,
    )
    H_IN = table_enthalpy(htable, T_IN)

    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model = VOF(cAlpha = 0.0, sigma = 0.0),
            phases = (
                Phase(rho = htable, mu = HelmholtzMu()),
                Phase(rho = 800.0, mu = 1.6e-4),   # inert secondary, alpha ≡ 1
            ),
            gravity = gravity,
        ),
        turbulence = RANS{KOmegaSST}(walls = (:Wall_Heated, :Wall_Unheated)),
        energy = Energy{HelmholtzEnthalpy}(Pr_t = PR_T),
        domain = mesh_dev,
    )

    BCs = assign(
        region = mesh_dev,
        (
            U = [
                Dirichlet(:Inlet, inletVelocity),
                Wall(:Wall_Heated,   noSlipVelocity),
                Wall(:Wall_Unheated, noSlipVelocity),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            p_rgh = [
                Zerogradient(:Inlet),
                Zerogradient(:Wall_Heated),
                Zerogradient(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Dirichlet(:Outlet, 0.0),
            ],
            alpha = [
                Dirichlet(:Inlet, 1.0),
                Zerogradient(:Wall_Heated),
                Zerogradient(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            h = [
                Dirichlet(:Inlet, H_IN),
                wall_heat_bc,
                Zerogradient(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            k = [
                Dirichlet(:Inlet, K_INLET),
                Dirichlet(:Wall_Heated, 0.0),
                Dirichlet(:Wall_Unheated, 0.0),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            omega = [
                Dirichlet(:Inlet, OMEGA_INLET),
                OmegaWallFunction(:Wall_Heated),
                OmegaWallFunction(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            nut = [
                Extrapolated(:Inlet),
                Dirichlet(:Wall_Heated, 0.0),
                Dirichlet(:Wall_Unheated, 0.0),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
        )
    )

    schemes = (
        U     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
        p     = Schemes(time = Euler, gradient   = Gauss,  laplacian = Linear),
        p_rgh = Schemes(time = Euler, gradient   = Gauss,  laplacian = Linear),
        alpha = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
        h     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
        k     = Schemes(divergence = Upwind),
        omega = Schemes(divergence = Upwind),
        y     = Schemes(),
    )

    solvers = (
        U = SolverSetup(
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.7, rtol = 0.0, atol = 1.0e-5,
        ),
        p_rgh = SolverSetup(
            solver = Cg(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.3, rtol = 0.0, atol = 1.0e-6, itmax = 25000,
        ),
        alpha = SolverSetup(
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 1.0, rtol = 0.0, atol = 1.0e-5,
        ),
        h = SolverSetup(    # enthalpy ~1e5 J/kg, so atol is scaled up
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-7, relax = 0.7, rtol = 1.0e-6, atol = 1.0,
        ),
        k = SolverSetup(
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-10, relax = 0.6, rtol = 1e-3,
        ),
        omega = SolverSetup(
            solver = Bicgstab(), preconditioner = Jacobi(),
            convergence = 1e-10, relax = 0.6, rtol = 1e-3,
        ),
        y = SolverSetup(
            solver = Cg(), preconditioner = Jacobi(),
            convergence = 1e-10, relax = 0.7, rtol = 1e-5, itmax = 5000,
        ),
    )

    adaptive = AdaptiveTimeStepping(
        maxCo      = Float64(adapt_cfg["maxCo"]),
        maxAlphaCo = Float64(adapt_cfg["maxAlphaCo"]),
        minShrink  = Float64(adapt_cfg["minShrink"]),
        maxGrow    = Float64(adapt_cfg["maxGrow"]),
    )

    runtime = Runtime(
        time_step      = Float64(run_cfg["time_step"]),
        write_interval = Int(write_interval),
        adaptive       = adaptive,
        t_end          = Float64(t_end),
    )

    config = Configuration(
        solvers = solvers, schemes = schemes, runtime = runtime,
        hardware = hardware, boundaries = BCs,
    )

    inits = (U = inletVelocity, T = T_IN, k = K_INLET,
             omega = OMEGA_INLET, nut = NUT_INLET)

    return model, config, htable, inits
end

function initialise_fields!(model, inits)
    initialise!(model.fluid.p_rgh,      0.0)
    initialise!(model.momentum.U,       inits.U)
    initialise!(model.fluid.alpha,      1.0)
    initialise!(model.energy.T,         inits.T)
    initialise!(model.turbulence.k,     inits.k)
    initialise!(model.turbulence.omega, inits.omega)
    initialise!(model.turbulence.nut,   inits.nut)
    return nothing
end

checkpoint_fields(model) = (
    U     = model.momentum.U,
    p_rgh = model.fluid.p_rgh,
    alpha = model.fluid.alpha,
    T     = model.energy.T,
    h     = model.energy.h,
    k     = model.turbulence.k,
    omega = model.turbulence.omega,
    nut   = model.turbulence.nut,
)
