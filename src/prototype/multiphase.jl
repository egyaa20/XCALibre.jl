using XCALibre
using CUDA




# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "pipe_coarse_mesh.unv"
# mesh_file = joinpath(grids_dir, grid)
# pipe_fine_mesh.unv
# grid = "pipe_coarse_mesh.unv"

# mesh = UNV2D_mesh(mesh_file)#, scale=0.001)

grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_pipe/")
mesh = FOAM3D_mesh(grids_dir)


backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]


# eos = PengRobinson(T_crit=33.145, p_crit=1.2964e6, omega=-0.216, M=2.016) # values for H2
eos = PengRobinson(T_crit=126.192, p_crit=3.3958e6, omega=0.037, M=28.02) # values for N2


# mu=Sutherland(mu_ref=1.8e-5, S=110.4)
# mu=Andrade(B = 1.732e-6, C = 1863.0)

# gravity = Gravity([0.0, 0.0, -9.81])

gravity = Gravity([0.0, 0.0, -0.98])
driftVelocity = DriftVelocity(
            gravity = gravity,
            d_p = 1.0e-5,
            drag = Drag_SchillerNaumann()
        )


# ConstMu <> 1.0e-3


model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = ( #first phase is liquid, second if vapour - common assumption
            # Phase(eosModel=eos, viscosityModel=ConstMu(1.8e-5)),       #air
            Phase(eosModel=ConstEos(1.225), viscosityModel=ConstMu(1.8e-5)),       #air
            Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #water
            # Phase(eosModel=eos, viscosityModel=ConstMu(1.0e-3))     #water
            # Phase(eosModel=HelmholtzEnergy(name=N2(), interpolationMode=false), viscosityModel=ConstMu(1.8e-5)),     #N2 solver is flawed
            # Phase(eosModel=HelmholtzEnergy(name=N2(), interpolationMode=false), viscosityModel=ConstMu(1.0e-3))    
            # Phase(eosModel=PerfectGas(rho=1.225, R=287.0), viscosityModel=Sutherland(mu_ref=1.8e-5, S=110.4)),       #air
            # Phase(eosModel=ConstEos(1000.0), viscosityModel=Andrade(B = 1.732e-6, C = 1863.0)),       #water
            # Phase(eos=ConstEos(1.225), mu=Sutherland(mu_ref=1.8e-5, S=110.4)),
            # Phase(eos=ConstEos(1.0), mu=ConstMu(1.8e-5)),       #air
            # Phase(eos=PerfectGas(rho=1.225, R=287.0), mu=ConstMu(1.8e-5)),       #air
            # Phase(eos=ConstEos(1000.0), mu=ConstMu(1.0e-3))     #water
            # Phase(eos=ConstEos(1000.0), mu=Andrade(B = 1.732e-6, C = 1863.0))     #water
        ),
        gravity = gravity,
        # driftVelocity = driftVelocity,
        # leeModel = LeeModel(evap_coeff=30.0, condens_coeff=30.0),
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    # energy = Energy{MultiphaseEnergy}(Pr_t=0.85),
    domain = mesh_dev
    )



inner_velocity = 1.0
outer_velocity = 2.0

inner_alpha = 1.0
outer_alpha = 0.0

Temp = 250.0

# model.fluid.phases[1].eosModel = HelmholtzEnergy(name=H2(), interpolationMode=true)

operating_pressure = 3.0e4

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:wall, noSlipVelocity),
            Dirichlet(:inlet_inner, [0.0, 0.0, inner_velocity]),
            Dirichlet(:inlet_outer, [0.0, 0.0, outer_velocity]),
            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        p = [
            Wall(:wall),
            Extrapolated(:inlet_inner),
            Extrapolated(:inlet_outer),
            Dirichlet(:outlet_inner, operating_pressure),
            Dirichlet(:outlet_outer, operating_pressure),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        alpha = [
            Extrapolated(:wall), # not sure..... Wall(:wall),
            Dirichlet(:inlet_inner, inner_alpha),
            Dirichlet(:inlet_outer, outer_alpha),
            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],

        T = [
            Dirichlet(:wall, Temp),
            Dirichlet(:inlet_inner, Temp),
            Dirichlet(:inlet_outer, Temp),
            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ]
    )
)

schemes = (
    U = Schemes(time=Euler, divergence=Upwind),
    p = Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence=Upwind),

    T = Schemes(time=Euler, divergence=Upwind, laplacian=Linear)
)


solvers = (
    U = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # ILU0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    ),
    p = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.2,
        rtol = 1e-3
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    ),
    
    T = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # Jacobi(), #NormDiagonal(), # DILU()
        convergence = 1e-8,
        relax       = 0.8,
        rtol = 1e-4,
        atol = 1e-5
    )
)

runtime = Runtime(
    iterations=100, time_step=1.0e-5, write_interval=25)
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, operating_pressure) # GAUGE vs ABSOLUTE ?
# initialise!(model.momentum.p, 0.1e6)
initialise!(model.fluid.alpha, 1.0)



# initialise!(model.energy.T, Temp)


residuals = run!(model, config) 