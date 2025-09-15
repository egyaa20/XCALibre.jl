using XCALibre
using CUDA




# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "pipe_coarse_mesh.unv"
# mesh_file = joinpath(grids_dir, grid)
# pipe_fine_mesh.unv
# grid = "pipe_coarse_mesh.unv"

# mesh = UNV2D_mesh(mesh_file)#, scale=0.001)

grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_pipe/")
mesh = FOAM3D_mesh(grids_dir, scale=0.001)


# backend = CUDABackend(); workgroup = 32
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]



gravity = Gravity([0.0, 0.0, -9.81])

driftVelocity = DriftVelocity(
            gravity = gravity,
            d_p = 1.0e-5,
            drag = Drag_SchillerNaumann()
        )

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #liquid
            Phase(eosModel=ConstEos(1.225), viscosityModel=ConstMu(1.8e-5)),        #vapour
        ),
        gravity = gravity,
        # csf = CSF(sigma=0.07),
        # artificialCompression = ArtificialCompression(compression_coeff=2.0)
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )



inner_velocity = 1.0
outer_velocity = 0.04

inner_alpha = 1.0
outer_alpha = 0.0

Temp = 250.0

# model.fluid.phases[1].eosModel = HelmholtzEnergy(name=H2(), interpolationMode=true)

operating_pressure = 3.0e4
# operating_pressure = 0.0

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:wall, noSlipVelocity),
            Dirichlet(:inlet_inner, [0.0, 0.0, inner_velocity]),
            Dirichlet(:inlet_outer, [0.0, 0.0, outer_velocity]),
            Zerogradient(:outlet_inner),
            Zerogradient(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        p = [
            Wall(:wall, 0.0),
            Zerogradient(:inlet_inner),
            Zerogradient(:inlet_outer),
            Dirichlet(:outlet_inner, 0.0),
            Dirichlet(:outlet_outer, 0.0),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        alpha = [
            Wall(:wall, 0.0), # not sure..... Wall(:wall),
            Dirichlet(:inlet_inner, inner_alpha),
            Dirichlet(:inlet_outer, outer_alpha),
            Zerogradient(:outlet_inner),
            Zerogradient(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ]#,

        # T = [
        #     Dirichlet(:wall, Temp),
        #     Dirichlet(:inlet_inner, Temp),
        #     Dirichlet(:inlet_outer, Temp),
        #     Extrapolated(:outlet_inner),
        #     Extrapolated(:outlet_outer),
        #     Symmetry(:sym_1),
        #     Symmetry(:sym_2)
        # ]
    )
)

schemes = (
    U =     Schemes(time=Euler, divergence=Upwind),
    p =     Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence=Upwind),

    # T = Schemes(time=Euler, divergence=Upwind, laplacian=Linear)
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
        rtol = 1e-3
    )#,
    
    # T = SolverSetup(
    #     solver      = Cg(), # Bicgstab(), Gmres()
    #     preconditioner = Jacobi(), # Jacobi(), #NormDiagonal(), # DILU()
    #     convergence = 1e-8,
    #     relax       = 0.8,
    #     rtol = 1e-4,
    #     atol = 1e-5
    # )
)

runtime = Runtime(
    iterations=2000, time_step=0.00001, write_interval=500)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

# initialise!(model.momentum.p, operating_pressure) # GAUGE vs ABSOLUTE ?
initialise!(model.momentum.p, 0.0)
# initialise!(model.fluid.alpha, 1.0)



# initialise!(model.energy.T, Temp)


residuals = run!(model, config) 