using XCALibre
using CUDA

# grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_dam/")
# mesh = FOAM3D_mesh(grids_dir, scale=1.0)

# scaling = 0.00001 # make sure the domain is 1x1 m * 0.01
scaling = 1.0 # make sure the domain is 1x1 m

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
grid = "damBreak_coarse.unv"
# grid = "damBreak_medium.unv"
# grid = "quad100.unv"[]
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=scaling)
    
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
# backend = CUDABackend(); workgroup=32

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]

gravity = Gravity([0.0, -9.81, 0.0]) # Define gravity direction and magnitude


model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(density=1000.0, mu=1.0e-3),       #liquid
            Phase(density=1.2, mu=1.8e-5),          #vapour
        ),
        gravity = gravity
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )

operating_pressure = 0.0


BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:leftWall, noSlipVelocity),
            Wall(:rightWall, noSlipVelocity),
            Zerogradient(:upperWall), 
            Wall(:lowerWall, noSlipVelocity),
            # Wall(:top, noSlipVelocity),
        ],
        p_rgh = [
            Zerogradient(:leftWall),
            Zerogradient(:rightWall),
            Zerogradient(:lowerWall),
            # Zerogradient(:top),
            Dirichlet(:upperWall, 0.0),
            # totalPressure(:top, 0.0),
        ],
        alpha = [
            Zerogradient(:leftWall),
            Zerogradient(:rightWall),
            Zerogradient(:lowerWall),
            Zerogradient(:upperWall),
            # Dirichlet(:top, 0.0),
        ]
    )
)


schemes = (
    U =     Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    p =     Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
    alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
)


solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # ILU0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol        = 0.0,
        atol        = 1.0e-5
    ),
    p_rgh = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol        = 0.0,
        atol        = 1.0e-5
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol        = 0.0,
        atol        = 1.0e-5
    )
)

runtime = Runtime(
    iterations=14000, time_step=1.0e-4, write_interval=500)
    # iterations=35000, time_step=2.5e-5, write_interval=500)
     
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, noSlipVelocity)
initialise!(model.fluid.alpha, 0.0)

min_corner_vec = [0.0, 0.0, -0.5] # column
max_corner_vec = [0.6, 0.3, 0.5] # column


setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=min_corner_vec, max_corner=max_corner_vec)

residuals = run!(model, config)