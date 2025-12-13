using XCALibre
using CUDA
using LinearAlgebra




# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "laplace_2d_mesh.unv"
# mesh_file = joinpath(grids_dir, grid)
# mesh = UNV2D_mesh(mesh_file) # scale????



grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_injection/")
mesh = FOAM3D_mesh(grids_dir)


# backend = CUDABackend(); workgroup = 32
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]



gravity = Gravity([0.0, -9.81, 0.0])



model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #liquid
            Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #liquid
        ),
        gravity = gravity,
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )



inner_velocity = 1.0
outer_velocity = 0.036

inner_alpha = 1.0
outer_alpha = 0.033


operating_pressure = 0.0



############################################################


BCs = assign(
    region = mesh_dev,
    (
        U = [
            Dirichlet(:inlet, [2.0, 0.0, 0.0]),
            # Wall(:outlet, [0.0, 0.0, 0.0]),
            Zerogradient(:outlet),
            Wall(:walls, [0.0, 0.0, 0.0]),
            Empty(:frontAndBack)
        ],
        p_rgh = [
            Zerogradient(:inlet),
            Dirichlet(:outlet, 0.0),
            # Extrapolated(:walls),
            Zerogradient(:walls, 0.0),
            Empty(:frontAndBack)
        ],
        alpha = [
            Dirichlet(:inlet, 1.0),
            Zerogradient(:outlet),
            Zerogradient(:walls),
            Empty(:frontAndBack)
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
    p = SolverSetup(
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
    iterations=10000, time_step=1.0e-7, write_interval=100)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()




initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, noSlipVelocity) #?????

initialise!(model.fluid.alpha, 0.0)

residuals = run!(model, config) 