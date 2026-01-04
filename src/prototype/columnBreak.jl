using XCALibre
using CUDA

# grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_dam/")
# mesh = FOAM3D_mesh(grids_dir, scale=1.0)

# scaling = 0.00001 # make sure the domain is 1x1 m * 0.01
scaling = 0.001 # make sure the domain is 1x1 m

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "quad40.unv"
grid = "quad100.unv"
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=scaling)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
# backend = CUDABackend(); workgroup=32

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]

gravity = Gravity([0.0, -9.81, 0.0]) # Define gravity direction and magnitude

AndradeModel = Andrade(B=2.414e-5, C=247.8) #water
SutherlandModel = Sutherland(mu_ref=1.716e-5, S=111.0) #air
PerfectGasModel = PerfectGas(rho=1.225, R=287.0) #air

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #liquid
            Phase(eosModel=ConstEos(1.2), viscosityModel=ConstMu(1.8e-5)),          #vapour
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
            Wall(:inlet, noSlipVelocity),
            Wall(:outlet, noSlipVelocity),
            Zerogradient(:top),
            Wall(:bottom, noSlipVelocity),
        ],
        p_rgh = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:bottom),
            Dirichlet(:top, operating_pressure),
        ],
        alpha = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:bottom),
            Zerogradient(:top, 0.0),
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
    iterations=10000, time_step=1.0e-4, write_interval=50)
    # iterations=30000, time_step=1.0e-9, write_interval=500)
     
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, noSlipVelocity)
initialise!(model.fluid.alpha, 0.0)

min_corner_vec = [0.0, 0.0, -0.5] # column
max_corner_vec = [0.3,0.4,0.5] # column




# setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=1.0, centre=[0.005,0.005], radius=0.0015)
setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=min_corner_vec, max_corner=max_corner_vec)

residuals = run!(model, config)