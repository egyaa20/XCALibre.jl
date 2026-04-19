using XCALibre
using Test
using LinearAlgebra
using CUDA

# scaling = 0.001
scaling = 1.0

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
grid = "fine_square.unv"
# grid = "multiphase_bubble_case_xfine.unv"
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=scaling)

# backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
backend = CUDABackend(); workgroup = 32

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]

gravity = Gravity([0.0, -9.81, 0.0]) # Define gravity direction and magnitude

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            water = Phase(rho=1000.0, mu=0.156),
            air = Phase(rho=100.0, mu=0.078),
            alpha = :water
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
            Zerogradient(:top),
            # Dirichlet(:top, 0.0),
        ],
        alpha = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:bottom),
            Zerogradient(:top),
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
        atol        = 1.0e-7
        
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jaco•bi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol        = 0.0,
        atol        = 1.0e-5
    )
)

adaptive = AdaptiveTimeStepping(
        maxCo=0.75,
        maxAlphaCo=0.5,
        minShrink=0.1,
        maxGrow=1.2
    )

runtime = Runtime(
    iterations=700, time_step=5.0e-4, write_interval=1)

hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.fluid.p_rgh, 0.0)
initialise!(model.momentum.U, noSlipVelocity)
initialise!(model.fluid.alpha, 1.0)

min_corner_vec = [0.0, 0.0, -0.5]
max_corner_vec = [1.0,0.5,0.5]

# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=min_corner_vec, max_corner=max_corner_vec)
setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=0.0, centre=[0.5, 0.65], radius=0.15)
setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=0.0, centre=[0.5, 0.35], radius=0.1)

@time residuals = run!(model, config, pref=0.0)
# @time residuals = run!(model, config)
