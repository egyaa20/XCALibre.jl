using XCALibre
using Test
using LinearAlgebra
using CUDA

scaling = 1.0

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
grid = "quad40.unv"
grid = "multiphase_bubble_case.unv"
grid = "multiphase_bubble_case_finer.unv"
# grid = "single_bubble_28k.unv"
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=scaling)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
# backend = CUDABackend(); workgroup = 32

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

gravity = Gravity([0.0, -9.80, 0.0]) # Define gravity direction and magnitude

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(rho=100.0, mu=0.07071),
            Phase(rho=1.0, mu=7.071e-4),
        ),
        gravity = gravity
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )

operating_pressure = 0.0

noSlipVelocity = [0.0, 0.0, 0.0]

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
    alpha = Schemes(time=Euler, gradient=Gauss, divergence=Upwind, laplacian=Linear),
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
        atol        = 1.0e-7,
        itmax       = 5000, 
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
        maxCo=0.5,
        maxAlphaCo=0.6,
        minShrink=0.1,
        maxGrow=1.2
    )

runtime = Runtime(
    iterations=4800, time_step=2.0e-4, write_interval=50)#, adaptive=adaptive)

# 4800 with dt=2e-4 returns the required T=6

    # try finer time step

hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.fluid.p_rgh, 0.0)
initialise!(model.momentum.U, noSlipVelocity)
initialise!(model.fluid.alpha, 1.0)

setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=0.0, centre=[0.5, 0.37], radius=0.25)

@time residuals = run!(model, config, pref=0.0)