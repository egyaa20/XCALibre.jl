# This test case is the mini dam-break which ensures correct development of alpha field with time under the influence of gravity

using XCALibre
using Test

scaling = 0.001*0.05

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "quad40.unv"
grid = "quad100.unv"
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=scaling)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]

gravity = Gravity([0.0, -9.81, 0.0]) # Define gravity direction and magnitude

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            # rho !
            Phase(density=1000.0, mu=1.0e-3), # Higher density phase always comes first
            Phase(density=1.2, mu=1.8e-5),
        ),
        # solver = VOF, # dispatch thingy
        # explicit = false,
        #explicit = true by default; otherwise explicict = false
        gravity = gravity
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )

# if typeof ===  free check for performance

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
            Dirichlet(:top, 0.0),
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
    iterations=5000, time_step=5.0e-5, write_interval=50)

hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, noSlipVelocity)
initialise!(model.fluid.alpha, 0.0)

min_corner_vec = [0.0, 0.0, -0.5] # column
max_corner_vec = [0.35, 0.6, 0.5] # column (was 0.3, 0.4 previously)

# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=min_corner_vec, max_corner=max_corner_vec)
setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=1.0, centre=[0.5*0.05, 0.5*0.05], radius=0.15*0.05)

alpha_avg_before = boundary_average(:bottom, model.fluid.alpha, config.boundaries.alpha, config)
@time residuals = run!(model, config)
alpha_avg_after = boundary_average(:bottom, model.fluid.alpha, config.boundaries.alpha, config)
alpha_top = boundary_average(:top, model.fluid.alpha, config.boundaries.alpha, config)

@test isapprox(alpha_avg_before, 0.3; rtol=0.01)
@test alpha_avg_after > 0.75
@test alpha_top < 1.0e-4

println(alpha_avg_after)


#SURFACE TENSION TEST : BUBBLE