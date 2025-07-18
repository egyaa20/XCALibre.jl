# TO-DO:

# 1) Construct the required eqns.
# 2) Sort the algorithm
# 3) Possibly fluid properties?



using XCALibre
using CUDA

# backwardFacingStep_2mm, 5mm or 10mm
# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "backwardFacingStep_10mm.unv"


# grids_dir = pkgdir(XCALibre, "prototype")
grids_dir = pkgdir(XCALibre, "src", "prototype", "damBreak_mesh")

# grid = "damBreak_mesh/"

# mesh_file = joinpath(grids_dir, grid)

mesh = FOAM3D_mesh(grids_dir)

# backend = CUDABackend(); workgroup = 32
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

velocity = [0.0, 0.0, 0.0]
nu = 1e-3
Re = velocity[1]*0.1/nu

model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(fluid = :hydrogen),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )




include("setField_utility.jl")

region_to_initialize = "(0 0 -1) (0.1461 0.292 1)"
alpha_val = 1.0

num_modified = setField!(mesh, model.fluid.alpha, alpha_val, region_to_initialize)
# num_modified = setField!(mesh, model.fluid.alpha, alpha_val, region_to_initialize)


# show(IOContext(stdout, :displaysize => (10000, 80)), model.momentum.U.x.values)





BCs = assign(
    region = mesh_dev,
    (
        U = [
            # Dirichlet(:leftWall, [5.0, -5.0, 0.0]),
            Wall(:leftWall, [0.0, 0.0, 0.0]),
            Wall(:rightWall, [0.0, 0.0, 0.0]),
            Wall(:lowerWall, [0.0, 0.0, 0.0]),
            Zerogradient(:atmosphere),
            Empty(:defaultFaces), #?
        ],
        p = [
            Extrapolated(:leftWall),
            Extrapolated(:rightWall),
            Extrapolated(:lowerWall),
            Extrapolated(:atmosphere),
            Empty(:defaultFaces), #?
        ],
        alpha = [
            Zerogradient(:leftWall),
            Zerogradient(:rightWall),
            Zerogradient(:lowerWall),
            Zerogradient(:atmosphere),
            Empty(:defaultFaces), #?
        ]
    )
)

schemes = (
    # U = Schemes(divergence = Linear, limiter=MFaceBased(model.domain)),
    # U = Schemes(divergence = Linear),
    U = Schemes(time=Euler, divergence = Upwind),
    p = Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence = Upwind)
    # p = Schemes(limiter=FaceBased(model.domain))
    # p = Schemes(limiter=MFaceBased(model.domain))
)





# julia> config.schemes
# (U = Schemes(SteadyState, Upwind, Linear, Gauss, nothing), p = Schemes(SteadyState, Linear, Linear, Gauss, nothing))

# julia> config.solvers
# (U = SolverSetup{Float64, Int64, Bicgstab, Nothing, Jacobi}(Bicgstab(), nothing, Jacobi(), 1.0e-7, 0.8, nothing, 1000, 8.161992717227193e-15, 0.01), p = SolverSetup{Float64, Int64, Cg, Nothing, Jacobi}(Cg(), nothing, Jacobi(), 1.0e-7, 0.2, nothing, 1000, 8.161992717227193e-15, 0.001))


# julia> config.schemes.U
# Schemes(SteadyState, Upwind, Linear, Gauss, MFaceBased())




solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # ILU0GPU, Jacobi, DILU
        # smoother=JacobiSmoother(domain=mesh_dev, loops=8, omega=1),
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    ),
    p = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        # smoother=JacobiSmoother(domain=mesh_dev, loops=8, omega=1),
        convergence = 1e-7,
        relax       = 0.2,
        rtol = 1e-3
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        # smoother=JacobiSmoother(domain=mesh_dev, loops=8, omega=1),
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    )
)

runtime = Runtime(
    iterations=20, time_step=0.05, write_interval=1)
    # iterations=1, time_step=1, write_interval=1)

# hardware = Hardware(backend=CUDABackend(), workgroup=32)
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)
# config = adapt(CUDABackend(), config)

GC.gc()

initialise!(model.momentum.p, 0.0)
# initialise!(model.momentum.U, velocity)

# initialise!(model.momentum.U, velocity)
# initialise!(model.momentum.p, 0.0)
# initialise!(model.momentum.p, 1000.0)

@time residuals = run!(model, config) # 1106 iterations!


# Profiling now 
# GC.gc()

# initialise!(model.momentum.p, 100.0)
# initialise!(model.momentum.p, 0.0)

# # @profview residuals = run!(model, config)
# @profview_allocs residuals = run!(model, config) sample_rate=0.00025

# @time residuals = run!(model, config)

# GC.gc()

# initialise!(model.momentum.U, velocity)
# initialise!(model.momentum.p, 0.0)

# @profview residuals = run!(model, config)

# using Plots
# iterations = runtime.iterations
# plot(yscale=:log10, ylims=(1e-9,1e-1))
# plot!(1:iterations, residuals.Ux, label="Ux")
# plot!(1:iterations, residuals.Uy, label="Uy")
# plot!(1:iterations, residuals.Uz, label="Uz")
# plot!(1:iterations, residuals.p, label="p")
# # plot!(1:iterations, residuals.alpha, label="alpha")
