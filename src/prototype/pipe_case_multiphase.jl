using XCALibre
using CUDA
grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_pipe")
mesh = FOAM3D_mesh(grids_dir, scale=0.001)


backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

velocity = [0.0, 0.0, 0.0]




### QUESTION : How to turn Const() into just "=10.0" ?
model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            liquid(name=Water(), transport=ConstMu(5.1e-3)),
            gas(name=Air(), eos=PerfectGas(1.0), transport=Sutherland()),
        )
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )


# fluid = Fluid{Multiphase}(
#     phases = (
#         Fluid(name=N2(), eos=Helmholtz()) #takes care about liquid / vapour on its own
#         #OR
#         Liquid(name=Water(), eos=Const(), mu=Const()),    #alpha = 1
#         Gas(name=Air(), eos=PerfectGas(), mu=Sutherland())     #alpha = 0
#         IF NAME IS NOT PASSED, then more info inside eos, mu, etc is required (e.g. custom fluid)
#     ),
    
#     physicsProperties = (
#         gravity=Gravity(axis=Y(), g=9.81), #Axis, magnitude (-ve by default)
#         slipVelocity(...) #define 'd' etc
#         interphase(...)
#     )
# )



U_inner = 1.0
U_outer = 0.036

alpha_inner = 1.0
alpha_outer = 0.0
# alpha_outer = 0.91667

operating_pressure = 0.0 # OR 0.0; please initialise at this value!
# operating_pressure = 100000.0 # OR 0.0; please initialise at this value!

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:wall, [0.0, 0.0, 0.0]), #[0.0, 0.0, 0.0]

            Dirichlet(:inlet_inner, [0.0, 0.0, U_inner]),
            Dirichlet(:inlet_outer, [0.0, 0.0, U_outer]),

            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),

            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        p = [
            Extrapolated(:wall),  
            # Wall(:wall), #[0.0, 0.0, 0.0]

            Extrapolated(:inlet_inner),
            Extrapolated(:inlet_outer),

            Dirichlet(:outlet_inner, operating_pressure),
            Dirichlet(:outlet_outer, operating_pressure),

            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        alpha = [
            Extrapolated(:wall),  
            # Wall(:wall), #[0.0, 0.0, 0.0]

            Dirichlet(:inlet_inner, alpha_inner),
            Dirichlet(:inlet_outer, alpha_outer),

            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),

            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ]
    )
)

schemes = (
    U = Schemes(time=Euler, divergence = BoundedUpwind),
    p = Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence = BoundedUpwind)
)







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
    iterations=1, time_step=1.0e-6, write_interval=100)
    # iterations=1, time_step=1, write_interval=1)

# hardware = Hardware(backend=CUDABackend(), workgroup=32)
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)
# config = adapt(CUDABackend(), config)

GC.gc()

initialise!(model.momentum.p, operating_pressure)

initialise!(model.fluid.alpha, 1.0) # Not sure

residuals = run!(model, config)

# println(model.fluid.alpha.values)


# region_to_initialize = "(0 0 -1) (0.05 0.1 1)"
# alpha_val = 0.0

# num_modified = setField!(mesh, model.fluid.alpha, alpha_val, region_to_initialize)


# initialise!(model.momentum.U, velocity)

# initialise!(model.momentum.U, velocity)
# initialise!(model.momentum.p, 0.0)
# initialise!(model.momentum.p, 1000.0)

# residuals = run!(model, config) # 1106 iterations!

# show(IOContext(stdout, :displaysize => (10000, 80)), model.fluid.alpha.values)



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
