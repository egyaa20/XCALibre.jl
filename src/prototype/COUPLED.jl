# using Revise
# using Adapt
using XCALibre

using CUDA

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "backwardFacingStep_10mm.unv"
# grid = "summer_2d_5x10.unv"
grid = "solid_mesh.unv"
# grid = "summer_3d_extruded_pipe.unv"
mesh_file = joinpath(grids_dir, grid)

mesh = UNV2D_mesh(mesh_file, scale=0.001)
# mesh = UNV3D_mesh(mesh_file, scale=0.001)

backend = CPU(); workgroup = 1024; activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)

mesh_dev_1 = adapt(backend, mesh)


grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "backwardFacingStep_10mm.unv"
grid = "fluid_mesh.unv"
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=0.001)
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev_2 = adapt(backend, mesh)


# model = Physics(
#     time = Transient(),
#     medium = Solid{Coupled}(k=15.0), #add this
#     energy = Energy{CryogenicConduction}(material = :Steel, rho = 8000.0),
#     domain = mesh_dev
#     )

velocity = [1.5, 0.0, 0.0]
nu = 1e-3
Re = velocity[1]*0.1/nu




solid_model = Physics(
    time = Transient(),
    medium = Solid{Uniform}(k = 16.2),
    energy = Energy{CryogenicConduction}(material = :Aluminium, rho = 8000.0),
    domain = mesh_dev_1
    )

fluid_model = Physics(
    time = Transient(),
    medium = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{KOmega}(),
    # turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev_2
    )






BCs_solid = assign(
    region = mesh_dev_1,
    (
        T = [
            Zerogradient(:walls),    
            Zerogradient(:interface)
            # Dirichlet(:walls, 50)      
        ],
    )
)

BCs_fluid = assign(
    region = mesh_dev_2,
    (
        U = [
            Dirichlet(:inlet, velocity),
            Extrapolated(:outlet),
            # Zerogradient(:outlet),
            # Dirichlet(:wall, [0.0, 0.0, 0.0]),
            # Dirichlet(:top, [0.0, 0.0, 0.0]),
            Wall(:wall, [0.0, 0.0, 0.0]),

            # Wall(:top, [0.0, 0.0, 0.0])
            Zerogradient(:interface)
        ],
        p = [
            # Neumann(:inlet, 0.0),
            # Zerogradient(:inlet),
            Extrapolated(:inlet),
            Dirichlet(:outlet, 0.0),
            Wall(:wall),
            # Neumann(:top, 0.0),

            # Wall(:top)
            Zerogradient(:interface)
        ]
    )
)



# mp = MultiPhysics(
#     Coupling1 = Coupling(solid_model, fluid_model, BCs_solid.T[1], BCs_fluid.U[1]) #..[1] needs to be the new interfacing BC! OR does it? I think so...
# #   Coupling2 = Coupling(model3, model4, interface34),
# )

mp = MultiPhysics(
    (
        Coupling1 = Coupling(solid_model, fluid_model, BCs_solid.T[1], BCs_fluid.U[1]),
    )
)

solvers_solid = (
    T = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # Jacobi(), #NormDiagonal(), # DILU()
        convergence = 1e-8,
        relax       = 0.8,
        rtol = 1e-4,
        atol = 1e-5
    )
)


solvers_fluid = (
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
    )
)

schemes_solid = (
    T = Schemes(laplacian = Linear)
)

schemes_fluid = (
    U = Schemes(divergence = Upwind),
    p = Schemes()
)

runtime = Runtime(
    iterations=10, write_interval=1, time_step=0.1) #0.1 * 10 = 1 sec




# It gets tricky here....

config_solid = Configuration(
    solvers=solvers_solid, schemes=schemes_solid, runtime=runtime, hardware=hardware, boundaries=BCs_solid)
    
config_fluid = Configuration(
    solvers=solvers_fluid, schemes=schemes_fluid, runtime=runtime, hardware=hardware, boundaries=BCs_fluid)

configs = [config_solid, config_fluid]

GC.gc(true)

initialise!(solid_model.energy.T, 100.0)



# residuals = run!(model, config)

residuals = multi_run!(mp, configs)