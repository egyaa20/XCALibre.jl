# using Revise
# using Adapt
using XCALibre

using CUDA

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "backwardFacingStep_10mm.unv"
# grid = "summer_2d_5x10.unv"
grid = "finer_mesh_laplace.unv"
grid = "laplace_2d_mesh.unv"
# grid = "summer_3d_extruded_pipe.unv"
mesh_file = joinpath(grids_dir, grid)

mesh = UNV2D_mesh(mesh_file, scale=0.001)
# mesh = UNV3D_mesh(mesh_file, scale=0.001)

backend = CPU(); workgroup = 1024; activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)

mesh_dev = adapt(backend, mesh)


model = Physics(
    time = Transient(),
    solid = Solid{Uniform}(k=10.0),
    energy = Energy{LaplaceEnergy}(),
    domain = mesh_dev
    )

BCs = assign(
    region = mesh_dev,
    (
        T = [     
            Dirichlet(:left_wall, 0),
            Zerogradient(:right_wall),
            Dirichlet(:bottom_wall, 1),
            Zerogradient(:upper_wall)
            # Zerogradient(:outlet),    
            # Zerogradient(:walls)
            # Dirichlet(:walls, 50)      
        ],
    )
)


# solvers = (
#     T = SolverSetup(
#         solver      = Cg(), # Bicgstab(), Gmres()
#         preconditioner = Jacobi(), # Jacobi(), #NormDiagonal(), # DILU()
#         convergence = 1e-8,
#         relax       = 0.8,
#         rtol = 1e-4,
#         atol = 1e-5
#     )
# )

solvers = (
    T = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = DILU(), # Jacobi(), #NormDiagonal(), # DILU()
        convergence = 1e-8,
        relax       = 0.8,
        rtol = 1e-4,
        atol = 1e-5
    )
)

schemes = (
    T = Schemes(laplacian = Linear)
)


runtime = Runtime(
    iterations=50, write_interval=1, time_step=0.1) #0.1 * 100 = 10 sec

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc(true)

initialise!(model.energy.T, 100.0)

residuals = run!(model, config)