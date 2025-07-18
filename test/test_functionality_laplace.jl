using XCALibre

using CUDA

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")

grid = "finer_mesh_laplace.unv"
# grid = "laplace_2d_mesh.unv"

mesh_file = joinpath(grids_dir, grid)

mesh = UNV2D_mesh(mesh_file, scale=0.001)

backend = CPU(); workgroup = 1024; activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)

mesh_dev = adapt(backend, mesh)


model = Physics(
    time = Steady(),
    solid = Solid{Uniform}(k=54.0), #Steady state requires only k; Transient needs k, cp, rho
    energy = Energy{LaplaceEnergy}(),
    # energy = Energy{CryogenicConduction}(material = :Steel, rho = 7850.0), #We still need to define cp and rho in Solid{}
    # CryogenicConduction crashes for Temps below 4K, so won't work for 0 and 1 test
    # On this geometry,the default (LaplaceEnergy) yields straight line division between two regions
    # CryogenicConduction however shows non-linear line; it's nonlinearity increases at lower temperatures
    #                   thus indicating that k and cp vary increasingly more at cryo temps => thus I assume it works well
    domain = mesh_dev
    )


BCs = assign(
    region = mesh_dev,
    (
        T = [     
            Dirichlet(:left_wall, 200.0),
            Dirichlet(:right_wall, 100.0),
            Zerogradient(:bottom_wall),
            Zerogradient(:upper_wall)
        ],
    )
)




    # EXTRAPOLATED BC ACTS VERY WEIRD!!!!!!
# BCs = assign(
#     region = mesh_dev,
#     (
#         T = [     
#             Dirichlet(:left_wall, 0),
#             Extrapolated(:right_wall),
#             Dirichlet(:bottom_wall, 1),
#             Extrapolated(:upper_wall)
#         ],
#     )
# )



solvers = (
    T = SolverSetup(
        solver      = Gmres(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # Jacobi(), #NormDiagonal(), # DILU()
        convergence = 1e-16,
        relax       = 0.5,
        rtol = 1e-10,
        atol = 1e-10
    )
)

schemes = (
    T = Schemes(laplacian = Linear)
)


runtime = Runtime(
    iterations=100, write_interval=100, time_step=1)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc(true)

initialise!(model.energy.T, 290)

residuals = run!(model, config)

# Results can be seen in figures "grid_laplace.png" and "testing_laplace.png"