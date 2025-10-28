using XCALibre
using CUDA

# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# # grid = "laplace_2d_mesh.unv"
# grid = "finer_mesh_laplace.unv"
# mesh_file = joinpath(grids_dir, grid)
# mesh = UNV2D_mesh(mesh_file)


grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_hydrostatic/")
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
            Phase(eosModel=ConstEos(1.2), viscosityModel=ConstMu(1.8e-5)),        #vapour
        ),
        gravity = gravity
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )


operating_pressure = 0.0


############################################################


# BCs = assign(
#     region = mesh_dev,
#     (
#         U = [
#             Zerogradient(:left_wall),
#             # Dirichlet(:left_wall, [0.001, 0.0, 0.0]),
#             Zerogradient(:right_wall), #slip
#             # Zerogradient(:top), #pressureInletOutletVelocity 0 0 0
#             Zerogradient(:upper_wall), #dirichlet 0 0 0
#             Wall(:bottom_wall, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
#             # Empty(:frontAndBack)
#         ],
#         p_rgh = [
#             # Zerogradient(:left), #Symmetry
#             Zerogradient(:left_wall), #Symmetry
#             Zerogradient(:right_wall), #Symmetry
#             Zerogradient(:bottom_wall), #Zerogradient
#             Dirichlet(:upper_wall, 0.0), #Zerogradient
#             # Zerogradient(:upper_wall), #Zerogradient
#             # Empty(:frontAndBack)
#         ],
#         alpha = [
#             # Zerogradient(:left), #Symmetry
#             Zerogradient(:left_wall), #Symmetry
#             Zerogradient(:right_wall), #Symmetry
#             Zerogradient(:bottom_wall), #Zerogradient
#             Dirichlet(:upper_wall, 0.0),
#             # Empty(:frontAndBack)
#         ]
#     )
# )







BCs = assign(
    region = mesh_dev,
    (
        U = [
            # Zerogradient(:left),
            # Dirichlet(:left_wall, [0.001, 0.0, 0.0]),
            # Zerogradient(:right), #slip
            # Zerogradient(:top), #pressureInletOutletVelocity 0 0 0
            # Zerogradient(:top), #dirichlet 0 0 0
            Wall(:left, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Wall(:right, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Zerogradient(:top), #dirichlet 0 0 0
            Wall(:bottom, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Empty(:frontAndBack)
        ],
        p_rgh = [
            # Zerogradient(:left), #Symmetry
            Zerogradient(:left), #Symmetry
            Zerogradient(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Dirichlet(:top, 11.76), #Zerogradient
            # Zerogradient(:top), #Zerogradient
            # Zerogradient(:upper_wall), #Zerogradient
            Empty(:frontAndBack)
        ],
        alpha = [
            # Zerogradient(:left), #Symmetry
            Zerogradient(:left), #Symmetry
            Zerogradient(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Zerogradient(:top), #Zerogradient
            Empty(:frontAndBack)
        ]
    )
)




schemes = (
    U =     Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    p =     Schemes(time=Euler, gradient=Midpoint, laplacian=Linear),
    p_rgh = Schemes(time=Euler, gradient=Midpoint, laplacian=Linear),
    alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
)


solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # ILU0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol = 1e-2
    ),
    p_rgh = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol = 1e-7
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol = 1e-2
    )
)

runtime = Runtime(
    iterations=1000, time_step=1.0e-5, write_interval=10)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, 0.0)
initialise!(model.momentum.U, [0.0, 0.0, 0.0]) #?????

initialise!(model.fluid.alpha, 0.0)


setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[-5.0, 0.0, -0.5], max_corner=[5.0,0.25,0.5])

# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.1, 0.25, -0.5], max_corner=[5.0,0.3,0.5])
# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.15, 0.3, -0.5], max_corner=[5.0,0.35,0.5])
# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.2, 0.35, -0.5], max_corner=[5.0,0.4,0.5])
# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.25, 0.4, -0.5], max_corner=[5.0,0.45,0.5])
# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.3, 0.45, -0.5], max_corner=[5.0,0.5,0.5])




# residuals = run!(model, config, pref=0.0)
residuals = run!(model, config)