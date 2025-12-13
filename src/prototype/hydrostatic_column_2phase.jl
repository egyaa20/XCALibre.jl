using XCALibre
using CUDA



# grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_hydrostatic/")
# mesh = FOAM3D_mesh(grids_dir)

scaling = 1.0
# scaling = 0.0025
grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_square_fine/")
# grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_column_fine/")
mesh = FOAM3D_mesh(grids_dir, scale=scaling)

# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "finer_mesh_laplace.unv"
# grid = "laplace_2d_mesh.unv"


# mesh_file = joinpath(grids_dir, grid)

# scaling = 0.025

# mesh = UNV2D_mesh(mesh_file, scale=scaling)

# backend = CUDABackend(); workgroup = 32;
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]


gravity = Gravity([0.0, -9.81, 0.0])
# gravity = Gravity([0.0, 0.0, 0.0])



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
            # Dirichlet(:left, [0.1, 0.0, 0.0]),
            Wall(:left, [0.0, 0.0, 0.0]),
            Wall(:right, [0.0, 0.0, 0.0]),
            # Wall(:top, [0.0, 0.0, 0.0]),
            # Wall(:upper_wall, [0.0, 0.0, 0.0]),
            Zerogradient(:top),
            # Dirichlet(:top, [0.0, -0.001, 0.0]),
            # Zerogradient(:bottom),
            Dirichlet(:bottom, [0.0, 5.0, 0.0]),
            Empty(:frontAndBack)
        ],
        p_rgh = [
            Zerogradient(:left, 0.0), #Symmetry
            Zerogradient(:right, 0.0), #Symmetry
            Zerogradient(:bottom, 0.0), #Zerogradient
            # Zerogradient(:upper_wall), #Zerogradient
            # Zerogradient(:top, 0.0), #Zerogradient
            Dirichlet(:top, 0.0), #Zerogradient
            Empty(:frontAndBack)
        ],
        alpha = [
            # Zerogradient(:left), #Symmetry
            Zerogradient(:left), #Symmetry
            Zerogradient(:right), #Symmetry
            Dirichlet(:bottom, 1.0), #Zerogradient
            Zerogradient(:top), #Zerogradient
            # Dirichlet(:top, 0.0), #Zerogradient
            Empty(:frontAndBack)
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
    p = SolverSetup(
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
    iterations=3000, time_step=1.0e-5, write_interval=50)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()



initialise!(model.momentum.p, 0.0)
initialise!(model.momentum.U, [0.0, 0.0, 0.0])

initialise!(model.fluid.alpha, 0.0)

# min_corner_vec = [-5.0, 0.0, -0.5] * scaling
# max_corner_vec = [5.0,3.0,0.5] * scaling


min_corner_vec = [0.7, 0.0, -0.5] * scaling
max_corner_vec = [1.5,0.5,0.5] * scaling


# setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=1.0, centre=[0.5, 0.5]*scaling, radius=0.15*scaling)
# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=min_corner_vec, max_corner=max_corner_vec)




residuals = run!(model, config)