using XCALibre
using CUDA


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


BCs = assign(
    region = mesh_dev,
    (
        U = [
            # Dirichlet(:left, [0.1, 0.0, 0.0]),
            Symmetry(:left),
            Symmetry(:right), #slip
            # Zerogradient(:top), #pressureInletOutletVelocity 0 0 0
            Dirichlet(:top, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Dirichlet(:bottom, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Empty(:frontAndBack)
        ],
        p_rgh = [
            # Zerogradient(:left), #Symmetry
            Symmetry(:left), #Symmetry
            Symmetry(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Dirichlet(:top, 0.0), #Zerogradient
            # Zerogradient(:top), #Zerogradient
            Empty(:frontAndBack)
        ],
        alpha = [
            # Zerogradient(:left), #Symmetry
            Symmetry(:left), #Symmetry
            Symmetry(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Zerogradient(:top),
            Empty(:frontAndBack)
        ]
    )
)


schemes = (
    U =     Schemes(time=Euler, divergence=Upwind),
    p =     Schemes(time=Euler, gradient=Gauss),
    p_rgh = Schemes(time=Euler, gradient=Gauss),
    alpha = Schemes(time=Euler, divergence=Upwind),
)


solvers = (
    U = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # ILU0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    ),
    p_rgh = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.2,
        rtol = 1e-3
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    )
)

runtime = Runtime(
    iterations=250, time_step=1.0e-5, write_interval=10)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, noSlipVelocity) #?????

initialise!(model.fluid.alpha, 0.0)
setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[-5.0, 0.0, -0.5], max_corner=[5.0,1.0,0.5])



residuals = run!(model, config) 