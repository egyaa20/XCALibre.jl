using XCALibre
using CUDA
using LinearAlgebra




# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "laplace_2d_mesh.unv"
# mesh_file = joinpath(grids_dir, grid)
# mesh = UNV2D_mesh(mesh_file) # scale????



grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_hydrostatic/")
mesh = FOAM3D_mesh(grids_dir)


# backend = CUDABackend(); workgroup = 32
backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]



gravity = Gravity([0.0, -9.81, 0.0])

# driftVelocity = DriftVelocity(
#             gravity = gravity,
#             d_p = 1.0e-5,
#             drag = Drag_SchillerNaumann()
        # )




model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #liquid
            Phase(eosModel=ConstEos(1.2), viscosityModel=ConstMu(1.8e-5)),        #vapour
            # Phase(eosModel=ConstEos(1.225), viscosityModel=ConstMu(1.8e-5)),        #vapour
            # Phase(eosModel=ConstEos(1000.0), viscosityModel=ConstMu(1.0e-3)),       #liquid
            # Phase(eosModel=ConstEos(1.225), viscosityModel=ConstMu(1.8e-5)),        #vapour
        ),
        gravity = gravity,
        # csf = CSF(sigma=0.072),
        # artificialCompression = ArtificialCompression(compression_coeff=2.0)
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )



inner_velocity = 1.0
outer_velocity = 0.036

inner_alpha = 1.0
outer_alpha = 0.033

# Temp = 250.0

# model.fluid.phases[1].eosModel = HelmholtzEnergy(name=H2(), interpolationMode=true)

# operating_pressure = 3.0e4

operating_pressure = 0.0



############################################################


BCs = assign(
    region = mesh_dev,
    (
        U = [
            # Dirichlet(:left, [0.1, 0.0, 0.0]),
            Wall(:left, [0.0, 0.0, 0.0]),
            Wall(:right, [0.0, 0.0, 0.0]),
            Wall(:top, [0.0, 0.0, 0.0]),
            Wall(:bottom, [0.0, 0.0, 0.0]),
            # Zerogradient(:top), #pressureInletOutletVelocity 0 0 0
            # Dirichlet(:bottom, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Empty(:frontAndBack)
        ],
        p = [
            Zerogradient(:left), #Symmetry
            Zerogradient(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Dirichlet(:top, 0.0), #Zerogradient
            Empty(:frontAndBack)
        ],
        p_rgh = [
            Zerogradient(:left), #Symmetry
            Zerogradient(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Dirichlet(:top, 0.0), #Zerogradient
            Empty(:frontAndBack)
        ],
        alpha = [
            # Zerogradient(:left), #Symmetry
            Zerogradient(:left), #Symmetry
            Zerogradient(:right), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Dirichlet(:top, 0.0),
            Empty(:frontAndBack)
        ]
    )
)


schemes = (
    U =     Schemes(time=Euler, divergence=Upwind),
    p =     Schemes(time=Euler, gradient=Gauss),
    p_rgh =     Schemes(time=Euler, gradient=Gauss),
    alpha = Schemes(time=Euler, divergence=Upwind),
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

# runtime = Runtime(
#     iterations=100000, time_step=0.0001, write_interval=10000)

    
runtime = Runtime(
    iterations=5000, time_step=1.0e-5, write_interval=50)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()




initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, noSlipVelocity) #?????

initialise!(model.fluid.alpha, 0.0)
setField_Circle2D!(mesh=mesh, field=model.fluid.alpha, value=1.0, centre=[0.5, 1.0], radius=0.35)



# initialise!(model.energy.T, Temp)


residuals = run!(model, config) 