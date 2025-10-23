using XCALibre
using CUDA

function setField_Box!(; mesh, field, value::F, min_corner::V, max_corner::V) where {F <: AbstractFloat, V <: AbstractVector}

    @assert length(min_corner) == 3 "`min_corner` must have exactly 3 elements"
    @assert length(max_corner) == 3 "`max_corner` must have exactly 3 elements"

    cells_in_region = Int[]
    
    for (id, cell) in enumerate(mesh.cells) #check that X, Y, Z coords are within the box
        center = cell.centre
        if (min_corner[1] <= center[1] <= max_corner[1] &&
            min_corner[2] <= center[2] <= max_corner[2] &&
            min_corner[3] <= center[3] <= max_corner[3])
            
            push!(cells_in_region, id)
            # Warning: if outer cell boundary is outside the box but its centre is within the box, this cell will count too
        end
    end

    if !isempty(cells_in_region)
        field.values[cells_in_region] .= value # Reassign values inside the field
    end
    
    return length(cells_in_region)
end




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
            Dirichlet(:left, [0.001, 0.0, 0.0]),
            Zerogradient(:right), #slip
            Zerogradient(:top), #pressureInletOutletVelocity 0 0 0
            Wall(:bottom, [0.0, 0.0, 0.0]), #dirichlet 0 0 0
            Empty(:frontAndBack)
        ],
        # p = [
        #     # Zerogradient(:left), #Symmetry
        #     Symmetry(:left), #Symmetry
        #     Symmetry(:right), #Symmetry
        #     Zerogradient(:bottom), #Zerogradient
        #     # Dirichlet(:top, 0.0), #Zerogradient
        #     Zerogradient(:top), #Zerogradient
        #     Empty(:frontAndBack)
        # ],
        p_rgh = [
            # Zerogradient(:left), #Symmetry
            Zerogradient(:left), #Symmetry
            Dirichlet(:right, 0.0), #Symmetry
            Zerogradient(:bottom), #Zerogradient
            Zerogradient(:top), #Zerogradient
            # Zerogradient(:top), #Zerogradient
            Empty(:frontAndBack)
        ],
        alpha = [
            # Zerogradient(:left), #Symmetry
            Dirichlet(:left, 1.0), #Symmetry
            Dirichlet(:right, 1.0), #Symmetry
            Dirichlet(:bottom, 1.0), #Zerogradient
            Dirichlet(:top, 1.0),
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
        relax       = 0.8,
        rtol = 1e-2
    ),
    p = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.2,
        rtol = 1e-3
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
    iterations=50, time_step=0.0001, write_interval=1)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()




initialise!(model.momentum.p, operating_pressure)
initialise!(model.momentum.U, [0.0001, 0.0, 0.0]) #?????

initialise!(model.fluid.alpha, 1.0)
# setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[-5.0, 0.0, -0.5], max_corner=[5.0,1.0,0.5])



# initialise!(model.energy.T, Temp)


residuals = run!(model, config) 