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








grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_pipe/")
mesh = FOAM3D_mesh(grids_dir, scale=0.001)


backend = CUDABackend(); workgroup = 32
# backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]



gravity = Gravity([0.0, 0.0, -9.81])

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
            Phase(eosModel=ConstEos(1.225), viscosityModel=ConstMu(1.8e-5)),        #vapour
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

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:wall, noSlipVelocity),
            Dirichlet(:inlet_inner, [0.0, 0.0, 0.0]),
            Dirichlet(:inlet_outer, [0.0, 0.0, 0.0]),
            Zerogradient(:outlet_inner),
            Zerogradient(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        p = [
            Wall(:wall),
            Zerogradient(:inlet_inner),
            Zerogradient(:inlet_outer),
            Dirichlet(:outlet_inner, operating_pressure),
            Dirichlet(:outlet_outer, operating_pressure),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        alpha = [
            Wall(:wall),
            Zerogradient(:inlet_inner),
            Zerogradient(:inlet_outer),
            Zerogradient(:outlet_inner),
            Zerogradient(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ]
    )
)

schemes = (
    U =     Schemes(time=Euler, divergence=Upwind),
    p =     Schemes(time=Euler),
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
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jacobi, DILU
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-2
    )
)

runtime = Runtime(
    iterations=200, time_step=1.0e-6, write_interval=50)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()




initialise!(model.momentum.p, operating_pressure)

initialise!(model.fluid.alpha, 0.0)
setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.0, 0.0, 0.0], max_corner=[1.0,1.0,0.5])



# initialise!(model.energy.T, Temp)


residuals = run!(model, config) 