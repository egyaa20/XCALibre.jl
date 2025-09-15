using XCALibre
using CUDA



grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_dam/")
mesh = FOAM3D_mesh(grids_dir, scale=1.0)


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
            Phase(eosModel=ConstEos(1.225), viscosityModel=ConstMu(1.8e-5)),        #vapour
        ),
        gravity = gravity,
        csf = CSF(sigma=0.07),
        artificialCompression = ArtificialCompression(compression_coeff=2.0)
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )



inner_velocity = 1.0
outer_velocity = 2.0

inner_alpha = 1.0
outer_alpha = 0.0

Temp = 250.0


operating_pressure = 3.0e4

BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:leftWall, noSlipVelocity),
            Wall(:rightWall, noSlipVelocity),
            Wall(:lowerWall, noSlipVelocity),
            Zerogradient(:atmosphere),
            Empty(:defaultFaces)
        ],
        p = [
            Zerogradient(:leftWall),
            Zerogradient(:rightWall),
            Zerogradient(:lowerWall),
            Dirichlet(:atmosphere, 0.0),
            Empty(:defaultFaces)
        ],
        alpha = [
            Zerogradient(:leftWall),
            Zerogradient(:rightWall),
            Zerogradient(:lowerWall),
            Dirichlet(:atmosphere, 0.0),
            Empty(:defaultFaces)
        ]
    )
)

schemes = (
    U = Schemes(time=Euler, divergence=LUST),
    p = Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence=LUST),
)


solvers = (
    U = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres()
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


"""
    setField_Box!(; mesh, field, value::F, min_corner::V, max_corner::V) where {F <: AbstractFloat, V <: AbstractVector}

Sets field values equal to `value` argument only for cells located within a box defined by min and max corners

# Input arguments

- `mesh` reference to `domain` inside a `Physics` model defined by the user.
- `field` field to be modified.
- `value` new values for the chosen field.
- `min_corner` minimum coordinate values of the box in the format of [x, y, z].
- `max_corner` maximum coordinate values of the box in the format of [x, y, z].
"""
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



runtime = Runtime(
    iterations=7000, time_step=5e-6, write_interval=500)
    
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.fluid.alpha, 0.0)
setField_Box!(mesh=mesh, field=model.fluid.alpha, value=1.0, min_corner=[0.0,0.1,-1.0], max_corner=[0.1461,0.292,1.0])


residuals = run!(model, config) 






