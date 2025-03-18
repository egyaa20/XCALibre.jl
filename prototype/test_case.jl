
using XCALibre
# using CUDA # uncomment to run on GPU

# grids_dir = pkgdir(XCALibre, "examples\\testing_grids\\3D")
# grid = "3D" 

grids_dir = pkgdir(XCALibre, "examples/testing_grids")
# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "messy_3d.unv"
grid = "4by4by4.unv"
# grid = "unstructured_2d.unv"


#Define 2 grids:
# 2D Unstructured
# 3D Unstructured


mesh_file = joinpath(grids_dir, grid)

mesh = UNV3D_mesh(mesh_file, scale=0.001)

sum_volumes = 0.0

for cell ∈ mesh.cells
    sum_volumes += cell.volume

    faces_pointer = cell.faces_range
    nodes_pointer = cell.nodes_range

    normal_signs_array = mesh.cell_nsign[faces_pointer]  #indices to access nsigns of the cell
    faces_array = mesh.cell_faces[faces_pointer] #indices to access faces of the cell
    nodes_array = mesh.cell_nodes[faces_pointer] #indices to access nodes of the cell


    # cell.volume: sum up and check against real solution
    # cell.centre: 
    # face.area
    # face.centre
    # face.weight
    # face.delta
    # face.ownerCells
    # face.normal
end

sum_volumes
analytical_volume = 0.1*0.1*0.1
error_volume = 100*(analytical_volume-sum_volumes)/analytical_volume

for face ∈ mesh.faces
    sum_volumes
end



mesh_dev = mesh
# mesh_dev = adapt(CUDABackend(), mesh) # uncomment to run on GPU

# Inlet conditions
velocity = [5, 0.0, 0.0]
noSlip = [0.0, 0.0, 0.0]
nu = 1e-3
Re = (0.2*velocity[1])/nu
display(Re)

model = Physics(
    time = Steady(),
    fluid = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(), 
    domain = mesh_dev
    )

@assign! model momentum U ( 
    Dirichlet(:inlet, velocity),
    Neumann(:outlet, 0.0),
    Wall(:cylinder, noSlip),
    Neumann(:bottom, 0.0),
    Neumann(:top, 0.0)
    # Neumann(:front, 0.0),
    # Neumann(:back, 0.0)
)

@assign! model momentum p (
    Neumann(:inlet, 0.0),
    Dirichlet(:outlet, 0.0),
    Neumann(:cylinder, 0.0),
    Neumann(:bottom, 0.0),
    Neumann(:top, 0.0)
    # Neumann(:front, 0.0),
    # Neumann(:back, 0.0)
    
)

solvers = (
    U = set_solver(
        model.momentum.U;
        solver      = BicgstabSolver, # BicgstabSolver, GmresSolver
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax       = 1.0,
        rtol = 1e-4,
        atol = 1e-5
    ),
    p = set_solver(
        model.momentum.p;
        solver      = CgSolver, # BicgstabSolver, GmresSolver
        preconditioner = Jacobi(), #NormDiagonal(),
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-4,
        atol = 1e-5
    )
)

schemes = (
    U = set_schemes(divergence=Upwind, gradient=Midpoint),
    p = set_schemes(gradient=Midpoint)
    
)


runtime = set_runtime(iterations=500, write_interval=20, time_step=1) 


hardware = set_hardware(backend=CPU(), workgroup=1024)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware)

GC.gc(true)

initialise!(model.momentum.U, velocity)
initialise!(model.momentum.p, 0.0)

residuals = run!(model, config)