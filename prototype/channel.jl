
# Step 0. Load libraries
using XCALibre
# using CUDA # Uncomment to run on NVIDIA GPUs
# using AMDGPU # Uncomment to run on AMD GPUs


grids_dir = pkgdir(XCALibre, "examples/testing_grids")
# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "messy_3d.unv"
# grid = "channel.unv"
grid = "cylinder_3d_netgen.unv"
# grid = "unstructured_2d.unv"


#Define 2 grids:
# 2D Unstructured
# 3D Unstructured


mesh_file = joinpath(grids_dir, grid)

# Step 1. Define mesh
# mesh_file = "channel.unv"

mesh = UNV3D_mesh(mesh_file, scale=0.001) # convert mesh

# Step 2. Select backend and setup hardware
backend = CPU()
# backend = CUDABackend() # ru non NVIDIA GPUs
# backend = ROCBackend() # run on AMD GPUs

hardware = set_hardware(backend=backend, workgroup=1024)
# hardware = set_hardware(backend=backend, workgroup=32) # use for GPU backends

mesh_dev = mesh # use this line to run on CPU
# mesh_dev = adapt(backend, mesh)  # Uncomment to run on GPU 

# Step 3. Flow conditions
velocity = [5, 0.0, 0.0]
nu = 1e-3
Re = velocity[1]*0.05/nu

# Step 4. Define physics
model = Physics(
    time = Steady(),
    fluid = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )

# Step 5. Define boundary conditions
@assign! model momentum U (
    Dirichlet(:inlet, velocity),
    Neumann(:outlet, 0.0),
    Wall(:cylinder, 0.0),
    Neumann(:bottom, 0.0),
    Neumann(:front, 0.0),
    Neumann(:back, 0.0),
    Neumann(:top, 0.0)
)

@assign! model momentum p (
    Neumann(:inlet, 0.0),
    Dirichlet(:outlet, 0.0),
    Neumann(:cylinder, 0.0),
    Neumann(:bottom, 0.0),
    Neumann(:front, 0.0),
    Neumann(:back, 0.0),
    Neumann(:top, 0.0)
)

# Step 6. Choose discretisation schemes
schemes = (
    U = set_schemes(divergence = Upwind),
    p = set_schemes() # no input provided (will use defaults)
)

# Step 7. Set up linear solvers and preconditioners
solvers = (
    U = set_solver(
        model.momentum.U;
        solver      = BicgstabSolver, # Options: GmresSolver
        preconditioner = Jacobi(), # Options: NormDiagonal()
        convergence = 1e-7,
        relax       = 0.7,
        rtol = 1e-1,
    ),
    p = set_solver(
        model.momentum.p;
        solver      = CgSolver, # Options: CgSolver, BicgstabSolver, GmresSolver
        preconditioner = Jacobi(), # Options: NormDiagonal()
        convergence = 1e-7,
        relax       = 0.3,
        rtol = 1e-2
    )
)

# Step 8. Specify runtime requirements
runtime = set_runtime(iterations=2000, time_step=1, write_interval=500)
# runtime = set_runtime(iterations=1, time_step=1, write_interval=-1) # hide

# Step 9. Construct Configuration object
config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware)

# Step 10. Initialise fields (initial guess)
initialise!(model.momentum.U, velocity)
initialise!(model.momentum.p, 0.0)

# Step 11. Run simulation
residuals = run!(model, config); #ncorrectors=4

# Step 12. Post-process
pwd() # find active directory where the file "iteration_002000.vtk" was saved