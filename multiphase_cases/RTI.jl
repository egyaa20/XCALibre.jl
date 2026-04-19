# This is the two-phase hydrostatic column under gravitational effect test case, it verifies velocity field coupling and convergence

using XCALibre
using Test
using LinearAlgebra
using CUDA

# scaling = 0.001
scaling = 1.0

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
grid = "RTI_coarse.unv"
# grid = "RTI_medium.unv"
# grid = "multiphase_bubble_case_xfine.unv"
mesh_file = joinpath(grids_dir, grid)
mesh = UNV2D_mesh(mesh_file, scale=scaling)

backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)
# backend = CUDABackend(); workgroup = 32

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]

# gravity = Gravity([0.0, -9.81, 0.0]) # Define gravity direction and magnitude

# API CHANGE from:
        # phases = (
        #     Phase(rho=1000.0, mu=1.0e-3),
        #     Phase(rho=1.2, mu=1.8e-5),
        # ),
# To:
        # phases = (
        #     water = Phase(rho=1000.0, mu=1.0e-3),
        #     air = Phase(rho=1.2, mu=1.8e-5),
        #     alpha = :water
        # ),

gravity = Gravity([0.0, -1.0, 0.0]) # Define gravity direction and magnitude
# gravity = Gravity([0.0, 0.0, 0.0]) # Define gravity direction and magnitude
model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(rho=1.8, mu=0.001),
            Phase(rho=1.0, mu=0.001)
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
            Wall(:inlet, noSlipVelocity),
            Wall(:outlet, noSlipVelocity),
            Zerogradient(:top), 
            Wall(:bottom, noSlipVelocity),
        ],
        p_rgh = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:bottom),
            Zerogradient(:top),
            # Dirichlet(:top, 0.0),
        ],
        alpha = [
            Zerogradient(:inlet),
            Zerogradient(:outlet),
            Zerogradient(:bottom),
            Zerogradient(:top),
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
        atol        = 1.0e-7
        
    ),
    alpha = SolverSetup(
        solver      = Bicgstab(), # Bicgstab(), Gmres(), Cg()
        preconditioner = Jacobi(), # IC0GPU, Jaco•bi, DILU
        convergence = 1e-7,
        relax       = 1.0,
        rtol        = 0.0,
        atol        = 1.0e-5
    )
)

adaptive = AdaptiveTimeStepping(
        maxCo=0.6,
        maxAlphaCo=0.4,
        minShrink=0.1,
        maxGrow=1.2
    )

runtime = Runtime(
    iterations=5000, time_step=1.0e-3, write_interval=50, adaptive=adaptive)

hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.fluid.p_rgh, 0.0)
initialise!(model.momentum.U, noSlipVelocity)
initialise!(model.fluid.alpha, 0.0)

# Domain dimensions
Lx     = 1.0
Ly     = 2.0

# With 5000 equal cells on a 1×2 domain:
# Equal cell size → dx = dy, so nx/ny = Lx/Ly = 1/2
# nx * ny = 5000 → nx = 50, ny = 100
n_cols = 50
dx     = Lx / n_cols

# Interface parameters — y = 1.0 + 0.1·cos(2π·x)
A      = 0.1
λ      = Lx          # wavelength = domain width (single cosine period)
z_min  = -0.5
z_max  =  0.5

for i in 1:n_cols
    x_left  = (i - 1) * dx
    x_right =  i      * dx
    x_mid   = 0.5 * (x_left + x_right)

    # Interface position for this column
    y_interface = 1.0 + A * cos(2π * x_mid / λ)
    y_interface = clamp(y_interface, 0.0, Ly)

    # Set heavy fluid (alpha = 1) above the interface
    setField_Box!(
        mesh       = mesh,
        field      = model.fluid.alpha,
        value      = 1.0,
        min_corner = [x_left,  y_interface, z_min],
        max_corner = [x_right, Ly,          z_max]
    )
end

@time residuals = run!(model, config, pref=0.0)