using XCALibre
using CUDA




# grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "pipe_coarse_mesh.unv"
# mesh_file = joinpath(grids_dir, grid)
# pipe_fine_mesh.unv
# grid = "pipe_coarse_mesh.unv"

# mesh = UNV2D_mesh(mesh_file)#, scale=0.001)

grids_dir = pkgdir(XCALibre, "src", "prototype", "polyMesh_pipe/")
mesh = FOAM3D_mesh(grids_dir)


backend = CPU(); workgroup = AutoTune(); activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)
mesh_dev = adapt(backend, mesh)

noSlipVelocity = [0.0, 0.0, 0.0]


model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = (
            Phase(eos=PerfectGas(rho=10.0, R=287.0), mu=ConstMu(1.0e-3)),
            Phase(eos=ConstEos(50.0), mu=ConstMu(5.0e-3))
        ),
        gravity = Gravity([0.0, -9.81, 0.0]),
        surfaceTension = ConstSurfaceTension(0.07),
        leeModel = LeeModel(evap_coeff=10.0, condens_coeff=20.0)
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )



inner_velocity = 1.0
outer_velocity = 2.0

inner_alpha = 1.0
outer_alpha = 0.0

operating_pressure = 0.0


BCs = assign(
    region = mesh_dev,
    (
        U = [
            Wall(:wall, noSlipVelocity),
            Dirichlet(:inlet_inner, [0.0, 0.0, inner_velocity]),
            Dirichlet(:inlet_outer, [0.0, 0.0, outer_velocity]),
            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        p = [
            Wall(:wall),
            Extrapolated(:inlet_inner),
            Extrapolated(:inlet_outer),
            Dirichlet(:outlet_inner, operating_pressure),
            Dirichlet(:outlet_outer, operating_pressure),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ],
        alpha = [
            Extrapolated(:wall), # not sure..... Wall(:wall),
            Dirichlet(:inlet_inner, inner_alpha),
            Dirichlet(:inlet_outer, outer_alpha),
            Extrapolated(:outlet_inner),
            Extrapolated(:outlet_outer),
            Symmetry(:sym_1),
            Symmetry(:sym_2)
        ]
    )
)

schemes = (
    U = Schemes(time=Euler, divergence = Upwind),
    p = Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence = Upwind)
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
    iterations=25, time_step=0.1, write_interval=1)
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc()

initialise!(model.momentum.p, 0.0)
initialise!(model.fluid.alpha, 1.0)

 residuals = run!(model, config) 