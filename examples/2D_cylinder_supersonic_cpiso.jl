using XCALibre
using CUDA # Run this if using NVIDIA GPU
# using AMDGPU # Run this if using AMD GPU

# Mach 1.2 flow over a cylinder using the pressure-based transient compressible
# solver (CPISO). A detached bow shock forms ahead of the cylinder.
#
# Notes for pressure-based shock capturing:
#  - `p` MUST use `divergence = Upwind`. The compressible pressure equation carries a
#    pressure-convection (transonic-correction) term; central differencing (`Linear`,
#    the default) oscillates and diverges at shocks. Upwind is essential here.
#  - Solver `limit` clamps on p and he keep the solution bounded through the strong
#    startup transient (standard practice for pressure-based compressible solvers).
#  - This is the pressure-based route; for strongly supersonic flow the density-based
#    `SupersonicFlow`/Godunov solver (see 2D_cylinder_supersonic.jl) is more robust.

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
mesh_file = joinpath(grids_dir, "cylinder_d10mm_25mm.unv")
mesh = UNV2D_mesh(mesh_file, scale=0.001)

# backend = CPU(); workgroup = 1024
backend = CUDABackend(); workgroup = 32
mesh_dev = adapt(backend, mesh)

# Freestream conditions at Mach 1.2
gamma = 1.4
cp = 1005.0
Pr = 0.7
nu = 1e-5
R = cp*(1.0 - 1.0/gamma)
cv = cp - R
T_inf = 300.0
p_inf = 101325.0
Tref = 0.0
Mach = 1.2
U_inf = Mach*sqrt(gamma*R*T_inf)
velocity = [U_inf, 0.0, 0.0]

model = Physics(
    time = Transient(),
    fluid = Fluid{Compressible}(nu=nu, cp=cp, gamma=gamma, Pr=Pr),
    turbulence = RANS{Laminar}(),
    # energy = Energy{SensibleEnthalpy}(Tref=Tref),
    energy = Energy{InternalEnergy}(Tref=Tref),
    domain = mesh_dev
    )

boundaries = assign(
    region = mesh,
    (
        U = [
            Dirichlet(:inlet, velocity),
            Zerogradient(:outlet),
            Slip(:cylinder),
            Slip(:top),
            Slip(:bottom)
        ],
        p = [
            Dirichlet(:inlet, p_inf),
            Zerogradient(:outlet),
            Zerogradient(:cylinder),
            Zerogradient(:top),
            Zerogradient(:bottom)
        ],
        he = [
            # FixedTemperature(:inlet, T=T_inf, Enthalpy(cp=cp, Tref=Tref)),
            FixedTemperature(:inlet, T=T_inf, IEnergy(cv=cv, Tref=Tref)),
            Zerogradient(:outlet),
            Zerogradient(:cylinder),
            Zerogradient(:top),
            Zerogradient(:bottom)
        ]
    )
)

solvers = (
    U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(), convergence=1e-7, relax=0.7, rtol=1e-2),
    p = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(), convergence=1e-7, relax=0.3, limit=(0.02*p_inf, 50*p_inf), rtol=1e-2),
    he = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(), convergence=1e-7, relax=0.4, limit=(50.0, 6000.0), rtol=1e-2)
)

schemes = (
    U = Schemes(divergence=Upwind, gradient=Gauss, time=Euler),
    p = Schemes(divergence=Upwind, gradient=Gauss, time=Euler),  # Upwind is required (see notes)
    he = Schemes(divergence=Upwind, gradient=Gauss, time=Euler)
)

runtime = Runtime(iterations=20000, write_interval=100, time_step=5e-8)
hardware = Hardware(backend=backend, workgroup=workgroup)

config = Configuration(;
    solvers, schemes, runtime, hardware, boundaries)

GC.gc(true)

initialise!(model.momentum.U, velocity)
initialise!(model.momentum.p, p_inf)
initialise!(model.energy.T, T_inf)

residuals = run!(model, config, output=VTK(), inner_loops=2)
