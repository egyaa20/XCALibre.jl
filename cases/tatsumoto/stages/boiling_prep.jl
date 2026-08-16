# =============================================================================
# Two-phase boiling case — STAGE 2: FLOW DEVELOPMENT (warmup)
# -----------------------------------------------------------------------------
# Develops the single-phase liquid-N2 turbulent flow at 0.5 MPa with NO wall
# heat flux, then checkpoints the developed field for boiling_heated.jl.
#
# Constant per-phase properties (subcritical -> no HelmholtzTable needed).
# Energy = MultiphaseTemperature (the formulation the boiling sources require).
#
# Usage:  julia --project=. cases/tatsumoto/stages/boiling_prep.jl [path/to/boiling.toml]
# =============================================================================

using XCALibre
using TOML
using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "stage_common.jl"))
const CASE_FILE, CFG, RUN_DIR = load_case(ARGS; default="boiling.toml")

mesh_cfg = CFG["mesh"];   hw_cfg   = CFG["hardware"]
flow_cfg = CFG["flow"];   energy_cfg = CFG["energy"]
liq_cfg  = CFG["phase"]["liquid"]; vap_cfg = CFG["phase"]["vapour"]
run_cfg  = CFG["run"];    adapt_cfg = run_cfg["adaptive"]

# --- Hardware ----------------------------------------------------------------
if uppercase(hw_cfg["backend"]) == "CUDA"
    using CUDA
    backend   = CUDABackend()
    workgroup = Int(hw_cfg["workgroup"])
else
    backend   = CPU()
    workgroup = AutoTune()
    activate_multithread(backend)
end
hardware = Hardware(backend = backend, workgroup = workgroup)

# --- Mesh --------------------------------------------------------------------
mesh_file = resolve_mesh(CFG)
mesh      = UNV3D_mesh(mesh_file, scale = Float64(mesh_cfg["scale"]))
mesh_dev  = adapt(backend, mesh)

# --- Inlet quantities --------------------------------------------------------
U_IN  = Float64(flow_cfg["U_in"])
T_IN  = Float64(flow_cfg["T_in"])
PR_T  = Float64(energy_cfg["Pr_t"])

TI    = Float64(flow_cfg["turb_intensity"])
D_HYD = Float64(flow_cfg["D_hyd"])
ELL   = 0.07 * D_HYD
K_INLET     = 1.5 * (TI * U_IN)^2
OMEGA_INLET = sqrt(K_INLET) / (ELL * 0.09^0.25)
NUT_INLET   = K_INLET / OMEGA_INLET

inletVelocity  = [0.0, 0.0, U_IN]
noSlipVelocity = [0.0, 0.0, 0.0]
gravity        = Gravity([0.0, 0.0, -9.81])

# --- Physics: constant-property two-phase mixture (no heat in warmup) --------
model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        model  = Mixture(diameter = 1.0e-7),
        phases = (
            Phase(rho = Float64(liq_cfg["rho"]), mu = Float64(liq_cfg["mu"]),
                  cp  = Float64(liq_cfg["cp"]),  k  = Float64(liq_cfg["k"])),   # 1: liquid
            Phase(rho = Float64(vap_cfg["rho"]), mu = Float64(vap_cfg["mu"]),
                  cp  = Float64(vap_cfg["cp"]),  k  = Float64(vap_cfg["k"])),   # 2: vapour
        ),
        gravity = gravity,
    ),
    turbulence = RANS{KOmegaSST}(walls = (:Wall_Heated, :Wall_Unheated)),
    energy = Energy{MultiphaseTemperature}(T_init = T_IN, Pr_t = PR_T),
    domain = mesh_dev,
)

# --- Boundary conditions (UNHEATED, single-phase liquid) ---------------------
BCs = assign(
    region = mesh_dev,
    (
        U = [
            Dirichlet(:Inlet, inletVelocity),
            Wall(:Wall_Heated,   noSlipVelocity),
            Wall(:Wall_Unheated, noSlipVelocity),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        p_rgh = [
            Wall(:Inlet), Wall(:Wall_Heated), Wall(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Dirichlet(:Outlet, 0.0),
        ],
        alpha = [
            Dirichlet(:Inlet, 1.0),                 # pure liquid in
            Zerogradient(:Wall_Heated), Zerogradient(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        T = [
            Dirichlet(:Inlet, T_IN),
            Neumann(:Wall_Heated, 0.0), Neumann(:Wall_Unheated, 0.0),  # adiabatic
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        k = [
            Dirichlet(:Inlet, K_INLET),
            KWallFunction(:Wall_Heated), KWallFunction(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        omega = [
            Dirichlet(:Inlet, OMEGA_INLET),
            OmegaWallFunction(:Wall_Heated), OmegaWallFunction(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        nut = [
            Extrapolated(:Inlet),
            NutWallFunction(:Wall_Heated), NutWallFunction(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
    )
)

# --- Numerics ----------------------------------------------------------------
schemes = (
    U     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
    p     = Schemes(time = Euler, gradient   = Gauss,  laplacian = Linear),
    p_rgh = Schemes(time = Euler, gradient   = Gauss,  laplacian = Linear),
    alpha = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
    T     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
    k     = Schemes(divergence = Upwind),
    omega = Schemes(divergence = Upwind),
    y     = Schemes(),
)

solvers = (
    U = SolverSetup(solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.7, rtol = 0.0, atol = 1.0e-5),
    p_rgh = SolverSetup(solver = Cg(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.3, rtol = 0.0, atol = 1.0e-6, itmax = 25000),
    alpha = SolverSetup(solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 1.0, rtol = 0.0, atol = 1.0e-5),
    T = SolverSetup(solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.7, rtol = 0.0, atol = 1.0e-7),
    k = SolverSetup(solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-10, relax = 0.6, rtol = 1e-3),
    omega = SolverSetup(solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-10, relax = 0.6, rtol = 1e-3),
    y = SolverSetup(solver = Cg(), preconditioner = Jacobi(),
        convergence = 1e-10, relax = 0.7, rtol = 1e-5, itmax = 5000),
)

adaptive = AdaptiveTimeStepping(
    maxCo      = Float64(adapt_cfg["maxCo"]),
    maxAlphaCo = Float64(adapt_cfg["maxAlphaCo"]),
    minShrink  = Float64(adapt_cfg["minShrink"]),
    maxGrow    = Float64(adapt_cfg["maxGrow"]),
)

runtime = Runtime(
    time_step      = Float64(run_cfg["time_step"]),
    write_interval = Int(run_cfg["write_interval_prep"]),
    adaptive       = adaptive,
    t_end          = Float64(run_cfg["warmup_end"]),
)

config = Configuration(
    solvers = solvers, schemes = schemes, runtime = runtime,
    hardware = hardware, boundaries = BCs,
)

# --- Initialise & run --------------------------------------------------------
GC.gc()
initialise!(model.fluid.p_rgh,      0.0)
initialise!(model.momentum.U,       inletVelocity)
initialise!(model.fluid.alpha,      1.0)
initialise!(model.energy.T,         T_IN)
initialise!(model.turbulence.k,     K_INLET)
initialise!(model.turbulence.omega, OMEGA_INLET)
initialise!(model.turbulence.nut,   NUT_INLET)

@info "Starting boiling-case warmup -> t_end = $(run_cfg["warmup_end"]) s ..."
@time residuals = run_solver!(model, config;
    inner_loops        = Int(run_cfg["inner_loops"]),
    n_outer_correctors = Int(run_cfg["n_outer_correctors"]),
    outer_tol          = Float64(run_cfg["outer_tol"]),
)

# --- Checkpoint --------------------------------------------------------------
WARMUP_END = Float64(run_cfg["warmup_end"])
checkpoint = joinpath(RUN_DIR, "warmup_t$(WARMUP_END).jld2")
save_checkpoint(checkpoint; time = WARMUP_END,
    U = model.momentum.U, p_rgh = model.fluid.p_rgh, alpha = model.fluid.alpha,
    T = model.energy.T, k = model.turbulence.k, omega = model.turbulence.omega,
    nut = model.turbulence.nut,
)
@info "Boiling-case warmup complete -> $(checkpoint)"
