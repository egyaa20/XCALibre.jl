# =============================================================================
# Two-phase boiling case — STAGE 3: HEATED RUN (subcooled flow boiling)
# -----------------------------------------------------------------------------
# Restarts the developed liquid flow and applies a ramped wall heat flux
# q(t) = Q0 e^{t/tau} on Wall_Heated. Boiling is modelled by:
#   * RPI (Kurul-Podowski) wall heat-flux partitioning at the heated wall
#     (q_w = q_c convection + q_q quenching + q_e evaporation) -> captures the
#     boiling heat-transfer enhancement that sets T_w for a given q.
#   * Lee bulk phase change in the subcooled core (latent heat sink/source).
#
# Energy = MultiphaseTemperature (the boiling source dispatcher only fires for
# this model). Subcritical 0.5 MPa -> properties smooth, T-formulation stable.
#
# CAVEAT (first-pass): the phase-change alpha (void-fraction) source is not yet
# wired into the alpha equation, so the VOID FRACTION won't grow — but the WALL
# HEAT TRANSFER (what the q-ΔT curve measures) is captured via the RPI partition
# and Lee latent heat. Validate the q-ΔT_L curve up to ~DNB.
#
# Usage:  julia --project=. cases/tatsumoto/stages/boiling_heated.jl [path/to/boiling.toml]
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
sat_cfg  = CFG["saturation"]; boil_cfg = CFG["boiling"]
heat_cfg = CFG["heating"]; run_cfg  = CFG["run"]; adapt_cfg = run_cfg["adaptive"]

CASE_NAME   = CFG["case"]["name"]
WARMUP_END  = Float64(run_cfg["warmup_end"])
CHECKPOINT  = warmup_checkpoint(CFG, RUN_DIR)
RESULTS_DIR = joinpath(RUN_DIR, "vtk")
isfile(CHECKPOINT) || error("Checkpoint not found: $(CHECKPOINT)\nRun boiling_prep.jl first.")

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
mesh_file = resolve_grid(mesh_cfg["grid"])
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

inletVelocity  = [0.0, 0.0, U_IN]
noSlipVelocity = [0.0, 0.0, 0.0]
gravity        = Gravity([0.0, 0.0, -9.81])

# --- Ramped wall heat flux  q(t) = Q0 e^{t/tau} ------------------------------
const Q0  = Float64(heat_cfg["Q0"])
const TAU = Float64(heat_cfg["tau"])
heated_flux(coords, time, fID) = Q0 * exp(time / TAU)

# --- Saturation properties (corrected N2 @ 0.5 MPa; paper T_sat = 94.02 K) ----
sat = SatProps(
    pressure = Float64(CFG["thermo"]["pressure"]),
    T_sat = Float64(sat_cfg["T_sat"]), h_fg = Float64(sat_cfg["h_fg"]),
    rho_l = Float64(sat_cfg["rho_l"]), rho_v = Float64(sat_cfg["rho_v"]),
    cp_l  = Float64(sat_cfg["cp_l"]),  mu_l  = Float64(sat_cfg["mu_l"]),
    k_l   = Float64(sat_cfg["k_l"]),   sigma = Float64(sat_cfg["sigma"]),
)

# --- Boiling models (built on the device mesh) -------------------------------
H_C   = Float64(boil_cfg["h_c"])
wall_boiling = init_rpi_wall_state(RPI(), sat, heated_flux, H_C, mesh_dev, [:Wall_Heated])
phase_change = Lee(C_e = Float64(boil_cfg["C_e"]), C_c = Float64(boil_cfg["C_c"]))

# --- Physics: constant-property two-phase + boiling --------------------------
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
        gravity      = gravity,
        sat_props    = sat,
        phase_change = phase_change,
        wall_boiling = wall_boiling,
    ),
    turbulence = RANS{KOmegaSST}(walls = (:Wall_Heated, :Wall_Unheated)),
    energy = Energy{MultiphaseTemperature}(T_init = T_IN, Pr_t = PR_T),
    domain = mesh_dev,
)

# --- Boundary conditions (HEATED: full ramped flux on Wall_Heated) -----------
# The NeumannFunction applies the FULL q_w; RPI subtracts the evaporation
# portion q_e via model.energy.T_source each corrector.
BCs = assign(
    region = mesh_dev,
    (
        U = [
            Dirichlet(:Inlet, inletVelocity),
            Wall(:Wall_Heated, noSlipVelocity), Wall(:Wall_Unheated, noSlipVelocity),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        p_rgh = [
            Wall(:Inlet), Wall(:Wall_Heated), Wall(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Dirichlet(:Outlet, 0.0),
        ],
        alpha = [
            Dirichlet(:Inlet, 1.0),
            Zerogradient(:Wall_Heated), Zerogradient(:Wall_Unheated),
            Symmetry(:Symmetry1), Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        T = [
            Dirichlet(:Inlet, T_IN),
            NeumannFunction(:Wall_Heated, heated_flux),   # ramped Joule flux
            Neumann(:Wall_Unheated, 0.0),                 # adiabatic
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
        convergence = 1e-7, relax = 0.3, rtol = 0.0, atol = 1.0e-6, itmax = 20000),
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
    write_interval = Int(run_cfg["write_interval_heated"]),
    adaptive       = adaptive,
    t_end          = Float64(run_cfg["heated_end"]),
)

config = Configuration(
    solvers = solvers, schemes = schemes, runtime = runtime,
    hardware = hardware, boundaries = BCs,
)

# --- Restore developed flow --------------------------------------------------
GC.gc()
ckpt_time = load_checkpoint!(CHECKPOINT;
    U = model.momentum.U, p_rgh = model.fluid.p_rgh, alpha = model.fluid.alpha,
    T = model.energy.T, k = model.turbulence.k, omega = model.turbulence.omega,
    nut = model.turbulence.nut,
)
@info "Restored checkpoint (warmup t = $(ckpt_time) s) -> heating from t = 0"

# --- Run in a clean results dir ----------------------------------------------
rm(RESULTS_DIR; recursive = true, force = true)
mkpath(RESULTS_DIR)

@info "Starting boiling heated run -> t_end = $(run_cfg["heated_end"]) s  (Q0=$(Q0), tau=$(TAU), T_sat=$(sat.T_sat))"
residuals = cd(RESULTS_DIR) do
    @time run_solver!(model, config;
        inner_loops        = Int(run_cfg["inner_loops"]),
        n_outer_correctors = Int(run_cfg["n_outer_correctors"]),
        outer_tol          = Float64(run_cfg["outer_tol"]),
    )
end

@info "Boiling heated run complete — results in $(RESULTS_DIR)"
