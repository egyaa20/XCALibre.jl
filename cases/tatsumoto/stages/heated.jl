# =============================================================================
# Stage 3 — HEATED RUN
# -----------------------------------------------------------------------------
# Restarts from the stage-2 checkpoint (developed flow) and applies a ramped
# wall heat flux  q(t) = Q0 * exp(t / tau)  on `Wall_Heated` until
# `run.heated_end`. VTU snapshots and the `times.csv` (iteration -> physical
# time) sidecar are written into  <case>_results/  for stage-4 postprocessing.
#
# Heating starts at the heated-run clock t = 0 (no offset). All parameters
# come from case.toml.
#
# Usage:
#     julia --project=. cases/tatsumoto/stages/heated.jl [path/to/case.toml] [resume]
# =============================================================================

using XCALibre
using TOML
using LinearAlgebra
using StaticArrays

# --- Load configuration ------------------------------------------------------
include(joinpath(@__DIR__, "stage_common.jl"))
const CASE_FILE, CFG, RUN_DIR = load_case(ARGS; default="supercritical.toml")

mesh_cfg   = CFG["mesh"];   hw_cfg   = CFG["hardware"]
flow_cfg   = CFG["flow"];   thermo   = CFG["thermo"]
energy_cfg = CFG["energy"]; run_cfg  = CFG["run"]
heat_cfg   = CFG["heating"]; adapt_cfg = run_cfg["adaptive"]

# --- Resume mode -------------------------------------------------------------
# `julia heated.jl case.toml resume` restarts from the pseudo-critical (stop_T)
# checkpoint instead of the warmup, and continues the heating ramp seamlessly
# (the heating clock is offset by the checkpoint time). Results land in a
# separate `<case>_resume_results/` dir so the original run is preserved.
const RESUME = length(ARGS) >= 2 && lowercase(ARGS[2]) == "resume"

# --- Paths (absolute, so the run can cd into the results dir) ----------------
CASE_NAME   = CFG["case"]["name"]
WARMUP_END  = Float64(run_cfg["warmup_end"])
STOP_T_K    = Float64(get(run_cfg, "stop_T", 124.0))
WARMUP_CKPT = warmup_checkpoint(CFG, RUN_DIR)
STOP_CKPT   = joinpath(RUN_DIR, "heated_stop$(STOP_T_K)K.jld2")
CHECKPOINT  = RESUME ? STOP_CKPT : WARMUP_CKPT
RESULTS_DIR = joinpath(RUN_DIR, RESUME ? "vtk_resume" : "vtk")
isfile(CHECKPOINT) || error("Checkpoint not found: $(CHECKPOINT)\n" *
    (RESUME ? "Run heated.jl (without 'resume') first to create the stop_T checkpoint." :
              "Run prep.jl first."))

# --- Hardware backend --------------------------------------------------------
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

# --- Flow / turbulence inlet quantities --------------------------------------
U_IN  = Float64(flow_cfg["U_in"])
T_IN  = Float64(flow_cfg["T_in"])
P_OP  = Float64(thermo["pressure"])
PR_T  = Float64(energy_cfg["Pr_t"])
ENERGY_FORM = lowercase(get(energy_cfg, "formulation", "temperature"))
ENERGY_FORM in ("temperature", "enthalpy") ||
    error("energy.formulation must be \"temperature\" or \"enthalpy\", got $(ENERGY_FORM)")

TI    = Float64(flow_cfg["turb_intensity"])
D_HYD = Float64(flow_cfg["D_hyd"])
ELL   = 0.07 * D_HYD
K_INLET     = 1.5 * (TI * U_IN)^2
OMEGA_INLET = sqrt(K_INLET) / (ELL * 0.09^0.25)

inletVelocity  = [0.0, 0.0, U_IN]
noSlipVelocity = [0.0, 0.0, 0.0]
gravity        = Gravity([0.0, 0.0, -9.81])

# --- Heated wall flux: q(t) = Q0 * exp((t0 + t) / tau) -----------------------
# `HEAT_T0` offsets the heating clock so a resumed run continues the ramp from
# the checkpoint time rather than restarting it at 0. Stays 0 for a fresh run.
const Q0  = Float64(heat_cfg["Q0"])
const TAU = Float64(heat_cfg["tau"])
const HEAT_T0 = Ref(0.0)
heated_flux(coords, time, fID) = Q0 * exp((HEAT_T0[] + time) / TAU)

# --- Real-fluid table (one handle: drives properties AND seeds enthalpy) ------
htable = HelmholtzTable(
    N2();
    p_operating = P_OP,
    T_range     = (Float64(thermo["T_table_lo"]), Float64(thermo["T_table_hi"])),
    p_range     = (Float64(thermo["p_lo"]), Float64(thermo["p_hi"])),
    n_T         = Int(thermo["n_points"]),
    n_p         = Int(thermo["n_p"]),
    T_snapshot  = T_IN,
)
H_IN = Ttoh_helmholtz(htable, T_IN)   # inlet enthalpy for the h-formulation

# --- Energy model (temperature or enthalpy formulation, per case.toml) --------
energy_model = if ENERGY_FORM == "enthalpy"
    Energy{VariableSensibleEnthalpy}(
        T_init = T_IN, Pr_t = PR_T,
        Pr_t_model = Symbol(get(energy_cfg, "Pr_t_model", "constant")),
        prop_relax = Float64(get(energy_cfg, "prop_relax", 1.0)),
    )
else
    Energy{MultiphaseTemperature}(
        T_init = T_IN, Pr_t = PR_T,
        Pr_t_model = Symbol(get(energy_cfg, "Pr_t_model", "constant")),
        prop_relax = Float64(get(energy_cfg, "prop_relax", 1.0)),
    )
end

# --- Physics model (same as prep.jl) -----------------------------------------
model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        model = Mixture(diameter = 1.0e-7),
        phases = (
            Phase(rho = 712.0, mu = 76.0e-6, cp = 2200.0, k = 0.105),
            Phase(rho = 712.0, mu = 76.0e-6, cp = 2200.0, k = 0.105),
        ),
        gravity = gravity,
        fluid_properties = htable,
    ),
    turbulence = RANS{KOmegaSST}(walls = (:Wall_Heated, :Wall_Unheated)),
    energy = energy_model,
    domain = mesh_dev,
)

# --- Boundary conditions (HEATED: ramped flux on Wall_Heated) ----------------
# In the enthalpy formulation the ramped flux drives `h` (the solved variable)
# and `T` is derived (face BCs only, no flux). In the temperature formulation
# the flux drives `T` directly. The wall-flux function is identical either way
# (NeumannFunction asserts the physical heat flux regardless of coefficient).
energy_BCs = if ENERGY_FORM == "enthalpy"
    (
        T = [
            Dirichlet(:Inlet, T_IN),
            Neumann(:Wall_Heated,   0.0),
            Neumann(:Wall_Unheated, 0.0),
            Symmetry(:Symmetry1),
            Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        h = [
            Dirichlet(:Inlet, H_IN),
            NeumannFunction(:Wall_Heated, heated_flux),   # ramped Joule flux on h
            Neumann(:Wall_Unheated, 0.0),                 # adiabatic
            Symmetry(:Symmetry1),
            Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
    )
else
    (
        T = [
            Dirichlet(:Inlet, T_IN),
            NeumannFunction(:Wall_Heated, heated_flux),   # ramped Joule flux on T
            Neumann(:Wall_Unheated, 0.0),                 # adiabatic
            Symmetry(:Symmetry1),
            Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
    )
end

BCs = assign(
    region = mesh_dev,
    merge(
        (
            U = [
                Dirichlet(:Inlet, inletVelocity),
                Wall(:Wall_Heated,   noSlipVelocity),
                Wall(:Wall_Unheated, noSlipVelocity),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            p_rgh = [
                Wall(:Inlet),
                Wall(:Wall_Heated),
                Wall(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Dirichlet(:Outlet, 0.0),
            ],
            alpha = [
                Dirichlet(:Inlet, 1.0),
                Zerogradient(:Wall_Heated),
                Zerogradient(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            k = [
                Dirichlet(:Inlet, K_INLET),
                Dirichlet(:Wall_Heated, 0.0), #Dirichlet
                Dirichlet(:Wall_Unheated, 0.0),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            omega = [
                Dirichlet(:Inlet, OMEGA_INLET),
                OmegaWallFunction(:Wall_Heated),
                OmegaWallFunction(:Wall_Unheated),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
            nut = [
                Extrapolated(:Inlet),
                Dirichlet(:Wall_Heated, 0.0), #Dirichlet
                Dirichlet(:Wall_Unheated, 0.0),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
        ),
        energy_BCs,
    )
)

# --- Numerics (same as prep.jl) ----------------------------------------------
schemes = (
    U     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
    p     = Schemes(time = Euler, gradient   = Gauss,  laplacian = Linear),
    p_rgh = Schemes(time = Euler, gradient   = Gauss,  laplacian = Linear),
    alpha = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
    T     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),
    h     = Schemes(time = Euler, divergence = Upwind, laplacian = Linear),  # enthalpy formulation
    k     = Schemes(divergence = Upwind),
    omega = Schemes(divergence = Upwind),
    y     = Schemes(),
)

solvers = (
    U = SolverSetup(
        solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.7, rtol = 0.0, atol = 1.0e-5,
    ),
    p_rgh = SolverSetup(
        solver = Cg(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.3, rtol = 0.0, atol = 1.0e-6, itmax = 20000,
    ),
    alpha = SolverSetup(
        solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 1.0, rtol = 0.0, atol = 1.0e-5,
    ),
    T = SolverSetup(
        solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.7, rtol = 0.0, atol = 1.0e-7,
    ),
    h = SolverSetup(    # enthalpy ~1e5 J/kg, so atol is scaled up from T's
        solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-7, relax = 0.7, rtol = 1.0e-6, atol = 1.0,
    ),
    k = SolverSetup(
        solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-10, relax = 0.6, rtol = 1e-3,
    ),
    omega = SolverSetup(
        solver = Bicgstab(), preconditioner = Jacobi(),
        convergence = 1e-10, relax = 0.6, rtol = 1e-3,
    ),
    y = SolverSetup(
        solver = Cg(), preconditioner = Jacobi(),
        convergence = 1e-10, relax = 0.7, rtol = 1e-5, itmax = 5000,
    ),
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

# --- Restore developed flow from checkpoint ----------------------------------
GC.gc()
ckpt_time = load_checkpoint!(CHECKPOINT;
    U     = model.momentum.U,
    p_rgh = model.fluid.p_rgh,
    alpha = model.fluid.alpha,
    T     = model.energy.T,
    k     = model.turbulence.k,
    omega = model.turbulence.omega,
    nut   = model.turbulence.nut,
)
if RESUME
    HEAT_T0[] = ckpt_time      # continue the q-ramp from the checkpoint time
    @info "RESUME: restored stop checkpoint -> heating ramp continues" ckpt_time q_now=Q0*exp(ckpt_time/TAU)
else
    @info "Restored checkpoint (warmup t = $(ckpt_time) s) -> heating from t = 0"
end

# --- Run inside a clean results dir (VTU + times.csv land there) -------------
rm(RESULTS_DIR; recursive = true, force = true)
mkpath(RESULTS_DIR)

# --- Pseudo-critical checkpoint (one-shot, run continues) ----------------------
# The FIRST time the hottest cell reaches `stop_T` (just past the pseudo-critical
# line) dump a fully restartable checkpoint — every solved field + the physical
# time — then carry on to heated_end. Restart with `heated.jl <case> resume`.
# (STOP_T_K / STOP_CKPT are defined in the paths block.) Disabled on a resume run
# so it doesn't overwrite the original checkpoint.
function save_stop_checkpoint(m, t)
    base = (
        U     = m.momentum.U,
        p_rgh = m.fluid.p_rgh,
        alpha = m.fluid.alpha,
        T     = m.energy.T,
        k     = m.turbulence.k,
        omega = m.turbulence.omega,
        nut   = m.turbulence.nut,
    )
    # In the enthalpy formulation `h` is the primary solved field — persist it too.
    fields = ENERGY_FORM == "enthalpy" ? merge(base, (h = m.energy.h,)) : base
    save_checkpoint(STOP_CKPT; time = t, fields...)
    @info "Stop checkpoint saved" file=STOP_CKPT time=t Q0=Q0 tau=TAU q_wall=Q0*exp(t/TAU)
end

@info "Starting heated run -> t_end = $(run_cfg["heated_end"]) s  (Q0=$(Q0), tau=$(TAU), stop_T=$(STOP_T_K) K)"
residuals = cd(RESULTS_DIR) do
    @time run_solver!(model, config;
        inner_loops        = Int(run_cfg["inner_loops"]),
        n_outer_correctors = Int(run_cfg["n_outer_correctors"]),
        outer_tol          = Float64(run_cfg["outer_tol"]),
        stop_T             = RESUME ? nothing : STOP_T_K,
        stop_save          = save_stop_checkpoint,
    )
end

@info "Heated run complete — results in $(RESULTS_DIR)"
