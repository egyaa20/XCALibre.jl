# =============================================================================
# Stage 2 — FLOW DEVELOPMENT (warmup)
# -----------------------------------------------------------------------------
# Runs the Tatsumoto case with NO wall heat flux until `run.warmup_end`, then
# writes a JLD2 checkpoint of the developed flow field for stage 3 (heated.jl)
# to restart from.
#
# All parameters come from case.toml — nothing physics-related is hard-coded.
#
# Usage:
#     julia --project=. cases/tatsumoto/stages/prep.jl [path/to/case.toml]
#
# (defaults to configs/supercritical.toml; outputs land in runs/<case.name>/)
# =============================================================================

using XCALibre
using TOML
using LinearAlgebra
using StaticArrays

# --- Load configuration ------------------------------------------------------
include(joinpath(@__DIR__, "stage_common.jl"))
const CASE_FILE, CFG, RUN_DIR = load_case(ARGS; default="supercritical.toml")

# Convenience accessors (keep the rest of the file readable).
mesh_cfg   = CFG["mesh"];   hw_cfg   = CFG["hardware"]
flow_cfg   = CFG["flow"];   thermo   = CFG["thermo"]
energy_cfg = CFG["energy"]; run_cfg  = CFG["run"]
heat_cfg   = CFG["heating"]; adapt_cfg = run_cfg["adaptive"]
Q0 = Float64(heat_cfg["Q0"])   # constant wall flux — matches heated.jl at t=0, eliminates discontinuity

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

# Inlet turbulence from intensity + hydraulic diameter (standard estimates).
TI    = Float64(flow_cfg["turb_intensity"])
D_HYD = Float64(flow_cfg["D_hyd"])
ELL   = 0.07 * D_HYD                              # turbulence length scale
K_INLET     = 1.5 * (TI * U_IN)^2
OMEGA_INLET = sqrt(K_INLET) / (ELL * 0.09^0.25)
NUT_INLET   = K_INLET / OMEGA_INLET

# Flow is axial along +z, opposing gravity.
inletVelocity  = [0.0, 0.0, U_IN]
noSlipVelocity = [0.0, 0.0, 0.0]
gravity        = Gravity([0.0, 0.0, -9.81])

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

# --- Physics model -----------------------------------------------------------
model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        model = Mixture(diameter = 1.0e-7),
        # Phase constants are placeholders — HelmholtzTable overwrites them at
        # T_snapshot and the live update then drives per-cell properties.
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

# --- Boundary conditions (constant Q0 on Wall_Heated) ------------------------
# Wall_Heated carries Q0 W/m² (same as heated.jl at t=0) so the checkpoint
# already has the near-wall thermal gradient — no step discontinuity on restart.
energy_BCs = if ENERGY_FORM == "enthalpy"
    (
        T = [
            Dirichlet(:Inlet, T_IN),
            Neumann(:Wall_Heated,   0.0),     # T is derived; flux goes on h
            Neumann(:Wall_Unheated, 0.0),
            Symmetry(:Symmetry1),
            Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
        h = [
            Dirichlet(:Inlet, H_IN),
            Neumann(:Wall_Heated,   Q0),      # constant Q0 W/m²
            Neumann(:Wall_Unheated, 0.0),
            Symmetry(:Symmetry1),
            Symmetry(:Symmetry2),
            Zerogradient(:Outlet),
        ],
    )
else
    (
        T = [
            Dirichlet(:Inlet, T_IN),
            Neumann(:Wall_Heated,   Q0),      # constant Q0 W/m²
            Neumann(:Wall_Unheated, 0.0),
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
                Dirichlet(:Wall_Heated, 0.0),
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
                Dirichlet(:Wall_Heated, 0.0),
                Dirichlet(:Wall_Unheated, 0.0),
                Symmetry(:Symmetry1),
                Symmetry(:Symmetry2),
                Zerogradient(:Outlet),
            ],
        ),
        energy_BCs,
    )
)

# --- Numerics ----------------------------------------------------------------
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
        convergence = 1e-7, relax = 0.3, rtol = 0.0, atol = 1.0e-6, itmax = 25000,
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
    write_interval = Int(run_cfg["write_interval_prep"]),
    adaptive       = adaptive,
    t_end          = Float64(run_cfg["warmup_end"]),   # exact stop at warmup_end
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

@info "Starting warmup -> t_end = $(run_cfg["warmup_end"]) s ..."
@time residuals = run_solver!(model, config;
    inner_loops        = Int(run_cfg["inner_loops"]),
    n_outer_correctors = Int(run_cfg["n_outer_correctors"]),
    outer_tol          = Float64(run_cfg["outer_tol"]),
)

# --- Save checkpoint for stage 3 ---------------------------------------------
WARMUP_END = Float64(run_cfg["warmup_end"])
checkpoint = joinpath(RUN_DIR, "warmup_t$(WARMUP_END).jld2")
save_checkpoint(checkpoint; time = WARMUP_END,
    U     = model.momentum.U,
    p_rgh = model.fluid.p_rgh,
    alpha = model.fluid.alpha,
    T     = model.energy.T,
    k     = model.turbulence.k,
    omega = model.turbulence.omega,
    nut   = model.turbulence.nut,
)
@info "Warmup complete — checkpoint saved -> $(checkpoint)"
