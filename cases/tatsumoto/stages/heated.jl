# =============================================================================
# Stage 3 — HEATED RUN
# -----------------------------------------------------------------------------
# Restarts from the stage-2 checkpoint (developed flow) and applies the ramped
# wall heat flux  q(t) = Q0 * exp(t / tau)  on Wall_Heated via HeatFluxFunction
# until `run.heated_end`. VTU snapshots + times.csv land in <run_dir>/vtk/ for
# stage-4 postprocessing; a final JLD2 checkpoint is saved on completion.
#
# Usage:
#     julia --project=. cases/tatsumoto/stages/heated.jl [path/to/case.toml] [resume]
#
# `resume` restarts from the final checkpoint of a previous heated run and
# continues the ramp from its saved time (results land in vtk_resume/).
# =============================================================================

using XCALibre
using TOML
using JLD2

include(joinpath(@__DIR__, "stage_common.jl"))
include(joinpath(@__DIR__, "case_model.jl"))

const CASE_FILE, CFG, RUN_DIR = load_case(ARGS; default="supercritical.toml")
const RESUME = length(ARGS) >= 2 && lowercase(ARGS[2]) == "resume"

const Q0  = Float64(CFG["heating"]["Q0"])
const TAU = Float64(CFG["heating"]["tau"])
HEATED_END = Float64(CFG["run"]["heated_end"])

FINAL_CKPT  = joinpath(RUN_DIR, "heated_final.jld2")
WARMUP_CKPT = warmup_checkpoint(CFG, RUN_DIR)
CHECKPOINT  = RESUME ? FINAL_CKPT : WARMUP_CKPT
RESULTS_DIR = joinpath(RUN_DIR, RESUME ? "vtk_resume" : "vtk")
isfile(CHECKPOINT) || error("Checkpoint not found: $(CHECKPOINT)\n" *
    (RESUME ? "Run heated.jl (without 'resume') first." : "Run prep.jl first."))

# Heating clock offset: 0 for a fresh run; the checkpoint time on resume so
# the exponential ramp continues instead of restarting.
const HEAT_T0 = Ref(0.0)
heated_flux(coords, time, i) = Q0 * exp((HEAT_T0[] + time) / TAU)

backend, hardware = build_backend(CFG["hardware"])
mesh_file = resolve_mesh(CFG)
mesh      = UNV3D_mesh(mesh_file, scale = Float64(CFG["mesh"]["scale"]))
mesh_dev  = adapt(backend, mesh)

model, config, htable, inits = build_case(CFG, mesh_dev, hardware;
    wall_heat_bc   = HeatFluxFunction(:Wall_Heated, heated_flux),
    write_interval = Int(CFG["run"]["write_interval_heated"]),
    t_end          = HEATED_END,
)

GC.gc()
initialise_fields!(model, inits)   # sane values under any field the checkpoint lacks
ckpt_time = load_checkpoint!(CHECKPOINT; checkpoint_fields(model)...)
if RESUME
    HEAT_T0[] = ckpt_time
    @info "RESUME: ramp continues" ckpt_time q_now=Q0*exp(ckpt_time/TAU)
else
    @info "Restored warmup checkpoint (t = $(ckpt_time) s) -> heating from t = 0"
end

rm(RESULTS_DIR; recursive = true, force = true)
mkpath(RESULTS_DIR)

@info "Starting heated run -> t_end = $(HEATED_END) s  (Q0=$(Q0), tau=$(TAU))"
residuals = cd(RESULTS_DIR) do
    @time run_solver!(model, config;
        inner_loops        = Int(CFG["run"]["inner_loops"]),
        n_outer_correctors = Int(get(CFG["run"], "n_outer_correctors", 3)),
        outer_tol          = Float64(get(CFG["run"], "outer_tol", 1.0e-4)),
        prop_relax         = Float64(get(CFG["run"], "prop_relax", 0.5)),
    )
end

t_final = RESUME ? HEAT_T0[] + HEATED_END : HEATED_END
save_checkpoint(FINAL_CKPT; time = t_final, checkpoint_fields(model)...)
@info "Heated run complete" results=RESULTS_DIR checkpoint=FINAL_CKPT t_final
