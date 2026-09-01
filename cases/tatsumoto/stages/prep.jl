# =============================================================================
# Stage 2 — FLOW DEVELOPMENT (warmup)
# -----------------------------------------------------------------------------
# Runs the case with the CONSTANT baseline wall flux Q0 (matches heated.jl at
# t=0, no restart discontinuity) until `run.warmup_end`, then writes a JLD2
# checkpoint of the developed flow for stage 3 to restart from.
#
# Usage:
#     julia --project=. cases/tatsumoto/stages/prep.jl [path/to/case.toml]
# =============================================================================

using XCALibre
using TOML
using JLD2

include(joinpath(@__DIR__, "stage_common.jl"))
include(joinpath(@__DIR__, "case_model.jl"))

const CASE_FILE, CFG, RUN_DIR = load_case(ARGS; default="supercritical.toml")

Q0 = Float64(CFG["heating"]["Q0"])
WARMUP_END = Float64(CFG["run"]["warmup_end"])

backend, hardware = build_backend(CFG["hardware"])
mesh_file = resolve_mesh(CFG)
mesh      = UNV3D_mesh(mesh_file, scale = Float64(CFG["mesh"]["scale"]))
mesh_dev  = adapt(backend, mesh)

model, config, htable, inits = build_case(CFG, mesh_dev, hardware;
    wall_heat_bc   = HeatFlux(:Wall_Heated, Q0),
    write_interval = Int(CFG["run"]["write_interval_prep"]),
    t_end          = WARMUP_END,
)

GC.gc()
initialise_fields!(model, inits)

PREP_RESULTS = joinpath(RUN_DIR, "vtk_prep")
rm(PREP_RESULTS; recursive = true, force = true)
mkpath(PREP_RESULTS)

@info "Starting warmup -> t_end = $(WARMUP_END) s ..."
residuals = cd(PREP_RESULTS) do
    @time run_solver!(model, config;
        inner_loops        = Int(CFG["run"]["inner_loops"]),
        n_outer_correctors = Int(get(CFG["run"], "n_outer_correctors", 3)),
        outer_tol          = Float64(get(CFG["run"], "outer_tol", 1.0e-4)),
        prop_relax         = Float64(get(CFG["run"], "prop_relax", 0.5)),
    )
end

checkpoint = joinpath(RUN_DIR, "warmup_t$(WARMUP_END).jld2")
save_checkpoint(checkpoint; time = WARMUP_END, checkpoint_fields(model)...)
@info "Warmup complete — checkpoint saved -> $(checkpoint)"
