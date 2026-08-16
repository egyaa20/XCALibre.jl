# =============================================================================
# Shared helpers for all tatsumoto stage scripts.
# Include AFTER `using XCALibre, TOML`.
# =============================================================================

# cases/tatsumoto/ — configs/, stages/, pipeline/, runs/ live under here.
const CAMPAIGN_ROOT = dirname(@__DIR__)

"""
    load_case(args; default) -> (case_file, cfg, run_dir)

Parse the case TOML given as ARGS[1] (falling back to `configs/<default>`).
Every output of a run (checkpoints, VTK, post) lives under its `run_dir`:
`[case] run_dir` if present (set by the pipeline submitter), otherwise
`cases/tatsumoto/runs/<case.name>/`.
"""
function load_case(args; default)
    case_file = length(args) >= 1 ? abspath(args[1]) :
                joinpath(CAMPAIGN_ROOT, "configs", default)
    cfg = TOML.parsefile(case_file)
    name = cfg["case"]["name"]
    run_dir = haskey(cfg["case"], "run_dir") ? abspath(cfg["case"]["run_dir"]) :
              joinpath(CAMPAIGN_ROOT, "runs", name)
    mkpath(run_dir)
    @info "Case '$name'" case_file run_dir
    return case_file, cfg, run_dir
end

"""
    warmup_checkpoint(cfg, run_dir) -> path

Path of the warmup checkpoint this case reads/writes. Batch variants that share
the baseline warmup carry `[run] warmup_from = "<other run name>"` (set by the
submitter) and read that run's checkpoint instead of their own.
"""
function warmup_checkpoint(cfg, run_dir)
    wend = Float64(cfg["run"]["warmup_end"])
    src  = get(cfg["run"], "warmup_from", nothing)
    dir  = src === nothing ? run_dir : joinpath(CAMPAIGN_ROOT, "runs", src)
    return joinpath(dir, "warmup_t$(wend).jld2")
end

"""
    resolve_grid(gridname) -> path

Find `examples/0_GRIDS/<gridname>.unv`; if only the committed `.unv.gz` exists
(the quarter_pipe grids are stored gzipped — xfine is 130 MB raw), decompress
it once with system gzip (present on any Linux/HPC; on Windows run
`gzip -dk <file>` from Git Bash manually if this step fails).
"""
function resolve_grid(gridname)
    dir = pkgdir(XCALibre, "examples/0_GRIDS")
    unv = joinpath(dir, gridname * ".unv")
    isfile(unv) && return unv
    gz = unv * ".gz"
    if isfile(gz)
        @info "Decompressing $(basename(gz)) (one-time)..."
        success(`gzip -dk $gz`) ||
            error("gzip -dk failed for $gz — decompress manually and rerun.")
        return unv
    end
    error("Grid '$gridname' not found: neither $unv nor $gz exists.")
end

"""
    resolve_mesh(cfg) -> path

Mesh file for this case, covering both sources:

  * generated — `[mesh.geometry]` present: the mesh was built by SALOME on the
    cluster and `submit.jl` recorded its cache path in `[mesh] generated_file`.
    Nothing is hashed or derived here; the materialised case.toml is the single
    source of truth so a run always uses exactly the mesh its DAG built.
  * named — legacy `[mesh] grid = "quarter_pipe_fine"` from examples/0_GRIDS.
"""
function resolve_mesh(cfg)
    m = cfg["mesh"]
    if haskey(m, "generated_file")
        f = abspath(m["generated_file"])
        isfile(f) || error("""
            Generated mesh missing: $f
            The mesh stage should have built it. Rebuild with:
              julia --project=. cases/tatsumoto/pipeline/submit.jl <case.toml> --stages mesh
            """)
        return f
    end
    haskey(m, "grid") || error(
        "[mesh] needs either `grid = \"<name>\"` or a [mesh.geometry] block")
    return resolve_grid(m["grid"])
end

"""
    run_solver!(model, config; kwargs...)

Call the multiphase solver with only the kwargs the INSTALLED solver supports.
The solver is being rebuilt feature-by-feature (PIMPLE outer correctors,
stop_T checkpointing, ... come and go); unsupported kwargs are dropped with a
warning instead of a MethodError, so the same case TOML works throughout.
"""
function run_solver!(model, config; kwargs...)
    f = XCALibre.Solvers.multiphase!
    T = Tuple{typeof(model), typeof(config)}
    good = Pair{Symbol,Any}[]
    dropped = Symbol[]
    for (k, v) in pairs(kwargs)
        v === nothing && continue
        hasmethod(f, T, (k,)) ? push!(good, k => v) : push!(dropped, k)
    end
    isempty(dropped) ||
        @warn "Installed solver does not support these options (ignored): $(dropped)"
    return f(model, config; good...)
end
