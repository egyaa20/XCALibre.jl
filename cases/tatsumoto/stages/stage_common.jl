# =============================================================================
# Shared helpers for all tatsumoto stage scripts.
# Include AFTER `using XCALibre, TOML, JLD2`.
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
    ckpt = get(cfg["run"], "warmup_checkpoint", nothing)
    ckpt !== nothing && return abspath(ckpt)
    wend = Float64(cfg["run"]["warmup_end"])
    src  = get(cfg["run"], "warmup_from", nothing)
    dir  = src === nothing ? run_dir : joinpath(CAMPAIGN_ROOT, "runs", src)
    return joinpath(dir, "warmup_t$(wend).jld2")
end

"""
    resolve_mesh(cfg) -> path

Mesh file for this case:

  * generated — `[mesh.geometry]` present: built by SALOME on the cluster,
    `submit.jl` recorded its cache path in `[mesh] generated_file`.
  * named — `[mesh] grid = "<name>"` from examples/0_GRIDS.
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
    unv = joinpath(pkgdir(XCALibre, "examples/0_GRIDS"), m["grid"] * ".unv")
    isfile(unv) && return unv
    gz = unv * ".gz"
    if isfile(gz)
        @info "Decompressing $(basename(gz)) (one-time)..."
        success(`gzip -dk $gz`) || error("gzip -dk failed for $gz")
        return unv
    end
    error("Grid '$(m["grid"])' not found: neither $unv nor $gz exists.")
end

"""
    save_checkpoint(path; time, fields...)

JLD2 checkpoint of the named fields' cell values plus the physical time.
"""
function save_checkpoint(path; time, fields...)
    data = Dict{String,Any}("time" => Float64(time))
    for (name, field) in pairs(fields)
        if field isa VectorField
            data[string(name)] = hcat(Array(field.x.values), Array(field.y.values),
                                      Array(field.z.values))
        else
            data[string(name)] = Array(field.values)
        end
    end
    JLD2.save(path, data)
    return path
end

"""
    load_checkpoint!(path; fields...) -> time

Restore the named fields from a JLD2 checkpoint written by `save_checkpoint`.
Fields absent from the file are skipped with a warning.
"""
function load_checkpoint!(path; fields...)
    data = JLD2.load(path)
    for (name, field) in pairs(fields)
        key = string(name)
        if !haskey(data, key)
            @warn "Checkpoint has no field '$key' — left at its initial value"
            continue
        end
        vals = data[key]
        if field isa VectorField
            copyto!(field.x.values, vals[:, 1])
            copyto!(field.y.values, vals[:, 2])
            copyto!(field.z.values, vals[:, 3])
        else
            copyto!(field.values, vals)
        end
    end
    return Float64(data["time"])
end

"""
    run_solver!(model, config; kwargs...)

Call the multiphase solver with only the kwargs the INSTALLED solver supports;
unsupported kwargs are dropped with a warning instead of a MethodError, so the
same case TOML works as solver features come and go.
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
