# =============================================================================
# DAG SLURM submission pipeline for the tatsumoto campaign.
#
#   julia --project=. cases/tatsumoto/pipeline/submit.jl <case.toml> [options]
#
# Options:
#   --hpc <hpc.toml>     cluster config          (default: configs/hpc.toml)
#   --dry-run            write job scripts + print the DAG, submit nothing
#   --local              run stages sequentially on this machine (no SLURM)
#   --stages a,b,c       subset of mesh,prep,heated,post (default: mesh,prep,
#                        heated [+post if enabled in hpc.toml / always locally])
#   --reprep             force warmup even if a checkpoint already exists
#   --job-report <path>  append "<name>-<stage> <jobid>" lines here as each
#                        sbatch call succeeds. Consumed by hpc/poller.sh's
#                        hook mode (CFD_JOB_REPORT) so a `git push` alone can
#                        drive submission with no interactive SSH/VPN needed
#                        -- see pipeline/on_update.sh. A no-op if omitted.
#
# Behaviour:
#   * No [batch] section  -> single variant, DAG: mesh -> prep -> heated (-> post).
#   * [mesh.geometry] present -> the mesh is BUILT by SALOME on the cluster
#     (singularity container, see [mesh] in hpc.toml) into a content-addressed
#     cache under runs/_meshes/, keyed by a hash of the mesh parameters. An
#     already-cached mesh is reused and no mesh job is submitted. Without
#     [mesh.geometry], [mesh] grid = "<name>" selects a committed grid as before.
#   * With [batch]        -> one variant per parameter combination
#     (mode="product") or per zipped tuple (mode="zip"). Variants whose varied
#     keys only touch heated-stage sections share ONE baseline warmup job;
#     variants touching mesh/flow/thermo/hardware/phase or run.warmup_end /
#     run.time_step get their own prep job. Dependencies via
#     `sbatch --dependency=afterok`.
#   * Each variant gets runs/<name>/: case.toml (materialised, self-contained),
#     slurm/ (job scripts + logs), warmup_t*.jld2, vtk/, post outputs.
#
# Julia stdlib only (TOML) — no Python needed to submit.
# =============================================================================

using TOML
using SHA

const CAMPAIGN = dirname(@__DIR__)                       # cases/tatsumoto
const STAGES   = joinpath(CAMPAIGN, "stages")
const RUNS     = joinpath(CAMPAIGN, "runs")
const MESHES   = get(ENV, "CFD_MESH_CACHE",              # shared mesh cache;
                     joinpath(RUNS, "_meshes"))          # env override makes it
                                                         # survive hook snapshots
const REPO     = dirname(dirname(CAMPAIGN))              # repository root
const MESH_PY  = joinpath(CAMPAIGN, "mesh", "quarter_pipe.py")

# Parameters understood by the SALOME generator, in the order it documents
# them. Anything else in [mesh.*] is a submitter concern (scale, grid, ...).
const MESH_PARAM_KEYS = ("radius", "L_entrance", "L_heated", "L_exit", "core_ratio",
                         "n_quarter", "n_radial", "n_axial",
                         "n_ax_entrance", "n_ax_heated", "n_ax_exit", "axial_ratio",
                         "first_cell", "bl_growth", "bl_flip")

# Varied keys with these prefixes invalidate the shared warmup:
const WARMUP_PREFIXES = ("mesh.", "flow.", "thermo.", "hardware.", "phase.")
const WARMUP_KEYS     = ("run.warmup_end", "run.time_step", "run.adaptive.maxCo")
# ^ maxCo governs the prep (warmup) stage's own adaptive dt just as much as
# heated's -- varying it without giving each variant its own prep would mean
# all variants share ONE warmup at the baseline maxCo and only differ once
# heated starts, silently defeating a sweep whose entire point is comparing
# warmup-stage timestep growth across maxCo values.

# ----------------------------------------------------------------- utilities
function parse_args(argv)
    opts = Dict{String,Any}("stages" => nothing, "hpc" => joinpath(CAMPAIGN, "configs", "hpc.toml"),
                            "dry" => false, "local" => false, "reprep" => false, "case" => nothing,
                            "job_report" => nothing)
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--dry-run";      opts["dry"] = true
        elseif a == "--local";    opts["local"] = true
        elseif a == "--reprep";   opts["reprep"] = true
        elseif a == "--hpc";      opts["hpc"] = abspath(argv[i+=1])
        elseif a == "--stages";   opts["stages"] = strip.(split(argv[i+=1], ","))
        elseif a == "--job-report"; opts["job_report"] = abspath(argv[i+=1])
        elseif startswith(a, "--"); error("unknown option $a")
        else opts["case"] = abspath(a)
        end
        i += 1
    end
    opts["case"] === nothing && error(
        "usage: julia pipeline/submit.jl <case.toml> [--hpc hpc.toml] [--dry-run] [--local] [--stages prep,heated,post] [--reprep] [--job-report path]")
    return opts
end

setdeep!(d::AbstractDict, key::AbstractString, v) = begin
    ks = split(key, ".")
    cur = d
    for k in ks[1:end-1]
        cur = get!(cur, String(k), Dict{String,Any}())
    end
    cur[String(ks[end])] = v
end

sanitize(x) = replace(string(x), "." => "p", "-" => "m", "/" => "_", " " => "")
short_key(k) = String(last(split(k, ".")))

is_warmup_key(k) = any(startswith(k, p) for p in WARMUP_PREFIXES) || k in WARMUP_KEYS

# --------------------------------------------------------------- meshing
"True when [mesh] describes a mesh to generate rather than a committed grid."
is_generated_mesh(mesh_cfg) = haskey(mesh_cfg, "geometry")

"Flatten [mesh.geometry|resolution|boundary_layer] into the generator's params."
function mesh_params(mesh_cfg)
    p = Dict{String,Any}()
    for sect in ("geometry", "resolution", "boundary_layer")
        haskey(mesh_cfg, sect) || continue
        for (k, v) in mesh_cfg[sect]
            k in MESH_PARAM_KEYS ||
                error("unknown mesh parameter [mesh.$sect] $k — valid: $(join(MESH_PARAM_KEYS, ", "))")
            p[k] = v
        end
    end
    return p
end

"""
Canonical JSON for the SALOME generator.

Keys are emitted in a fixed order with a fixed float format so the text is a
function of the values alone — the cache hash is taken over this string, and
an incidental reordering in the TOML must not silently invalidate a mesh.
`first_cell` is emitted as null when absent, which is how the generator is
told to fall back to `bl_growth`.
"""
function mesh_params_json(p)
    fmt(v) = v isa Bool ? (v ? "true" : "false") :
             v isa Integer ? string(v) : string(float(v))
    lines = String[]
    for k in MESH_PARAM_KEYS
        haskey(p, k) || continue
        push!(lines, "  \"$k\": $(fmt(p[k]))")
    end
    haskey(p, "first_cell") || push!(lines, "  \"first_cell\": null")
    return "{\n" * join(lines, ",\n") * "\n}\n"
end

"Cache paths for a mesh spec: (unv, params.json, 12-char hash)."
function mesh_cache_paths(mesh_cfg)
    json = mesh_params_json(mesh_params(mesh_cfg))
    h = bytes2hex(sha256(json))[1:12]
    return (joinpath(MESHES, "quarter_pipe-$h.unv"),
            joinpath(MESHES, "quarter_pipe-$h.params.json"),
            h, json)
end

"Write the params file for a generated mesh and return its .unv cache path."
function stage_mesh_params(mesh_cfg)
    unv, jsonf, h, json = mesh_cache_paths(mesh_cfg)
    mkpath(MESHES)
    write(jsonf, json)
    return unv, jsonf, h
end

function write_mesh_job(hpc, name, run_dir, unv, jsonf)
    m = hpc["mesh"]
    s = hpc["slurm"]
    lines = String[
        "#!/bin/bash",
        "#SBATCH --job-name=$(name)-mesh",
        "#SBATCH --partition=$(s["partition"])",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=$(get(m, "cpus_per_task", 4))",
        "#SBATCH --mem=$(get(m, "mem", "16G"))",
        "#SBATCH --time=$(get(m, "time", "02:00:00"))",
        "#SBATCH --output=$(joinpath(run_dir, "slurm"))/%x-%j.out",
    ]
    isempty(s["account"]) || push!(lines, "#SBATCH --account=$(s["account"])")
    push!(lines, "")
    for mod in get(m, "modules", String[])
        push!(lines, "module load $mod")
    end
    append!(lines, [
        "set -euo pipefail",
        # Idempotent: the cache key is the mesh spec, so an existing file is
        # by definition the mesh this job would have produced. Re-running the
        # DAG (or a sibling variant) must not rebuild it.
        "if [[ -f '$unv' ]]; then echo \"mesh already built: $unv\"; exit 0; fi",
        "export SALOME_MESH_PARAMS='$jsonf'",
        "export SALOME_MESH_OUT='$unv.tmp'",
        "singularity exec '$(m["container"])' salome -t '$MESH_PY'",
        # Publish atomically: a half-written .unv left by a timeout would
        # otherwise look like a valid cache hit to every later run. The stats
        # sidecar needs no move — the generator derives its name by stripping
        # one extension, so writing to <unv>.tmp already lands it at
        # <unv>.mesh.json, next to where the mesh itself is about to appear.
        "mv '$unv.tmp' '$unv'",
        "echo \"mesh written: $unv\"",
    ])
    path = joinpath(run_dir, "slurm", "mesh.sbatch")
    write(path, join(lines, "\n") * "\n")
    return path
end

# ------------------------------------------------------- batch expansion
"Expand [batch] into (name_suffix, overrides::Dict, own_prep::Bool) tuples."
function expand_batch(cfg)
    haskey(cfg, "batch") || return [("", Dict{String,Any}(), true)]
    b = cfg["batch"]
    vary = get(b, "vary", Dict{String,Any}())
    isempty(vary) && return [("", Dict{String,Any}(), true)]
    keys_sorted = sort(collect(keys(vary)))
    lists = [vary[k] isa AbstractVector ? vary[k] : [vary[k]] for k in keys_sorted]
    mode = lowercase(get(b, "mode", "product"))

    combos = if mode == "zip"
        n = length(first(lists))
        all(length(l) == n for l in lists) ||
            error("[batch] mode=\"zip\" requires equal-length lists")
        [[l[i] for l in lists] for i in 1:n]
    else # product
        result = [Any[]]
        for l in lists
            result = [vcat(r, [v]) for r in result for v in l]
        end
        result
    end

    variants = Vector{Tuple{String,Dict{String,Any},Bool}}()
    for combo in combos
        ov = Dict{String,Any}(k => v for (k, v) in zip(keys_sorted, combo))
        suffix = join(["$(short_key(k))-$(sanitize(v))" for (k, v) in zip(keys_sorted, combo)], "__")
        own_prep = any(is_warmup_key(k) for k in keys_sorted)
        push!(variants, (suffix, ov, own_prep))
    end
    return variants
end

# --------------------------------------------------- variant materialisation
"Deep-copy cfg, apply overrides, write runs/<name>/case.toml. Returns its path."
function materialise(cfg, name, overrides, hpc; warmup_from=nothing)
    v = deepcopy(cfg)
    haskey(v, "batch") && delete!(v, "batch")
    for (k, val) in overrides
        setdeep!(v, k, val)
    end
    v["case"]["name"] = name
    run_dir = joinpath(RUNS, name)
    v["case"]["run_dir"] = run_dir
    if !get(hpc["slurm"], "gpu", false)
        setdeep!(v, "hardware.backend", "CPU")
    end
    warmup_from === nothing || setdeep!(v, "run.warmup_from", warmup_from)
    # Resolve the mesh AFTER overrides: a batch variant that varies a mesh
    # parameter must hash its own value, not the baseline's.
    if is_generated_mesh(v["mesh"])
        unv, _, _ = stage_mesh_params(v["mesh"])
        v["mesh"]["generated_file"] = unv
    end
    mkpath(joinpath(run_dir, "slurm"))
    path = joinpath(run_dir, "case.toml")
    open(path, "w") do io
        TOML.print(io, v)
    end
    return path
end

# ------------------------------------------------------------- job scripts
function sbatch_script(hpc, jobname, run_dir, cmd; time_override="")
    s = hpc["slurm"]; e = hpc["env"]
    threads = e["threads"] == 0 ? s["cpus_per_task"] : e["threads"]
    lines = String[
        "#!/bin/bash",
        "#SBATCH --job-name=$jobname",
        "#SBATCH --partition=$(s["partition"])",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=$(s["cpus_per_task"])",
        "#SBATCH --mem=$(s["mem"])",
        "#SBATCH --time=$(isempty(time_override) ? s["time"] : time_override)",
        "#SBATCH --output=$(joinpath(run_dir, "slurm"))/%x-%j.out",
    ]
    isempty(s["account"]) || push!(lines, "#SBATCH --account=$(s["account"])")
    get(s, "gpu", false) && push!(lines, "#SBATCH --gres=gpu:1")
    push!(lines, "")
    for m in e["modules"]
        push!(lines, "module load $m")
    end
    push!(lines, "export JULIA_NUM_THREADS=$threads")
    push!(lines, "set -euo pipefail")
    push!(lines, cmd)
    return join(lines, "\n") * "\n"
end

function write_job(hpc, kind, name, run_dir, stage_file, case_toml; time_override="")
    e = hpc["env"]
    project = isempty(e["project"]) ? REPO : e["project"]
    cmd = if kind == "post"
        "$(hpc["post"]["python"]) $(joinpath(STAGES, "post.py")) $case_toml"
    else
        "$(e["julia"]) --project=$project --threads=\$JULIA_NUM_THREADS $(joinpath(STAGES, stage_file)) $case_toml"
    end
    script = sbatch_script(hpc, "$(name)-$(kind)", run_dir, cmd; time_override=time_override)
    path = joinpath(run_dir, "slurm", "$kind.sbatch")
    write(path, script)
    return path
end

sbatch!(script; dep=nothing, dry=false) = begin
    args = ["--parsable"]
    dep === nothing || push!(args, "--dependency=afterok:$dep")
    push!(args, script)
    dry && return "DRY"
    readchomp(`sbatch $args`)
end

"""
    report_job!(opts, label, jobid, run_dir)

Append "`label` `jobid` `run_dir`" to `opts["job_report"]` (if set) -- one
line per submitted job. `hpc/poller.sh`'s hook mode reads this back to
publish job state and notifications, exactly like the spec-file path
already does for single-job runs; this is what lets a DAG (mesh -> prep ->
heated, possibly one per batch variant) report itself the same way.

`run_dir` is each job's own absolute output directory (e.g.
cases/tatsumoto/runs/<variant>/), not the poller's hook-trigger snapshot
dir -- without it, every job from one hook trigger used to get tagged with
the shared snapshot dir, so the poller's summary/ artifact publishing
(validation.png, etc.) was looking in the wrong place for every hook job.

Silently skipped for "DRY" ids (dry-run) and when job_report isn't set
(interactive/local use).
"""
function report_job!(opts, label, jobid, run_dir)
    path = get(opts, "job_report", nothing)
    (path === nothing || jobid == "DRY") && return
    open(path, "a") do io
        println(io, "$label $jobid $run_dir")
    end
end

"""
    ensure_mesh!(mesh_ids, hpc, opts, case_toml, name, run_dir) -> jobid | nothing

Make sure this case's mesh exists, returning a job id to depend on (or
`nothing` if the mesh is already on disk / was submitted for an earlier
variant sharing the same spec). `mesh_ids` maps cache hash -> job id, which is
what stops a sweep of N solver variants from launching N identical mesh jobs.
"""
function ensure_mesh!(mesh_ids, hpc, opts, case_toml, name, run_dir)
    vcfg = TOML.parsefile(case_toml)
    is_generated_mesh(vcfg["mesh"]) || return nothing
    unv, jsonf, h = stage_mesh_params(vcfg["mesh"])
    if isfile(unv)
        println("[$name]  mesh   cached ($(basename(unv)))")
        return nothing
    end
    if haskey(mesh_ids, h)
        println("[$name]  mesh   shared with an earlier variant ($(mesh_ids[h]))")
        return mesh_ids[h]
    end
    if opts["local"]
        run_local_mesh(hpc, unv, jsonf)
        return nothing
    end
    js = write_mesh_job(hpc, name, run_dir, unv, jsonf)
    id = sbatch!(js; dry=opts["dry"])
    mesh_ids[h] = id
    println("[$name]  mesh   $id  ($(basename(unv)))")
    report_job!(opts, "$name-mesh", id, run_dir)
    return id
end

# ------------------------------------------------------------------- main
function main(argv)
    opts = parse_args(argv)
    cfg = TOML.parsefile(opts["case"])
    hpc = TOML.parsefile(opts["hpc"])
    base_name = cfg["case"]["name"]
    prep_stage   = get(cfg["case"], "prep_stage", "prep.jl")
    heated_stage = get(cfg["case"], "heated_stage", "heated.jl")
    do_post = opts["stages"] !== nothing ? ("post" in opts["stages"]) :
              get(hpc["post"], "enable", false)
    want(st) = opts["stages"] === nothing || st in opts["stages"]

    variants = expand_batch(cfg)
    batch = length(variants) > 1 || variants[1][1] != ""
    share_warmup = batch && any(!own for (_, _, own) in variants)

    # One mesh job per distinct mesh spec, keyed by cache hash: variants that
    # differ only in solver settings share a mesh and must not each rebuild it.
    mesh_ids = Dict{String,Any}()

    println("="^78)
    println("Campaign: $base_name   variants: $(length(variants))" *
            (share_warmup ? "   (shared baseline warmup)" : ""))
    println("="^78)

    # Shared baseline warmup job (baseline parameters, prep only)
    base_prep_id = nothing
    if share_warmup && want("prep")
        btoml = materialise(cfg, base_name, Dict{String,Any}(), hpc)
        bdir = joinpath(RUNS, base_name)
        wend = Float64(cfg["run"]["warmup_end"])
        ckpt = joinpath(bdir, "warmup_t$(wend).jld2")
        if isfile(ckpt) && !opts["reprep"]
            println("[baseline]  warmup checkpoint exists -> shared prep skipped")
        elseif opts["local"]
            want("mesh") && ensure_mesh!(mesh_ids, hpc, opts, btoml, "baseline", bdir)
            run_local_stage(hpc, prep_stage, btoml)
        else
            bmesh = want("mesh") ?
                ensure_mesh!(mesh_ids, hpc, opts, btoml, "baseline", bdir) : nothing
            js = write_job(hpc, "prep", base_name, bdir, prep_stage, btoml;
                           time_override=hpc["slurm"]["time_prep"])
            base_prep_id = sbatch!(js; dep=bmesh, dry=opts["dry"])
            println("[baseline]  prep job $base_prep_id  ($js)" *
                    (bmesh === nothing ? "" : "  (afterok:$bmesh)"))
            report_job!(opts, "$base_name-prep", base_prep_id, bdir)
        end
    end

    for (suffix, overrides, own_prep) in variants
        name = isempty(suffix) ? base_name : "$(base_name)__$(suffix)"
        wfrom = (!own_prep && share_warmup) ? base_name : nothing
        vtoml = materialise(cfg, name, overrides, hpc; warmup_from=wfrom)
        vdir = joinpath(RUNS, name)
        wend = Float64(cfg["run"]["warmup_end"])
        own_ckpt = joinpath(vdir, "warmup_t$(wend).jld2")

        # --- mesh (stage 0): build this variant's mesh if it isn't cached ---
        mesh_id = want("mesh") ?
            ensure_mesh!(mesh_ids, hpc, opts, vtoml, name, vdir) : nothing

        if opts["local"]
            println("-"^78); println("[local] $name")
            if own_prep && want("prep") && (!isfile(own_ckpt) || opts["reprep"])
                run_local_stage(hpc, prep_stage, vtoml)
            end
            want("heated") && run_local_stage(hpc, heated_stage, vtoml)
            do_post && run_local_post(hpc, vtoml)
            continue
        end

        prep_id = base_prep_id
        if own_prep && want("prep")
            if isfile(own_ckpt) && !opts["reprep"]
                println("[$name]  warmup checkpoint exists -> prep skipped")
                prep_id = nothing
            else
                js = write_job(hpc, "prep", name, vdir, prep_stage, vtoml;
                               time_override=hpc["slurm"]["time_prep"])
                prep_id = sbatch!(js; dep=mesh_id, dry=opts["dry"])
                println("[$name]  prep   $prep_id" *
                        (mesh_id === nothing ? "" : "  (afterok:$mesh_id)"))
                report_job!(opts, "$name-prep", prep_id, vdir)
            end
        end
        heat_id = nothing
        if want("heated")
            js = write_job(hpc, "heated", name, vdir, heated_stage, vtoml;
                           time_override=hpc["slurm"]["time_heated"])
            # Skipped prep still leaves a mesh dependency to honour: heated
            # reads the mesh too, so it must not start before the mesh lands.
            hdep = prep_id === nothing ? mesh_id : prep_id
            heat_id = sbatch!(js; dep=hdep, dry=opts["dry"])
            println("[$name]  heated $heat_id" *
                    (hdep === nothing ? "" : "  (afterok:$hdep)"))
            report_job!(opts, "$name-heated", heat_id, vdir)
        end
        if do_post
            js = write_job(hpc, "post", name, vdir, "", vtoml;
                           time_override=hpc["slurm"]["time_post"])
            pid = sbatch!(js; dep=heat_id, dry=opts["dry"])
            println("[$name]  post   $pid" *
                    (heat_id === nothing ? "" : "  (afterok:$heat_id)"))
            report_job!(opts, "$name-post", pid, vdir)
        end
    end

    println("="^78)
    opts["dry"] && println("DRY RUN — job scripts written under runs/<name>/slurm/, nothing submitted.")
    opts["local"] && println("Local run complete.")
end

function run_local_stage(hpc, stage_file, case_toml)
    e = hpc["env"]
    project = isempty(e["project"]) ? REPO : e["project"]
    jl = e["julia"]
    println("  \$ $jl --project=$project $(joinpath(STAGES, stage_file)) $case_toml")
    run(`$jl --project=$project $(joinpath(STAGES, stage_file)) $case_toml`)
end

function run_local_mesh(hpc, unv, jsonf)
    m = hpc["mesh"]
    println("  \$ singularity exec $(m["container"]) salome -t $MESH_PY")
    withenv("SALOME_MESH_PARAMS" => jsonf, "SALOME_MESH_OUT" => unv * ".tmp") do
        run(`singularity exec $(m["container"]) salome -t $MESH_PY`)
    end
    mv(unv * ".tmp", unv; force=true)   # stats sidecar already lands correctly
end

function run_local_post(hpc, case_toml)
    py = hpc["post"]["python"]
    run(`$py $(joinpath(STAGES, "post.py")) $case_toml`)
end

main(ARGS)
