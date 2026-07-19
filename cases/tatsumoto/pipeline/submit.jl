# =============================================================================
# DAG SLURM submission pipeline for the tatsumoto campaign.
#
#   julia --project=. cases/tatsumoto/pipeline/submit.jl <case.toml> [options]
#
# Options:
#   --hpc <hpc.toml>     cluster config          (default: configs/hpc.toml)
#   --dry-run            write job scripts + print the DAG, submit nothing
#   --local              run stages sequentially on this machine (no SLURM)
#   --stages a,b,c       subset of prep,heated,post (default: prep,heated
#                        [+post if enabled in hpc.toml / always offered locally])
#   --reprep             force warmup even if a checkpoint already exists
#
# Behaviour:
#   * No [batch] section  -> single variant, DAG: prep -> heated (-> post).
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

const CAMPAIGN = dirname(@__DIR__)                       # cases/tatsumoto
const STAGES   = joinpath(CAMPAIGN, "stages")
const RUNS     = joinpath(CAMPAIGN, "runs")
const REPO     = dirname(dirname(CAMPAIGN))              # repository root

# Varied keys with these prefixes invalidate the shared warmup:
const WARMUP_PREFIXES = ("mesh.", "flow.", "thermo.", "hardware.", "phase.")
const WARMUP_KEYS     = ("run.warmup_end", "run.time_step")

# ----------------------------------------------------------------- utilities
function parse_args(argv)
    opts = Dict{String,Any}("stages" => nothing, "hpc" => joinpath(CAMPAIGN, "configs", "hpc.toml"),
                            "dry" => false, "local" => false, "reprep" => false, "case" => nothing)
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--dry-run";      opts["dry"] = true
        elseif a == "--local";    opts["local"] = true
        elseif a == "--reprep";   opts["reprep"] = true
        elseif a == "--hpc";      opts["hpc"] = abspath(argv[i+=1])
        elseif a == "--stages";   opts["stages"] = strip.(split(argv[i+=1], ","))
        elseif startswith(a, "--"); error("unknown option $a")
        else opts["case"] = abspath(a)
        end
        i += 1
    end
    opts["case"] === nothing && error(
        "usage: julia pipeline/submit.jl <case.toml> [--hpc hpc.toml] [--dry-run] [--local] [--stages prep,heated,post] [--reprep]")
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
            run_local_stage(hpc, prep_stage, btoml)
        else
            js = write_job(hpc, "prep", base_name, bdir, prep_stage, btoml;
                           time_override=hpc["slurm"]["time_prep"])
            base_prep_id = sbatch!(js; dry=opts["dry"])
            println("[baseline]  prep job $base_prep_id  ($js)")
        end
    end

    for (suffix, overrides, own_prep) in variants
        name = isempty(suffix) ? base_name : "$(base_name)__$(suffix)"
        wfrom = (!own_prep && share_warmup) ? base_name : nothing
        vtoml = materialise(cfg, name, overrides, hpc; warmup_from=wfrom)
        vdir = joinpath(RUNS, name)
        wend = Float64(cfg["run"]["warmup_end"])
        own_ckpt = joinpath(vdir, "warmup_t$(wend).jld2")

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
                prep_id = sbatch!(js; dry=opts["dry"])
                println("[$name]  prep   $prep_id")
            end
        end
        heat_id = nothing
        if want("heated")
            js = write_job(hpc, "heated", name, vdir, heated_stage, vtoml;
                           time_override=hpc["slurm"]["time_heated"])
            heat_id = sbatch!(js; dep=prep_id, dry=opts["dry"])
            println("[$name]  heated $heat_id" *
                    (prep_id === nothing ? "" : "  (afterok:$prep_id)"))
        end
        if do_post
            js = write_job(hpc, "post", name, vdir, "", vtoml)
            pid = sbatch!(js; dep=heat_id, dry=opts["dry"])
            println("[$name]  post   $pid" *
                    (heat_id === nothing ? "" : "  (afterok:$heat_id)"))
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

function run_local_post(hpc, case_toml)
    py = hpc["post"]["python"]
    run(`$py $(joinpath(STAGES, "post.py")) $case_toml`)
end

main(ARGS)
