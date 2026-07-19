# Tatsumoto (2010) cryogenic N2 heat-transfer campaign

Validation of XCALibre's multiphase solver against Tatsumoto et al. (2010):
heated vertical quarter-pipe, liquid/supercritical N2.

| Case | Config | Physics |
|---|---|---|
| Supercritical | `configs/supercritical.toml` | 3.5 MPa, T_in = 78 K, through the pseudo-critical line |
| Subcooled boiling | `configs/boiling.toml` | 0.5 MPa, T_in = 77.84 K, RPI wall partition + Lee phase change |

## Layout

```
cases/tatsumoto/
├── configs/            case definitions (TOML) + hpc.toml (SLURM/cluster config)
├── stages/             prep.jl → heated.jl → post.py (+ boiling_* variants,
│                       steady_state_analysis.py, stage_common.jl helpers)
├── pipeline/submit.jl  DAG submitter: single run or parameter sweep, SLURM or local
└── runs/<name>/        ALL outputs (gitignored): case.toml (materialised),
                        warmup_t*.jld2, heated_stop*.jld2, vtk/, validation.png,
                        slurm/ (job scripts + logs)
```

Every run — single or sweep variant — owns one folder under `runs/`. Nothing is
written anywhere else; deleting `runs/` resets the campaign.

## Run locally

```bash
# whole chain (prep skipped automatically if its checkpoint exists):
julia --project=. cases/tatsumoto/pipeline/submit.jl cases/tatsumoto/configs/supercritical.toml --local

# individual stages:
julia --project=. cases/tatsumoto/stages/prep.jl   cases/tatsumoto/configs/supercritical.toml
julia --project=. cases/tatsumoto/stages/heated.jl cases/tatsumoto/configs/supercritical.toml
julia --project=. cases/tatsumoto/stages/heated.jl runs/<name>/case.toml resume   # from stop_T checkpoint
python cases/tatsumoto/stages/post.py              cases/tatsumoto/configs/supercritical.toml
python cases/tatsumoto/stages/steady_state_analysis.py <case.toml or --dir runs/<name>/vtk>
```

## Run on HPC (SLURM)

One-time setup (login node):
```bash
# Julia 1.11.3 via juliaup (the module tree only has 1.10):
curl -fsSL https://install.julialang.org | sh -s -- --yes
juliaup add 1.11.3 && juliaup default 1.11.3

git clone --branch AA/hpc --single-branch https://github.com/egyaa20/XCALibre.jl.git
cd XCALibre.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# fill in partition (and account, modules if needed):
$EDITOR cases/tatsumoto/configs/hpc.toml
```

Submit:
```bash
# inspect first — writes runs/<name>/slurm/*.sbatch + prints the DAG, submits nothing:
julia --project=. cases/tatsumoto/pipeline/submit.jl cases/tatsumoto/configs/supercritical.toml --dry-run

# submit for real:
julia --project=. cases/tatsumoto/pipeline/submit.jl cases/tatsumoto/configs/supercritical.toml

# monitor / cancel:
squeue -u $USER
scancel <jobid>
tail -f cases/tatsumoto/runs/<name>/slurm/*-<jobid>.out
```

The DAG per variant is `prep → heated (→ post)` chained with
`sbatch --dependency=afterok`. Post-processing is off by default on the cluster
(`[post] enable=false` in hpc.toml — it needs python3 + pyvista); the usual flow
is syncing results home and running `post.py` locally:
```bash
rsync -av hpc:XCALibre.jl/cases/tatsumoto/runs/ cases/tatsumoto/runs/
python cases/tatsumoto/stages/post.py cases/tatsumoto/runs/<name>/case.toml
```

## Parameter sweeps (batch mode)

Uncomment/add a `[batch]` section in the case TOML:

```toml
[batch]
mode = "product"                      # every combination; "zip" pairs lists 1:1
[batch.vary]
"energy.prop_relax"  = [0.3, 0.4, 0.6]
"run.adaptive.maxCo" = [0.9, 5.0]
```

`submit.jl` expands this into variants named
`<case>__prop_relax-0p3__maxCo-0p9`, …, each with a materialised, self-contained
`runs/<variant>/case.toml`. Variants whose varied keys only touch heated-stage
sections (`energy`, `heating`, `run`, `boiling`, `saturation`) **share one
baseline warmup job**; variants touching `mesh/flow/thermo/hardware/phase` or
`run.warmup_end`/`run.time_step` get their own prep, all wired automatically
into the SLURM dependency graph. `--local` runs the same expansion sequentially.

## Notes / current limitations

- The multiphase solver is being rebuilt feature-by-feature. Stage scripts drop
  solver options the installed solver doesn't support yet (with a warning), so
  these configs stay valid throughout. Currently pending re-add: the enthalpy
  (`VariableSensibleEnthalpy`) energy formulation, live table-driven properties,
  PIMPLE outer correctors, and the `stop_T` checkpoint hook — use
  `formulation = "temperature"` until then. The V&V platform
  (`test/verification/`) guards each re-addition.
- Grids: the `quarter_pipe_*` UNV files are committed gzipped
  (`examples/0_GRIDS/*.unv.gz`); stages auto-decompress on first use (needs
  `gzip`, present on any Linux/HPC).
- `runs/` is gitignored — checkpoints and VTU series never enter git.
