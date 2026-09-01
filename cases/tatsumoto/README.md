# Tatsumoto (2010) cryogenic N2 heat-transfer campaign

Validation of the multiphase solver (`Solvers_5_Multiphase.jl`) against
Tatsumoto et al. (2010): heated vertical quarter-pipe, supercritical N2 at
3.5 MPa, T_in = 78 K, U_in = 3.7 m/s, through the pseudo-critical line.

Physics stack: `HelmholtzTable` real-fluid properties (rho/mu/k/cp/h = f(T) at
p_ref, with T<->h inversion), `Energy{HelmholtzEnthalpy}` enthalpy equation,
compressible-form continuity, `HeatFlux`/`HeatFluxFunction` wall BCs (exact
q*A injection — the applied flux never depends on local lambda), KOmegaSST.

## Layout

```
cases/tatsumoto/
├── configs/            supercritical.toml + hpc.toml (SLURM/cluster config)
├── mesh/               quarter_pipe.py (SALOME, runs on the cluster)
├── stages/             prep.jl → heated.jl → post.py
│                       (stage_common.jl + case_model.jl shared helpers)
├── pipeline/submit.jl  DAG submitter: single run or parameter sweep, SLURM or local
└── runs/<name>/        ALL outputs (gitignored): case.toml (materialised),
                        warmup_t*.jld2, heated_final.jld2, vtk/, validation.png,
                        slurm/ (job scripts + logs)
```

## Stages

1. **mesh** — SALOME builds the quarter-pipe grid on the cluster from
   `[mesh.geometry|resolution|boundary_layer]` into a content-addressed cache.
2. **prep** — flow development with the constant baseline flux Q0 until
   `run.warmup_end`; saves `warmup_t*.jld2`.
3. **heated** — restarts from the warmup checkpoint, applies the ramp
   `q(t) = Q0*exp(t/tau)` on `Wall_Heated` via `HeatFluxFunction` until
   `run.heated_end`; writes `vtk/time_*.vtu` + `times.csv` and saves
   `heated_final.jld2` on completion (reusable via `heated.jl <case> resume`).
4. **post** — `post.py` extracts peak wall ΔT per snapshot and plots q vs ΔT_L
   against the digitised Tatsumoto 3.5 MPa reference into
   `summary/validation.png` (published to pipeline-status automatically).

## Run locally

```bash
julia --project=. cases/tatsumoto/pipeline/submit.jl cases/tatsumoto/configs/supercritical.toml --local

# individual stages:
julia --project=. cases/tatsumoto/stages/prep.jl   cases/tatsumoto/configs/supercritical.toml
julia --project=. cases/tatsumoto/stages/heated.jl cases/tatsumoto/configs/supercritical.toml
julia --project=. cases/tatsumoto/stages/heated.jl runs/<name>/case.toml resume
python cases/tatsumoto/stages/post.py              cases/tatsumoto/configs/supercritical.toml
```

## Run on HPC (SLURM)

```bash
# inspect first — writes runs/<name>/slurm/*.sbatch + prints the DAG, submits nothing:
julia --project=. cases/tatsumoto/pipeline/submit.jl cases/tatsumoto/configs/supercritical.toml --dry-run

# submit for real (or push to the queue branch and let the poller do it):
julia --project=. cases/tatsumoto/pipeline/submit.jl cases/tatsumoto/configs/supercritical.toml
```

The DAG per variant is `mesh → prep → heated → post` chained with
`sbatch --dependency=afterok`. Resources come from `configs/hpc.toml`
(16 CPU cores per solver job).

## Parameter sweeps (batch mode)

Add a `[batch]` section to the case TOML:

```toml
[batch]
mode = "product"                      # every combination; "zip" pairs lists 1:1
[batch.vary]
"run.adaptive.maxCo" = [0.5, 0.9]
```

Variants whose varied keys only touch heated-stage sections share one baseline
warmup job; variants touching `mesh/flow/thermo/hardware` or
`run.warmup_end`/`run.time_step`/`run.adaptive.maxCo` get their own prep.

## Notes

- `test/smoke.jl` is the validate.sh gate: table build at 3.5 MPa, HeatFlux BC,
  JLD2 checkpoint roundtrip, short multiphase solve with t_end.
- Solver options unsupported by the installed solver are dropped by
  `run_solver!` with a warning, so configs stay valid as features evolve.
- `runs/` is gitignored — checkpoints and VTU series never enter git.
