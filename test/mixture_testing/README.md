# Mixture solver validation suite

Systematic validation of `src/Solvers/Solvers_2_MULTIPHASE.jl` (Manninen
drift-flux mixture model) against analytical closures.

## Three-stage hierarchy

| Stage | File | What it validates | CFD? |
|-------|------|-------------------|------|
| 1 | `stage1_terminal_velocity.jl` | Drag closure, kernel correctness, field blending | No |
| 2 | `stage2_hindered_settling.jl` | Richardson-Zaki hindered settling `u_s = u_t·(1-α)^n` | Yes |
| 3 | `stage3_batch_sedimentation.jl` | 2D batch settling, interface descent rate | Yes |

## Reference physics

All three stages share the same canonical setup:

- **Continuous phase** (tracked α): water, ρ = 1000 kg/m³, μ = 1.0e-3 Pa·s
- **Dispersed phase**: glass beads, ρ = 2500 kg/m³, μ = 1.0e-3 Pa·s
- **Particle diameter**: d = 200 µm
- **Virtual-mass coefficient**: C_vm = 0.5
- **Gravity**: 9.81 m/s² along -y

**Hand calculations (dilute, α_d → 0)**
- Stokes: `u_t = g(ρ_d-ρ_c)d² / (18μ_c) ≈ 33 mm/s`
- Schiller-Naumann (Manninen buoyancy form): `u_t ≈ 21 mm/s`
- Schiller-Naumann (**solver's** buoyancy form `(ρ_d-ρ_m)/(ρ_m+C_vm·ρ_c)`):
  `u_t ≈ 35 mm/s` — note the ~1.7× factor vs standard Manninen.

## Required edits to `Solvers_2_MULTIPHASE.jl`

Near line ~180 of the solver, replace the default bubble diameter with the
200 µm particle diameter used by these tests:

```julia
rho1_val = phases[main].rho[1]
rho2_val = phases[secondary].rho[1]
mu1_val  = phases[main].mu[1]
diameter = 2.0e-4      # was 3.0e-3 (3 mm bubbles) — change to 200 µm
C_vm     = 0.5
Sc_t     = 0.7
```

Running the stage files **without** this edit will produce terminal velocities
consistent with 3 mm bubbles in water (~6× larger), which is the expected
behavior of the solver-as-written.

## Running

From the package root:

```julia
# Stage 1 — no CFD, fast (< 5 s)
julia --project=. test/mixture_testing/stage1_terminal_velocity.jl

# Stage 2 — short CFD run (minutes)
julia --project=. test/mixture_testing/stage2_hindered_settling.jl

# Stage 3 — longer CFD run (tens of minutes)
julia --project=. test/mixture_testing/stage3_batch_sedimentation.jl
```

## Sweeps and extensions

Stage 3 can be run multiple times with different `α_d0` and `d_part` to
verify scaling laws:

- `α_d0 ∈ {0.05, 0.10, 0.15, 0.20}` → verify `(1-α_d)^4.65` scaling
- `d ∈ {50, 100, 200, 500} µm`      → verify `d²` at small d; SN roll-off at large d

A reduction of `dH/dt` across the sweep that does **not** follow
Richardson-Zaki indicates either a bug in the slip-velocity path or — more
likely — the buoyancy-form discrepancy noted above.

## Interpretation notes

- If Stage 1 fails → kernel-level problem (drag closure, blending, gradient).
- If Stage 1 passes but Stage 2 produces `mean(α)` drift → α-equation /
  drift-flux divergence inconsistency.
- If Stages 1–2 pass but Stage 3 has wrong slope → either the buoyancy form
  (1.5–1.7× offset) or interface smearing from the upwind α scheme.
