# VOF solver verification suite

Systematic verification of the VOF branch of `src/Solvers/Solvers_2_MULTIPHASE.jl`
against the modified-equation analysis of the MULES-advected indicator field.

## Modified equation — what the scheme actually solves

The continuous α-transport equation is

```
∂α/∂t + ∇·(Uα) = 0
```

The discrete update used here is forward-Euler in time with a MULES-blended
flux: low-order upwind plus a limited anti-diffusive (compression) flux.
Taylor-expanding the discrete update on a uniform 1D grid gives the
*modified equation* the scheme actually solves, to leading order:

```
∂α/∂t + ∇·(Uα) = ½ |U| Δx (1 − C)(1 − λ_M) ∂²α/∂x²              (numerical
                                                                   diffusion)
                 − λ_M C_α |U| ∂/∂x [ α(1−α) n̂ ]                   (compression
                                                                   anti-diffusion)
                 + 𝒪(Δx², Δt²)
```

where
- `C = U Δt / Δx` is the cell Courant number,
- `λ_M ∈ [0, 1]` is the FCT (MULES) limiter,
- `C_α` is the compression coefficient (`cAlpha` in the API),
- `n̂ = ∇α / |∇α|` is the interface normal.

Concrete predictions:

1. **Without compression** (`cAlpha = 0`, `λ_M ≡ 0`): the scheme is plain
   upwind. Effective diffusivity `D_num = ½ |U| Δx (1 − C)`. A planar interface
   smears as an error-function profile with width `σ(t) = √(2 D_num · t)`.

2. **With balanced compression** (`cAlpha = 1`): inside the interface band
   (where `α(1−α) > 0`), the anti-diffusive flux formally cancels the
   leading-order numerical diffusion. The interface settles to a steady
   ~1.5–2 cell width and stays there.

3. **Without the FCT limiter**, `cAlpha > 1` would produce overshoots
   (negative effective diffusion → ill-posed). With MULES, `λ_M` dynamically
   suppresses overshoots while permitting sharper-than-balance compression.
   The result is interface preserved at ~1 cell with no overshoot.

4. **Outside the interface band** (`α ≡ 0` or `α ≡ 1`), the limiter forces
   `λ_M = 0`, so plain upwind diffusion is recovered. Smooth gas/liquid
   bulk regions don't get spuriously sharpened by compression.

5. **Convergence**: the scheme is between 1st-order (away from interface,
   pure upwind) and 2nd-order (in the interface band, where compression
   cancels the leading-order error). Empirically, mesh refinement gives a
   slope between 1 and 2 in `log(error)` vs `log(Δx)`.

This README states (1)–(5) as predictions; the stage files below verify
each empirically.

## Important caveat — XCALibre's MULES is *not* simple upwind+compression

The LLM's modified-equation analysis assumes the scheme blends pure upwind
with a compression flux. **XCALibre's MULES blends pure upwind with a
van-Leer HO flux** (and additionally adds a compression term inside that HO
flux when `cAlpha > 0`). At `cAlpha = 0` the scheme is therefore *not* plain
upwind — it's already 2nd-order in the interface band via the limited HO
contribution. The pure-upwind modified equation `σ ≈ √(|U|Δx t)` does
**not** apply directly to the solver.

The verification stages below are therefore framed around what the solver
*actually* does:

- Compare the measured interface width to the pure-upwind reference and
  show it is *much sharper* (HO is active).
- Show that adding `cAlpha = 1` further reduces width to ~1 Δx.
- Verify bound preservation (`α ∈ [0, 1]`) at all stages.
- Measure the empirical convergence rate, which should fall between 1 and 2.

For a clean test of the pure-upwind modified-equation prediction, the HO
flux assembly in `Solvers_2_MULTIPHASE.jl` would need a temporary toggle
that replaces `alphaf_HO` with `alphaf_upwind` (so MULES degenerates to
pure upwind). Not done in the current solver.

## Three-stage hierarchy

| Stage | File | Verifies | Quantity measured |
|-------|------|----------|--------------------|
| 1 | `stage1_upwind_diffusion.jl` | HO is active at `cAlpha = 0` | Interface width « pure-upwind prediction |
| 2 | `stage2_compression_balance.jl` | Compression sharpens interface to ~1 Δx | Interface width with `cAlpha = 1` |
| 3 | `stage3_convergence.jl` | Empirical convergence rate (1–2) | Error vs `Δx`, multiple meshes |

## Common setup

All stages use:

- **Mesh**: `quad40.unv` (40 × 40 uniform Cartesian, 1000 × 1000 m)
- **Phases**: same density, same viscosity (`ρ = 1.0`, `μ = 1e-6`) — uniform
  fluid so the only physics is advection
- **No gravity, no surface tension** — pure transport
- **Initial condition**: planar α step at x = 500 m, advected by `U = (1, 0, 0) m/s`
- **Boundary conditions**: Dirichlet U at inlet, Zerogradient elsewhere; α
  zerogradient everywhere

Treating the test as quasi-1D: every horizontal row of cells evolves
identically, so the 1D modified equation applies cell-by-cell along x.

## Required edits to `Solvers_2_MULTIPHASE.jl`

None. The stages set `cAlpha` and other knobs through the standard `VOF(...)`
API.

## Running

From the package root:

```julia
julia --project=. test/vof_testing/stage1_upwind_diffusion.jl
julia --project=. test/vof_testing/stage2_compression_balance.jl
julia --project=. test/vof_testing/stage3_convergence.jl
```

Each stage prints a table with measured vs predicted values and a
pass/fail indicator. Stage 3 requires multiple mesh files (`quad20.unv`,
`quad40.unv`, `quad80.unv`); currently only `quad40` is included — extend
the `MESH_FILES` list when more meshes become available.

## Interpretation notes

- Stage 1 failing → either the upwind LO flux assembly is wrong or the
  initial condition is being set inconsistently. Check `phiLf` per face.
- Stage 1 passing but Stage 2 failing → MULES limiter or compression flux
  has a bug. Check `phirf` magnitude and the FCT bounds.
- Stages 1–2 passing but Stage 3 showing slope > 2 or < 1 → a bug elsewhere
  (boundary handling, time-discretisation order, or mass leakage at walls).
