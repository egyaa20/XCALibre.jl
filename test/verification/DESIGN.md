# XCALibre Multiphase V&V Platform — Design (v2)

A tiered verification & validation suite for the multiphase solver, rerunnable
after every solver change. v2 reworks the weak tests of v1 around four rules
learned from this project's own history:

1. **Test invariants and fixed points, not trajectories.** Anything asserted
   "after N iterations at dt=X" is a lottery ticket. Assert properties that are
   independent of dt, duration, and mesh by construction — or assert the
   *independence itself*.
2. **No absolute magic tolerances.** Every tolerance is tied to round-off, the
   linear-solver tolerance, a convergence *rate*, or a dimensionless physical
   ratio. A tolerance you can't justify is a tolerance that will be tuned until
   the test passes.
3. **Two classes of test, clearly separated.**
   - **Class A — physics contracts** (conservation, thermodynamic consistency,
     reduction limits, discrete equilibria, analytic solutions). Written against
     the solver's *interface* (fields in → fields out). Survive any
     reformulation of the internals.
   - **Class B — implementation-pinned unit tests** (e.g. `compute_psi` vs
     finite differences). Cheap armor while that code exists; **deleted with the
     module, by design**. The platform's value must not depend on Class B.
     *Case in point: the ψ-split pressure equation was replaced by a lumped
     explicit ∂ρ/∂t within days. Under this split, only one Class-B file dies;
     every Class-A test below runs unchanged and still verifies the new form.*
4. **Validation ≠ CI.** Experiment-comparison campaigns (Tatsumoto) are manual,
   documented studies. The automated suite tops out at correlation checks and
   robustness pass/fails that run in minutes.

---

## Traceability: bugs this platform would have caught

| Past bug | Test |
|---|---|
| cp 1000× molar/mass unit error | 5.1 global energy balance (off by 1000×) |
| Missing ρ in energy convection flux | 5.1/5.2 energy balance (off by ~ρ) |
| ∇·U=0 continuity with ρ(T) → T checkerboard | 5.4 heated-duct mass-flux constancy |
| maxCo=10 vs 90 giving different physics | 3.4 / 5.7 Courant-invariance |
| NaN → `unsafe_trunc` → illegal GPU address | 0.3 data-contract NaN safety |
| β/entropy swapped in EOS functor returns | 0.2 table self-consistency (cp≠dh/dT) |
| Spurious currents at stratified interface | 3.1 discrete-equilibrium fixed point |

---

## Tier 0 — Property **data contracts** (pure functions, < 10 s, CPU)

File: `tier0_contracts.jl` — Class A. These test whatever property source the
solver consumes, through its public lookup interface. They do not care how the
table is built or which PDE terms use which column.

- **0.1 Interpolation fidelity.** Table lookup at off-grid (T,p) vs direct EOS
  call at the same state: ρ, cp, μ, k rel. error < 1e-4 away from T_c ± 2 K
  (< 1e-2 inside the spike, where linear interpolation error is the point).
- **0.2 Thermodynamic self-consistency (formulation-proof).** Consistency
  *between* exposed quantities, by finite differences **on the table itself**:
  cp ≈ ∂h/∂T|_p (1%); h strictly monotone in T per pressure column; T↔h
  round-trip < 1e-9 K; and *if* the source exposes ψ or ∂ρ/∂T, they must equal
  FD of its own ρ(T,p) (auto-skipped if the columns don't exist). If someone
  swaps return order or units anywhere in the chain, cross-derivative identities
  break loudly.
- **0.3 Poisoned-input safety.** Lookups and inversions with T,p ∈
  {NaN, ±Inf, −1e300}: finite edge-clamped outputs, no throw, no out-of-range
  index (guards the illegal-GPU-address class of failure).
- **0.4 Grid contract.** Strictly increasing T grid, no duplicates; refinement
  band densities match the constructor's promise (×5 within ±1 K, ×3 within
  ±10 K) — asserted as *ratios of local spacing*, not point counts.
- **0.5 Blend identities.** α∈{0,½,1} endpoints/midpoint for `blend_properties!`;
  ν·ρ_blend ≡ linear μ-blend.
- **0.6 Boiling closure identities.** Lee: zero on the wrong side of T_sat, hand
  value at T_sat+1 K; RPI: q_c+q_q+q_e ≡ q_w partition identity at a reference
  state; submodel formulas vs hand-evaluated numbers.

File: `unitB_helmholtz.jl` — **Class B, disposable.** `compute_psi` /
`compute_drho_dT_p` vs FD of the EOS, refined-grid internals, kernel index
helpers. Rewritten or deleted whenever the Helmholtz module changes shape.

## Tier 1 — Discrete operator verification (raw kernels, < 30 s)

File: `tier1_operators.jl` (pattern of existing `unit_test_laplace.jl`)

- **1.1** Gauss gradient exact for linear φ (machine precision, interior cells).
- **1.2** Laplacian of x²+y²: observed order ≥ 1.8 on quad40/80/160.
- **1.3** `div!` of a uniform prescribed face flux = 0 to round-off.
- **1.4** cell→face interpolation 2nd order on the quad ladder.

## Tier 2 — Single-equation verification: MMS & exact solutions (~1–2 min)

File: `tier2_equations.jl`

- **2.1 Energy MMS.** `T_source` is a user-injectable volumetric source → full
  Method of Manufactured Solutions: frozen uniform U, constant props,
  T* = T₀ + ΔT·sin(πx/L)cos(πy/L), inject S = ρcp U·∇T* − k∇²T*, steady solve
  on quad40/80/160, L2(T−T*) order: ≥ 1.8 diffusion-dominated, ≈ 1
  convection-dominated (Upwind — assert what the scheme *is*). Run with both
  `MultiphaseTemperature` and `VariableSensibleEnthalpy` (constant-cp table):
  same target ⇒ also proves T/h equivalence.
- **2.2 Transient conduction.** U=0 slab, wall step, erfc solution while the
  front is deep inside the domain; dt-halving ⇒ Euler temporal order ≈ 1.
  (Order asserted, not absolute error — mesh/dt-robust by construction.)
- **2.3 α step advection (1D, frozen flux).** Boundedness α∈[−1e-12, 1+1e-12]
  every step; conservation to 1e-12 rel; front within 1 cell of x=Ut; MULES
  interface width < pure-upwind width.
- **2.4 Zalesak-lite.** One solid-body revolution of a circle on quad100:
  bounded, conservative, L1 shape error bounded.

## Tier 3 — Pressure–velocity coupling, incompressible baseline (~2 min)

File: `tier3_coupling.jl`

- **3.1 Hydrostatic equilibrium as a discrete fixed point** (replaces the
  old duration/dt-sensitive spurious-currents test).
  Initialise the *exactly representable* equilibrium: interface aligned to a
  face row (box init on grid lines), U=0, p_rgh consistent. Then:
  - (a) **One-step spurious acceleration.** After a single step,
    a_sp = max|U|/dt is the residual force per unit mass and is
    **dt-independent to first order**. Assert the dimensionless
    a_sp/g < tol, with tol derived from the p_rgh linear-solve tolerance —
    and assert **dt-invariance itself**: a_sp at dt=1e-5 and dt=1e-4 agree
    within 2×. No iteration-count lottery.
  - (b) **No secular growth.** 200 steps at both dts:
    max|U|(end) ≤ 3 × max|U|(after 10 steps). Bounds the trend, not a magic
    absolute at a magic time.
- **3.2 Poiseuille channel.** Body-force-driven laminar channel to steady:
  parabolic profile L2 < 1% on quad80; identical mass flux through every
  cross-section (round-off).
- **3.3 Free-surface dynamics with analytic targets** (replaces the
  timed-once dam-break check):
  - (a) **Standing-wave (sloshing) frequency — primary.** Small staircase-tilted
    interface in a box; mode-1 gravity wave has ω² = gk·tanh(kh), k=π/L.
    Measure the period from zero-crossings of the wall water-column height over
    ≥3 periods; assert within 5%. Frequency is insensitive to run length, dt
    (CFL-stable), and interface smearing — exactly the robustness the old test
    lacked. Exact-conservation assert rides along for free.
  - (b) **Column collapse vs Martin–Moyce (1952) — secondary.** Track the
    dimensionless front x/a vs t√(2g/a) against the tabulated experimental
    curve (±15% band) until the front nears the far wall. Turns the dam break
    from "nothing went wrong" into a real validation.
  - (c) **RTI growth rate — stretch.** Grids exist (RTI_VOF_*): seed a
    single-mode perturbation, fit early-time exponential amplitude growth,
    compare to σ = √(A g k) (inviscid band). Optional; most expensive.
- **3.4 Courant invariance (regression).** Short transient at maxCo 0.5 vs 5
  with outer correctors: final fields agree < 0.5%.

## Tier 4 — Mixture (drift-flux) validation (~2 min)

File: `tier4_mixture.jl`

- **4.1 Stokes terminal velocity.** Quiescent dilute column, Re_p < 0.1:
  computed drift = Δρ g d²/(18 μ_c) within 2%. Validates `compute_Ur!` + drag.
- **4.2 Kynch sedimentation.** 1D settling column: clear-fluid interface
  descends at the analytic kinematic speed; front within 2 cells over the run.
- **4.3 Dispersed-phase conservation under drift** to round-off.
- **4.4 VOF limit.** d→0 ⇒ τ_d→0 ⇒ Ur→0: Mixture fields match a VOF run
  of the same α transport < 1e-10.

## Tier 5 — Energy + variable density: **formulation-independent contracts** (~3–5 min)

File: `tier5_energy_density.jl` — all Class A. These are the primary
verification of the compressible pressure equation, *whatever its current form*
(lumped explicit ∂ρ/∂t today; ψ-split, implicit, or anything else tomorrow):

- **5.1 Global energy balance, constant props.** Heated streamtube to steady:
  ṁ·cp·(T_out − T_in) = q·A within 1%.
- **5.2 Energy balance, h-formulation + real table.** ṁ(h_out − h_in) = q·A
  within 1%.
- **5.3 Constant-density reduction.** Variable-density path with a constant-ρ
  table ≡ incompressible path to near round-off. Any continuity formulation
  must degenerate exactly; this pins it.
- **5.4 Heated-duct continuity.** Strongly heated duct, N2 table, steady:
  ρ̄ŪA constant through every cross-section (< 1%); U_out/U_in ≈ ρ_in/ρ_out.
  The discrete statement of ∂ρ/∂t + ∇·(ρU) = 0 — the checkerboard root cause.
- **5.5 Transient mass audit.** Per step: d(ΣρV)/dt = ṁ_in − ṁ_out; bounded
  accumulated drift. Holds for *any* correct density-rate treatment.
- **5.6 T vs h equivalence** on the real table below pseudo-critical.
- **5.7 Courant invariance with energy on.** maxCo 5 vs 50 short heated run:
  sampled q–ΔT agree < 5%.
- **5.8 prop_relax consistency.** 0.4 vs 1.0 to steady: same steady fields
  < 0.1% (relaxation must vanish at convergence).

## Tier 6 — Boiling submodels in situ (~1 min)

File: `tier6_boiling.jl`

- **6.1** Lee in-solver: superheated box source = hand value; subcooled ≡ 0.
- **6.2** RPI in-solver: partition sums to applied q_w; T_source integral = −q_e·A.
- **6.3** Single-phase limit: below-onset q ⇒ q_e ≈ 0 and wall HTC within ±20%
  of Dittus–Boelter.

## Tier 7 — Feasible validation & robustness (minutes, opt-in; Tatsumoto goes manual)

File: `tier7_validation.jl`

The full Tatsumoto campaign (fine mesh, ramped q, hours of GPU, checkpoints) is
**not an automated test** — it is a documented manual study (`test/tatsumoto/`
runbook: run.ps1 stages, post.py, steady_state_analysis.py). The automated tier
keeps only what runs in minutes and has a defensible target:

- **7.1 Mini heated pipe vs Dittus–Boelter.** Coarse quarter-pipe, *constant*
  modest q (no ramp, no pseudo-critical crossing), short run to
  quasi-steady: wall HTC within ±25% of Nu = 0.023 Re^0.8 Pr^0.4. Exercises the
  full stack (turbulence + energy + table properties + wall flux) against an
  accepted correlation with an honest band — no experiment digitisation, no
  hours of runtime.
- **7.2 Pseudo-critical robustness (pass/fail, not accuracy).** Coarse mesh,
  aggressive constant q chosen to drive the hottest cells through T_pc quickly:
  assert the run completes with (i) no NaN/divergence-trap trip, (ii) T bounded
  by the table range, (iii) the near-wall std(T) speckle metric not growing
  monotonically (no checkerboard). This is the actual failure mode being
  fought, distilled into a minutes-long gate.
- **7.3 Golden-file drift detection.** Store 7.1/7.2 key outputs (JLD2);
  compare with tolerance on rerun. Catches silent behavior change without
  claiming experimental accuracy.

---

## Infrastructure

```
test/verification/
├── DESIGN.md
├── runtests.jl        tier selection: ENV["XCAL_VV_TIERS"]="0,1,2" or ARGS; CI default 0–3
├── common.jl          norms, observed_order, run_to_steady!, boundary audits,
│                      @vvtest macro (prints value/target/tol on failure)
├── tier0_contracts.jl … tier7_validation.jl        (Class A)
├── unitB_helmholtz.jl, unitB_<module>.jl …          (Class B, disposable)
└── golden/            JLD2 references for tier 7
```

- Boundary audits (`mass_flux(patch)`, `enthalpy_flux(patch)`, `wall_heat(patch)`)
  power every conservation contract — implement once, reuse everywhere.
- Hook into `test/runtests.jl` as `@testset "V&V"`; tiers 0–2 CI default;
  `XCAL_VV_BACKEND=CUDA` reruns tiers 3 & 5 on GPU (catches scalar indexing /
  `unsafe_trunc` classes).
- Tier 2.3/2.4 want a frozen-velocity α-only driver (prescribe `mdotf`, skip
  U/p) — a small solver hook worth adding for isolation and speed.

## Implementation order
1. Skeleton + `common.jl` audits + **Tier 0 contracts** (an afternoon).
2. **5.1 + 5.2 energy balances** (few lines; highest bug-catch rate).
3. **3.1 fixed-point hydrostatic** (replaces the flaky test immediately).
4. **Tier 2** MMS + α step.
5. **4.1 Stokes + 4.2 Kynch**; **3.3a sloshing frequency**.
6. Remaining Tier 5, then 6, then 7.
