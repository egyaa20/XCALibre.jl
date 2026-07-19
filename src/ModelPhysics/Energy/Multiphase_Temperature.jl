export MultiphaseTemperature, MultiphaseTemperatureModel

"""
    MultiphaseTemperature <: AbstractEnergyModel

Mixture-temperature energy model for two-phase flows. Solves

    ∂(ρ_m c_{p,m} T)/∂t + ∇·(ρ_m c_{p,m} U T) − ∇·(k_m ∇T) = ∂p/∂t + S_T

The `∂p/∂t` pressure-work source (enabled by `pressure_work=true`, the default)
is the low-Mach approximation of the exact term β_T·T·Dp/Dt with the
dimensionless group β_T·T taken as 1 — exact for an ideal gas (β_T = 1/T),
an over-estimate for liquids (β_T·T ≈ 0.45 for liquid N2 at 78 K). The
h-formulation carries the exact Dp/Dt with no such factor — prefer it where
pressure work matters. Disable with `pressure_work=false`.

where ρ_m, c_{p,m}, k_m are the α-blended mixture properties:

    c_{p,m}(α) = α·cp1 + (1−α)·cp2
    k_m(α)     = α·k1  + (1−α)·k2

α is the tracked-phase volume fraction (the existing `model.fluid.alpha`),
and `cp1, cp2, k1, k2` come from `Phase(rho=…, mu=…, cp=…, k=…)`.

### Fields
- `T`        : cell temperature [K]
- `Tf`       : face temperature [K]
- `cp_m`     : mixture cp at cells          [J/(kg·K)]
- `cpf_m`    : mixture cp at faces
- `k_m`      : mixture k  at cells          [W/(m·K)]
- `kf_m`     : mixture k  at faces
- `rho_cp`   : ρ·cp at cells (Time-term coefficient, ρcpⁿ⁺¹)
- `rho_cp_prev` : ρ·cp at the start of the time step (ρcpⁿ); RHS of the
                  conservative time term
- `mdot_cpf` : ρ·cp·U·n·A at faces (Divergence flux)
- `dpdt`     : pressure-work source ∂p/∂t [W/m³] (zeros when disabled)
- `p_old`    : pⁿ snapshot for the ∂p/∂t finite difference
- `coeffs`   : NamedTuple of constant coefficients (e.g. `T_init`)
"""
struct MultiphaseTemperature{S1,F1,S2,F2,S3,F3,S4,S6,F4,S5,S7,S8,C} <: AbstractEnergyModel
    T::S1
    Tf::F1
    cp_m::S2
    cpf_m::F2
    k_m::S3
    kf_m::F3
    rho_cp::S4
    rho_cp_prev::S6
    mdot_cpf::F4
    T_source::S5         # external T source (Lee latent heat etc.); zeros if unused
    dpdt::S7             # pressure-work source ∂p/∂t; zeros if pressure_work=false
    p_old::S8            # pⁿ snapshot (p = p_rgh + ρ·gh)
    coeffs::C
end
Adapt.@adapt_structure MultiphaseTemperature

struct MultiphaseTemperatureModel{E1, State}
    energy_eqn::E1
    state::State
end
Adapt.@adapt_structure MultiphaseTemperatureModel

# API constructor
#
# `Pr_t` is the turbulent Prandtl number, used to convert the eddy
# viscosity `ν_t` from the active RANS/LES turbulence model into a
# turbulent thermal diffusivity. Defaults to `0.85` (standard for
# gases and most liquids). Set to `Inf` to disable turbulent thermal
# diffusion entirely while keeping the turbulence model active for the
# momentum equation.
#
# `prop_relax ∈ (0, 1]` under-relaxes the time-term coefficient ρcp across
# outer (PIMPLE) correctors. `1.0` (default) takes the fresh table value
# each corrector (no relaxation). Lower it (e.g. `0.3–0.5`) for real-fluid
# runs where ρcp swings steeply through the pseudo-critical line and the
# conservative time term otherwise destabilises the outer loop. It is
# consistent — at outer convergence the relaxation vanishes.
Energy{MultiphaseTemperature}(; T_init=300.0, Pr_t=0.85, Pr_t_model::Symbol=:constant,
                                prop_relax=1.0, pressure_work::Bool=true) = begin
    Pr_t_model in (:constant, :kays_crawford) ||
        error("MultiphaseTemperature: Pr_t_model must be :constant or :kays_crawford, got $Pr_t_model")
    (0.0 < prop_relax <= 1.0) ||
        error("MultiphaseTemperature: prop_relax must be in (0, 1], got $prop_relax")
    coeffs = (T_init=T_init, Pr_t=Pr_t, Pr_t_model=Pr_t_model,
              prop_relax=float(prop_relax), pressure_work=pressure_work)
    ARG = typeof(coeffs)
    Energy{MultiphaseTemperature, ARG}(coeffs)
end

# Functor that allocates fields once the mesh is known.
(energy::Energy{EnergyModel, ARG})(mesh, fluid) where {EnergyModel<:MultiphaseTemperature, ARG} = begin
    T        = ScalarField(mesh)
    Tf       = FaceScalarField(mesh)
    cp_m     = ScalarField(mesh)
    cpf_m    = FaceScalarField(mesh)
    k_m      = ScalarField(mesh)
    kf_m     = FaceScalarField(mesh)
    rho_cp   = ScalarField(mesh)
    rho_cp_prev = ScalarField(mesh)
    mdot_cpf = FaceScalarField(mesh)
    T_source = ScalarField(mesh)   # zeros — populated by Lee/RPI updates each iter
    dpdt     = ScalarField(mesh)   # zeros — populated per step if pressure_work
    p_old    = ScalarField(mesh)
    coeffs   = energy.args
    MultiphaseTemperature(T, Tf, cp_m, cpf_m, k_m, kf_m, rho_cp, rho_cp_prev,
                          mdot_cpf, T_source, dpdt, p_old, coeffs)
end

"""
    initialise(energy::MultiphaseTemperature, model, mdotf, rho, peqn, config)

Builds the temperature transport equation and pre-allocates its solver
workspace. `peqn` is reused only for `peqn.equation` (provides the
ScalarEquation skeleton sized to the mesh).
"""
function initialise(
    energy::MultiphaseTemperature, model::Physics{T1,F,SO,M,Tu,E,D,BI},
    mdotf, rho, peqn, config
) where {T1,F,SO,M,Tu,E,D,BI}

    (; T, rho_cp, kf_m, mdot_cpf, T_source, dpdt, p_old) = energy
    (; solvers, schemes, boundaries) = config
    eqn = peqn.equation

    # `+ Source(T_source)` is for Lee phase-change latent heat (and any
    # other external source). The field starts at zero — only populated
    # if Lee/RPI is configured. With it zero, the equation is unchanged.
    # `+ Source(dpdt)` is the pressure-work term ∂p/∂t (zeros when disabled).
    energy_eqn = (
        Time{schemes.T.time}(rho_cp, T)
        + Divergence{schemes.T.divergence}(mdot_cpf, T)
        - Laplacian{schemes.T.laplacian}(kf_m, T)
        ==
        Source(ConstantScalar(0.0))
        + Source(T_source)
        + Source(dpdt)
    ) → eqn

    # Seed the pⁿ snapshot with the full pressure the solver maintains
    # (p = p_rgh + ρ·gh). Without this, the first step's ∂p/∂t would be the
    # spurious (p − 0)/Δt of the whole hydrostatic field.
    gh = model.fluid.physics_properties.gravity.gh
    @. p_old.values = model.fluid.p_rgh.values + rho.values * gh.values

    @reset energy_eqn.preconditioner = set_preconditioner(solvers.T.preconditioner, energy_eqn)
    @reset energy_eqn.solver         = _workspace(solvers.T.solver, _b(energy_eqn))

    init_residual = (:T, 1.0)
    init_converged = false
    state = ModelState(init_residual, init_converged)

    return MultiphaseTemperatureModel(energy_eqn, state)
end

"""
    energy!(energy::MultiphaseTemperatureModel, model, prev, mdotf, rho, time, config)

Per-iteration update:
  1. Blend cp and k from the two phases using α (cells and faces).
  2. Form the time-term coefficient ρ·cp and the divergence face flux ρ·cp·U·n·A.
  3. Discretise / solve the linear system for T.
  4. Refresh the face temperature.

Boundaries for T are taken from `config.boundaries.T`. The solver / schemes
must therefore include a `T` entry (e.g. `Schemes(time=Euler, divergence=Upwind, laplacian=Linear)`).
"""
function energy!(
    energy::MultiphaseTemperatureModel, model::Physics{T1,F,SO,M,Tu,E,D,BI},
    prev, mdotf, rho, time, config; outer=1
) where {T1,F,SO,M,Tu,E,D,BI}

    mesh = model.domain
    (; T, Tf, cp_m, cpf_m, k_m, kf_m, rho_cp, rho_cp_prev, mdot_cpf, dpdt, p_old, coeffs) = model.energy
    (; energy_eqn, state) = energy
    (; solvers, hardware, boundaries) = config
    (; alpha, alphaf) = model.fluid

    phases = model.fluid.phases
    Pr_t = coeffs.Pr_t
    # Kays-Crawford toggle: when true, kernel computes a local Pr_t from
    # ν_t/ν per cell using the Kays-Crawford correlation. Passed as Int
    # for GPU-bitstype safety (Bool isa bitstype too but Int is universal).
    use_kc = (coeffs.Pr_t_model === :kays_crawford) ? 1 : 0

    # 1) Mixture cp and k_eff at cells.
    #    k_eff = k_m + ρ·cp_m·ν_t/Pr_t
    #
    #    The molecular `k_m` is α-blended from per-phase values, which
    #    may be `ConstantScalar` (uniform) or `ScalarField` (per-cell,
    #    from a HelmholtzTable live update). On top of that we add the
    #    turbulent thermal diffusion `k_t = ρ·cp·ν_t/Pr_t` whenever the
    #    active turbulence model carries a `nut` field (KOmega, KOmegaSST,
    #    Smagorinsky, etc.). For `Laminar` there is no `nut` and the
    #    contribution is zero.
    nut_values = _turbulence_nut_values(model.turbulence)
    # Pass a scalar 0 for Laminar (no nut field) so the kernel can read a
    # uniform value; otherwise the per-cell nut array.
    nut_arg = isnothing(nut_values) ? zero(eltype(cp_m.values)) : nut_values

    # cp cap for the turbulent thermal conductivity (Reynolds analogy
    # breaks down through pseudo-critical peaks). Inf disables the cap.
    cp_kt_cap = convert(eltype(cp_m.values), Inf)

    (; backend, workgroup) = hardware
    ndrange = length(alpha.values)
    kernel! = _blend_cp_k!(_setup(backend, workgroup, ndrange)...)
    kernel!(
        cp_m.values, k_m.values, alpha.values,
        phases[1].cp.values, phases[2].cp.values,
        phases[1].k.values,  phases[2].k.values,
        rho.values, nut_arg, model.fluid.nu.values,
        Pr_t, cp_kt_cap, use_kc,
    )

    # 2) Mixture cp and k_eff at faces. Use the cell-level resolved
    #    mixture values and `interpolate!` to get face values.
    interpolate!(cpf_m, cp_m, config)
    interpolate!(kf_m,  k_m,  config)

    # 3) Divergence face flux (always the fresh value).
    @. mdot_cpf.values = mdotf.values * cpf_m.values

    # Time-term coefficient ρcp, with optional property under-relaxation
    # across outer correctors.
    #
    # `outer == 1` establishes ρcpⁿ: properties are still evaluated at Tⁿ,
    # so the fresh `rho·cp_m` IS ρcpⁿ. We take it un-relaxed and snapshot it
    # into `rho_cp_prev` as the frozen RHS coefficient for the conservative
    # time term (this also makes the very first time step correct, avoiding
    # a zero-initialised snapshot).
    #
    # Later correctors advance `rho_cp` toward ρcpⁿ⁺¹. Real-fluid ρcp swings
    # steeply through the pseudo-critical line and the conservative source
    # ∝ (ρcpⁿ⁺¹ − ρcpⁿ)/Δt amplifies any corrector-to-corrector jump (∝ 1/Δt).
    # Blending the fresh table value with the previous corrector's `rho_cp`
    # (factor `prop_relax`) damps that swing so the outer loop contracts —
    # the standard real-fluid remedy (cf. OpenFOAM `rho.relax()`). It is
    # consistent: at outer convergence the fresh and previous values
    # coincide, so the relaxation vanishes and the true ρcpⁿ⁺¹ is recovered.
    # `prop_relax == 1` reproduces the un-relaxed behaviour exactly.
    β = coeffs.prop_relax
    if outer == 1
        @. rho_cp.values      = rho.values * cp_m.values
        @. rho_cp_prev.values = rho_cp.values
    else
        @. rho_cp.values = (1 - β) * rho_cp.values + β * (rho.values * cp_m.values)
    end

    # 3b) Pressure-work source ∂p/∂t — approximates the exact β_T·T·Dp/Dt with
    #     β_T·T ≈ 1 (ideal gas; ~0.45 for liquid N2) and drops the convective
    #     U·∇p (low-Mach, consistent with neglecting kinetic-energy transport;
    #     the h-formulation carries the exact term). `p = p_rgh + ρ·gh` is refreshed by the solver before
    #     the energy step. Computed once per time step (outer == 1) against
    #     the pⁿ snapshot — same convention as the ρcp snapshot above; later
    #     outer correctors reuse it.
    if coeffs.pressure_work
        p = model.momentum.p
        if outer == 1
            dt_cpu = zeros(eltype(dpdt.values), 1)
            copyto!(dt_cpu, config.runtime.dt)
            @. dpdt.values = (p.values - p_old.values) / dt_cpu[1]
            @. p_old.values = p.values
        end
    end

    # 4) Solve the T equation (conservative time term via `rho_prev`).
    @. prev = T.values
    discretise!(energy_eqn, T, config; rho_prev=rho_cp_prev)
    apply_boundary_conditions!(energy_eqn, boundaries.T, nothing, time, config)
    implicit_relaxation_diagdom!(energy_eqn, T.values, solvers.T.relax, nothing, config)
    update_preconditioner!(energy_eqn.preconditioner, mesh, config)
    T_res = solve_system!(energy_eqn, solvers.T, T, nothing, config)

    # 5) Refresh face temperature for downstream interpolations.
    interpolate!(Tf, T, config)
    correct_boundaries!(Tf, T, boundaries.T, time, config)

    state.residuals = (:T, T_res)
    state.converged = T_res <= solvers.T.convergence
    return nothing
end

# Helpers — resolve a per-phase property handle to a cell value.
# Phase properties are either `ConstantScalar` (uniform) or `ScalarField`
# (per-cell, populated by HelmholtzTable's live update routine).
@inline _phase_cell_value(p::ConstantScalar, _i)              = p.values
@inline _phase_cell_value(p::ScalarField,    i::Integer)      = p.values[i]

# In-kernel uniform accessor: a property `.values` is either a scalar
# (ConstantScalar) or a per-cell array (ScalarField). Dispatch resolves
# at compile time, so this stays GPU-safe.
@inline _at(x::Number, _i)        = x
@inline _at(x::AbstractArray, i)  = @inbounds x[i]

# Kays-Crawford (1994) variable Pr_t correlation.
#   1/Pr_t = 0.5882 + 0.228·x − 0.0441·x²·[1 − exp(−5.165/x)]
# where x = ν_t / ν is the local turbulent-viscosity ratio.
# Behaviour: x → 0 (laminar) → Pr_t → 1/0.5882 ≈ 1.7
#            x → ∞ (fully turbulent) → Pr_t → ~0.85
# Numerically guarded against x = 0.
@inline function _prt_kays_crawford(x::T) where {T<:Real}
    xs   = max(x, eps(T))
    inv_Prt = T(0.5882) + T(0.228)*xs - T(0.0441)*xs*xs * (one(T) - exp(-T(5.165)/xs))
    return one(T) / max(inv_Prt, eps(T))
end

# Blend mixture cp and effective k (molecular + turbulent) per cell:
#   k_eff = k_m + ρ·min(cp, cp_kt_cap)·ν_t/Pr_t
# `cp_kt_cap` caps the cp value that feeds the turbulent thermal
# conductivity. The storage-term cp (used to assemble rho_cp afterwards)
# stays at the full value. Pass `cp_kt_cap = Inf` to disable the cap.
# `use_kc` selects the Pr_t model:
#   0 → constant scalar `Pr_t_const`
#   1 → Kays-Crawford local Pr_t(ν_t/ν)
@kernel inbounds=true function _blend_cp_k!(
        cp_m, k_m, alpha, cp1, cp2, k1, k2, rho, nut, nu,
        Pr_t_const, cp_kt_cap, use_kc)
    i = @index(Global)
    a      = alpha[i]
    cp_m_i = _at(cp1, i) * a + _at(cp2, i) * (1.0 - a)
    k_m_i  = _at(k1, i)  * a + _at(k2, i)  * (1.0 - a)
    cp_kt  = min(cp_m_i, cp_kt_cap)

    nut_i  = _at(nut, i)
    if use_kc == 1
        nu_i   = _at(nu, i)
        x      = nut_i / max(nu_i, eps(typeof(nu_i)))
        Pr_t_l = _prt_kays_crawford(x)
    else
        Pr_t_l = Pr_t_const
    end

    k_t_i  = rho[i] * cp_kt * nut_i / Pr_t_l
    cp_m[i] = cp_m_i
    k_m[i]  = k_m_i + k_t_i
end

# Turbulence `nut` accessor — returns the per-cell `nut.values` array
# when the active turbulence model carries an eddy-viscosity field, or
# `nothing` for `Laminar`. The energy! loop interprets `nothing` as
# "no turbulent thermal diffusion".
@inline _turbulence_nut_values(turb) =
    hasproperty(turb, :nut) ? getproperty(turb, :nut).values : nothing
