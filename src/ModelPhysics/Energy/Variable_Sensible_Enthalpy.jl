export VariableSensibleEnthalpy, VariableSensibleEnthalpyModel

"""
    VariableSensibleEnthalpy <: AbstractEnergyModel

Mixture **enthalpy** energy model for real-fluid / supercritical flows.
Solves the sensible-enthalpy transport equation

    ∂(ρ h)/∂t + ∇·(ṁ h) − ∇·((k_eff/cp) ∇h) = q_src

with all properties (ρ, cp, k) taken live from a `HelmholtzTable`, blended
across the two phases exactly like [`MultiphaseTemperature`](@ref). The
temperature is recovered each iteration by the **monotone inversion**
`T = T(h)` of the tabulated `h(T)` grid.

### Why enthalpy instead of temperature
The temperature equation advances T by dividing the deposited energy by the
heat capacity, `ΔT ≈ q·dt/(ρ cp V)`. Through the pseudo-critical line cp
spikes then falls, and a cell that overshoots the peak lands on the
falling side where lower cp → larger ΔT → runaway to the table ceiling.
Here the unknown is `h` (energy itself): the storage term carries no cp in
a denominator, `h(T)` is smooth and strictly monotone (the cp spike is just
a steep ramp), and `T(h)` is a bounded, single-valued inverse that
*plateaus* through pseudo-critical instead of running away. The cp spike
only reappears as the benign diffusion coefficient `k_eff/cp`.

### Fields
- `h`, `hf`     : sensible enthalpy at cells / faces [J/kg]
- `T`, `Tf`     : temperature recovered from `h` at cells / faces [K]
- `cp_m`        : mixture cp [J/(kg·K)] (sets the diffusion coefficient)
- `k_m`         : mixture effective conductivity k + k_t [W/(m·K)]
- `alpha_eff`   : `k_eff/cp` at cells (Laplacian coefficient)
- `alpha_eff_f` : `k_eff/cp` at faces (Laplacian flux)
- `rho_h`       : Time-term coefficient ρ (ρⁿ⁺¹; optionally under-relaxed)
- `rho_h_prev`  : ρⁿ snapshot — RHS of the conservative time term
- `T_source`    : external heat source [W/m³] (Lee/RPI); zeros if unused
- `coeffs`      : NamedTuple (`T_init`, `Pr_t`, `Pr_t_model`, `prop_relax`)
"""
struct VariableSensibleEnthalpy{S,F,C} <: AbstractEnergyModel
    h::S
    hf::F
    T::S
    Tf::F
    cp_m::S
    k_m::S
    alpha_eff::S
    alpha_eff_f::F
    rho_h::S
    rho_h_prev::S
    T_source::S
    coeffs::C
end
Adapt.@adapt_structure VariableSensibleEnthalpy

struct VariableSensibleEnthalpyModel{E1,State}
    energy_eqn::E1
    state::State
end
Adapt.@adapt_structure VariableSensibleEnthalpyModel

# API constructor — mirrors MultiphaseTemperature (same knobs), so a case
# can switch formulation by changing only the model type.
Energy{VariableSensibleEnthalpy}(; T_init=300.0, Pr_t=0.85, Pr_t_model::Symbol=:constant, prop_relax=1.0) = begin
    Pr_t_model in (:constant, :kays_crawford) ||
        error("VariableSensibleEnthalpy: Pr_t_model must be :constant or :kays_crawford, got $Pr_t_model")
    (0.0 < prop_relax <= 1.0) ||
        error("VariableSensibleEnthalpy: prop_relax must be in (0, 1], got $prop_relax")
    coeffs = (T_init=T_init, Pr_t=Pr_t, Pr_t_model=Pr_t_model, prop_relax=float(prop_relax))
    ARG = typeof(coeffs)
    Energy{VariableSensibleEnthalpy, ARG}(coeffs)
end

# Functor — allocate fields once the mesh is known.
(energy::Energy{EnergyModel, ARG})(mesh, fluid) where {EnergyModel<:VariableSensibleEnthalpy, ARG} = begin
    h           = ScalarField(mesh)
    hf          = FaceScalarField(mesh)
    T           = ScalarField(mesh)
    Tf          = FaceScalarField(mesh)
    cp_m        = ScalarField(mesh)
    k_m         = ScalarField(mesh)
    alpha_eff   = ScalarField(mesh)
    alpha_eff_f = FaceScalarField(mesh)
    rho_h       = ScalarField(mesh)
    rho_h_prev  = ScalarField(mesh)
    T_source    = ScalarField(mesh)
    coeffs      = energy.args
    VariableSensibleEnthalpy(h, hf, T, Tf, cp_m, k_m, alpha_eff, alpha_eff_f,
                             rho_h, rho_h_prev, T_source, coeffs)
end

# Fetch the HelmholtzTable the model depends on (errors clearly if absent —
# this model is meaningless without a real-fluid table).
_helmholtz_table(model) = begin
    fp = get(model.fluid.physics_properties, :fluid_properties, nothing)
    fp isa HelmholtzTable ||
        error("VariableSensibleEnthalpy requires `fluid_properties = HelmholtzTable(...)` on the Fluid.")
    fp
end

"""
    initialise(energy::VariableSensibleEnthalpy, model, mdotf, rho, peqn, config)

Build the enthalpy transport equation and seed `h` from the current `T`
field (so `h` starts consistent with the case's initial temperature).
"""
function initialise(
    energy::VariableSensibleEnthalpy, model::Physics{T1,F,SO,M,Tu,E,D,BI},
    mdotf, rho, peqn, config
) where {T1,F,SO,M,Tu,E,D,BI}

    (; h, T, rho_h, alpha_eff_f, T_source) = energy
    (; solvers, schemes, boundaries) = config
    eqn = peqn.equation

    energy_eqn = (
        Time{schemes.h.time}(rho_h, h)
        + Divergence{schemes.h.divergence}(mdotf, h)
        - Laplacian{schemes.h.laplacian}(alpha_eff_f, h)
        ==
        Source(ConstantScalar(0.0))
        + Source(T_source)
    ) → eqn

    @reset energy_eqn.preconditioner = set_preconditioner(solvers.h.preconditioner, energy_eqn)
    @reset energy_eqn.solver         = _workspace(solvers.h.solver, _b(energy_eqn))

    # Seed h from T (case sets T = T_init before the run). ρ-coefficient is
    # initialised to the current density so the first assembly is sane.
    table = _helmholtz_table(model)
    interp_h_from_T!(h, T, table, config)
    @. rho_h.values = rho.values

    init_residual = (:h, 1.0)
    state = ModelState(init_residual, false)
    return VariableSensibleEnthalpyModel(energy_eqn, state)
end

"""
    energy!(energy::VariableSensibleEnthalpyModel, model, prev, mdotf, rho, time, config; outer=1)

Per-corrector update: blend cp/k_eff, form the `k_eff/cp` diffusion
coefficient, assemble and solve the enthalpy equation (conservative time
term with optional ρ under-relaxation), then recover `T = T(h)` by the
monotone table inversion and refresh face fields.
"""
function energy!(
    energy::VariableSensibleEnthalpyModel, model::Physics{T1,F,SO,M,Tu,E,D,BI},
    prev, mdotf, rho, time, config; outer=1
) where {T1,F,SO,M,Tu,E,D,BI}

    mesh = model.domain
    (; h, hf, T, Tf, cp_m, k_m, alpha_eff, alpha_eff_f, rho_h, rho_h_prev, coeffs) = model.energy
    (; energy_eqn, state) = energy
    (; solvers, hardware, boundaries) = config
    (; alpha) = model.fluid

    phases = model.fluid.phases
    table  = _helmholtz_table(model)
    Pr_t   = coeffs.Pr_t
    use_kc = (coeffs.Pr_t_model === :kays_crawford) ? 1 : 0
    β      = coeffs.prop_relax

    # 1) Blend mixture cp and effective conductivity k_eff = k_m + ρ·cp·ν_t/Pr_t
    #    (same kernel as MultiphaseTemperature). No cp cap here — the cp peak
    #    is harmless in enthalpy form.
    nut_values = _turbulence_nut_values(model.turbulence)
    nut_arg    = isnothing(nut_values) ? zero(eltype(cp_m.values)) : nut_values
    cp_kt_cap  = convert(eltype(cp_m.values), Inf)

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

    # 2) Effective thermal diffusivity for the enthalpy Laplacian: k_eff/cp.
    #    (Dips at pseudo-critical — a small positive coefficient, never a
    #    division of the unknown, so harmless.)
    @. alpha_eff.values = k_m.values / cp_m.values
    interpolate!(alpha_eff_f, alpha_eff, config)

    # 3) Time-term coefficient ρ for the conservative ∂(ρh)/∂t, with optional
    #    under-relaxation across outer correctors (cf. MultiphaseTemperature).
    #    `outer == 1` establishes ρⁿ (properties still at Tⁿ) and snapshots it.
    if outer == 1
        @. rho_h.values      = rho.values
        @. rho_h_prev.values = rho_h.values
    else
        @. rho_h.values = (1 - β) * rho_h.values + β * rho.values
    end

    # 4) Assemble and solve the enthalpy equation.
    @. prev = h.values
    discretise!(energy_eqn, h, config; rho_prev=rho_h_prev)
    apply_boundary_conditions!(energy_eqn, boundaries.h, nothing, time, config)
    implicit_relaxation_diagdom!(energy_eqn, h.values, solvers.h.relax, nothing, config)
    update_preconditioner!(energy_eqn.preconditioner, mesh, config)
    h_res = solve_system!(energy_eqn, solvers.h, h, nothing, config)

    # 5) Recover T from h by the monotone inversion (runaway-free), refresh faces.
    invert_T_from_h!(T, h, table, config)
    interpolate!(hf, h, config)
    correct_boundaries!(hf, h, boundaries.h, time, config)
    interpolate!(Tf, T, config)
    correct_boundaries!(Tf, T, boundaries.T, time, config)

    state.residuals = (:h, h_res)
    state.converged = h_res <= solvers.h.convergence
    return nothing
end
