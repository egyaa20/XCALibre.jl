export Lee, apply_lee_sources!

"""
    Lee(; C_e, C_c)

Lee (1979) volumetric phase-change model.

For each cell, an empirical mass-transfer rate is computed:

    ṁ_lv  =  C_e · α_l · ρ_l · max((T − T_sat) / T_sat, 0)        [kg/(m³·s)]
    ṁ_vl  =  C_c · α_v · ρ_v · max((T_sat − T) / T_sat, 0)        [kg/(m³·s)]

`ṁ_lv` (liquid → vapour) acts in *superheated* cells; `ṁ_vl` (vapour →
liquid) acts in *subcooled* cells. Each generates a corresponding latent
heat source in the energy equation:

    S_α_v  =  (ṁ_lv − ṁ_vl) / ρ_v                                  [1/s]
    S_T    =  −h_fg · (ṁ_lv − ṁ_vl) / (ρ_m · cp_m)                 [K/s]
                                       ^ evap cools, cond heats

# Inputs
- `C_e` evaporation rate coefficient. Typical 0.01–100 (case-tuned).
- `C_c` condensation rate coefficient. Typical 0.01–100 (case-tuned).

# Notes
- The model is for *bulk* volumetric phase change. Wall nucleation is
  handled separately by the RPI wall BC mechanism (see `RPI.jl`).
- Coefficients are empirical: validation against the target experiment
  usually requires sweeping `C_e`/`C_c` in the 0.1–10 range.
- Both branches require α_l > 0 and α_v > 0 respectively to avoid
  generating phases that don't exist.

# Convention
This file follows the user's existing convention where **α tracks the
LIQUID volume fraction** (matches `Dirichlet(:Inlet, 1.0)` in HeatFlux.jl):
- α = α_l, α_v = 1 − α_l
- Evaporation:    S_α_l = −ṁ_lv / ρ_l  (liquid disappears)
- Condensation:   S_α_l = +ṁ_vl / ρ_l  (liquid appears)
"""
Base.@kwdef struct Lee{F<:AbstractFloat}
    C_e::F = 0.1
    C_c::F = 0.1
end
Adapt.@adapt_structure Lee

# Pass-through for multiphase build pipeline.
build_property(lee::Lee, mesh) = lee


# =============================================================================
# Bulk phase-change source kernel
# =============================================================================

"""
    _lee_sources!(S_alpha_l, S_T, T, alpha_l, T_sat, h_fg, rho_l, rho_v,
                  rho_m, cp_m, C_e, C_c)

Per-cell Lee phase-change source kernel. Fills:
- `S_alpha_l` : α_l equation source [1/s]  (added as +Source)
- `S_T`       : T  equation source [K/s]  (added as +Source)

All `T_sat`, `h_fg`, `rho_l`, `rho_v` are scalars (problem-setup
constants from the operating pressure); `T`, `alpha_l`, `rho_m`, `cp_m`
are per-cell ScalarFields.
"""
@kernel inbounds=true function _lee_sources!(
        S_alpha_l, S_T,
        T, alpha_l,
        T_sat, h_fg, rho_l, rho_v, C_e, C_c)
    i = @index(Global)

    a_l = alpha_l[i]
    a_v = one(eltype(alpha_l)) - a_l
    Ti  = T[i]

    # Bounded driving potentials
    super = max((Ti - T_sat) / T_sat, zero(eltype(T)))      # evap if T > T_sat
    sub   = max((T_sat - Ti) / T_sat, zero(eltype(T)))      # cond if T < T_sat

    # Mass-transfer rates  [kg/(m³·s)]
    mdot_lv = C_e * max(a_l, zero(eltype(alpha_l))) * rho_l * super
    mdot_vl = C_c * max(a_v, zero(eltype(alpha_l))) * rho_v * sub

    # α_l source [1/s]: liquid disappears when evaporating, appears when condensing
    # ∂α_l/∂t = (−ṁ_lv + ṁ_vl) / ρ_l   per unit volume
    S_alpha_l[i] = (-mdot_lv + mdot_vl) / rho_l

    # Energy source [W/m³]: evap cools (−), cond heats (+) — direct volumetric power.
    # The ρ_m·cp_m·∂T/∂t equation absorbs this naturally; no division by ρ·cp here.
    S_T[i] = -h_fg * (mdot_lv - mdot_vl)
end


"""
    apply_lee_sources!(S_alpha, S_T, lee::Lee, sat, model, config)

Driver that launches the Lee kernel and fills the source ScalarFields.
`sat` is a `SatProps` (see sat_props.jl) carrying T_sat, h_fg, rho_l, rho_v.
"""
function apply_lee_sources!(
        S_alpha_l::ScalarField, S_T::ScalarField,
        lee::Lee, sat,
        model, config)
    (; backend, workgroup) = config.hardware
    T_field     = model.energy.T
    alpha_field = model.fluid.alpha

    F = eltype(T_field.values)

    ndrange = length(T_field.values)
    kernel! = _lee_sources!(_setup(backend, workgroup, ndrange)...)
    kernel!(
        S_alpha_l.values, S_T.values,
        T_field.values, alpha_field.values,
        F(sat.T_sat), F(sat.h_fg), F(sat.rho_l), F(sat.rho_v),
        F(lee.C_e), F(lee.C_c),
    )
    return nothing
end

# No-op when phase_change is not set (graceful default).
apply_lee_sources!(_S_alpha, _S_T, ::Nothing, _sat, _model, _config) = nothing
