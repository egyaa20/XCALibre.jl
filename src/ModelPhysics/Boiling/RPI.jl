export RPI
export bubble_departure_diameter_TK
export nucleation_site_density_LC
export bubble_departure_frequency_Cole
export rpi_partition
export rpi_solve_wall_temperature

"""
    RPI(; ...)

Configuration of the RPI (Kurul & Podowski, 1990) wall heat-flux
partitioning model with Lemmert-Chawla nucleation site density,
Tolubinsky-Kostanchuk bubble departure diameter, and Cole bubble
departure frequency.

The full partitioning is

    q_w = q_c + q_q + q_e                               (1)

with

    q_c = (1 − A_b) · h_c · (T_w − T_l)                 single-phase convection
    q_q = 2 · A_b · √(k_l ρ_l c_{p,l} f / π) · (T_w − T_l)   quenching
                                              [Del Valle-Kenning 1985]
    q_e = N'' · f · (π/6) D_bw³ · ρ_v · h_fg            evaporation

where `A_b = min(1, K · π D_bw² N'' / 4)` is the area fraction
influenced by departing bubbles. The associated mass-transfer rate
into the vapour phase is

    ṁ_e = q_e / h_fg                                     (kg/m²/s)

# Submodels (with their default constants from the literature)
- Tolubinsky-Kostanchuk (1970)  bubble departure diameter
    `D_bw = min(D_max, D_ref · exp(-ΔT_sub / ΔT_0))`
    Defaults: `D_ref = 0.6 mm`, `ΔT_0 = 45 K`, `D_max = 1.4 mm`.

- Lemmert-Chawla (1977)         nucleation site density
    `N'' = (m · ΔT_sup)^p`
    Defaults: `m = 210`, `p = 1.805`.

- Cole (1960)                   bubble departure frequency
    `f = √(4 g (ρ_l − ρ_v) / (3 D_bw ρ_l))`

# Inputs (keyword constructor)
- `D_ref`, `ΔT0_TK`, `D_max`     Tolubinsky-Kostanchuk constants
- `m_LC`, `p_LC`                 Lemmert-Chawla constants
- `K_influence`                  bubble-influence multiplier in `A_b`
- `g`                            gravitational acceleration magnitude

# API
Configured at the multiphase fluid level alongside `gravity`:

```julia
fluid = Fluid{Multiphase}(
    model    = Mixture(diameter=1e-6),
    phases   = (Phase(...), Phase(...)),
    gravity  = Gravity([0, -9.81, 0]),
    RPI      = RPI(),                        # defaults from literature
    # or RPI(K_influence=2.0, ΔT0_TK=30.0)   # tuned
)
```

The solver reads it back as `model.fluid.physics_properties.RPI` and
calls the partitioning functions per heated-wall face.
"""
Base.@kwdef struct RPI{F<:AbstractFloat}
    # Tolubinsky-Kostanchuk
    D_ref::F  = 6.0e-4
    ΔT0_TK::F = 45.0
    D_max::F  = 1.4e-3
    # Lemmert-Chawla
    m_LC::F   = 210.0
    p_LC::F   = 1.805
    # Bubble influence in quenching/convective area split
    K_influence::F = 4.0
    # Gravity magnitude (positive)
    g::F      = 9.81
end

# Pass-through for the multiphase build pipeline. `RPI` is a
# configuration record — it doesn't need any per-mesh allocation, so
# it passes straight into `model.fluid.physics_properties.RPI`.
build_property(rpi::RPI, mesh) = rpi


# =============================================================================
# Submodels — pure functions
# =============================================================================

"""
    bubble_departure_diameter_TK(ΔT_sub, model::RPI) -> D_bw

Tolubinsky-Kostanchuk departure diameter (m) given liquid subcooling
ΔT_sub = T_sat − T_l (K, ≥ 0). Saturated or super-heated bulk
(`ΔT_sub ≤ 0`) clamps the exponential to 1, giving `D_bw = D_max`.
"""
@inline function bubble_departure_diameter_TK(ΔT_sub::F, m::RPI) where {F<:Real}
    sub  = max(ΔT_sub, zero(F))
    dbw  = m.D_ref * exp(-sub / m.ΔT0_TK)
    return min(dbw, m.D_max)
end

"""
    nucleation_site_density_LC(ΔT_sup, model::RPI) -> N''

Lemmert-Chawla nucleation-site density (sites / m²) given wall
superheat ΔT_sup = T_w − T_sat (K). Returns 0 for non-positive
superheat (no boiling).
"""
@inline function nucleation_site_density_LC(ΔT_sup::F, m::RPI) where {F<:Real}
    ΔT_sup <= zero(F) && return zero(F)
    return (m.m_LC * ΔT_sup)^m.p_LC
end

"""
    bubble_departure_frequency_Cole(D_bw, ρ_l, ρ_v, model::RPI) -> f

Cole (1960) bubble departure frequency (1/s) given departure
diameter `D_bw` (m) and liquid / vapour densities (kg/m³).
Falls back to 0 when `D_bw <= 0` or the density jump is
non-positive (no buoyant bubble release).
"""
@inline function bubble_departure_frequency_Cole(D_bw::F, ρ_l::F, ρ_v::F,
                                                  m::RPI) where {F<:Real}
    (D_bw > zero(F) && ρ_l > ρ_v) || return zero(F)
    return sqrt(F(4.0) * m.g * (ρ_l - ρ_v) / (F(3.0) * D_bw * ρ_l))
end


# =============================================================================
# Heat-flux partitioning
# =============================================================================

"""
    rpi_partition(T_w, T_l, T_sat, h_c, fluid, model)

Compute the three RPI heat-flux components and the bubble-induced
mass-transfer rate at a wall face. All inputs in SI units.

# Inputs
- `T_w`    Wall temperature (K).
- `T_l`    Local liquid temperature (K).
- `T_sat`  Saturation temperature at the wall pressure (K).
- `h_c`    Single-phase convective HTC at the wall (W/m²/K).
- `fluid`  NamedTuple with saturation properties:
           `(ρ_l, ρ_v, μ_l, cp_l, k_l, h_fg, σ)` — all SI. Only `μ_l`
           and `σ` are unused by this combination of submodels but
           are accepted so the same NamedTuple can serve more
           closures later. Set them to zero if you don't have them.
- `model`  `RPI` with submodel constants.

# Returns
NamedTuple of:
- `q_c, q_q, q_e`    [W/m²]   convective, quenching, evaporative components
- `q_total`           [W/m²]   q_c + q_q + q_e
- `mdot_e`            [kg/m²/s] vapour generation rate per face area
- `A_b`               [-]       bubble-influence area fraction
- `D_bw`              [m]
- `N_pp`              [sites/m²]
- `f`                 [1/s]
- `ΔT_sub, ΔT_sup`    [K]
"""
function rpi_partition(T_w::F, T_l::F, T_sat::F, h_c::F,
                       fluid::NamedTuple, model::RPI) where {F<:AbstractFloat}
    ΔT_sup = T_w  - T_sat            # wall superheat
    ΔT_sub = T_sat - T_l             # liquid subcooling

    # Submodels
    D_bw = bubble_departure_diameter_TK(ΔT_sub, model)
    N_pp = nucleation_site_density_LC(ΔT_sup, model)
    f    = bubble_departure_frequency_Cole(D_bw, fluid.ρ_l, fluid.ρ_v, model)

    # Bubble-influence area fraction
    A_b = min(one(F), F(model.K_influence) * F(π) * D_bw^2 * N_pp / F(4.0))

    # Convective: portion of the wall not covered by bubble influence
    q_c = (one(F) - A_b) * h_c * (T_w - T_l)

    # Quenching: transient conduction during bubble departure rewetting
    # (Del Valle-Kenning 1985 form used in Kurul & Podowski 1990).
    quench_coeff = F(2.0) * A_b * sqrt(fluid.k_l * fluid.ρ_l * fluid.cp_l * f / F(π))
    q_q = quench_coeff * (T_w - T_l)

    # Evaporation: latent heat carried by departing bubbles
    V_b = F(π) / F(6.0) * D_bw^3
    q_e = N_pp * f * V_b * fluid.ρ_v * fluid.h_fg

    q_total = q_c + q_q + q_e
    mdot_e  = fluid.h_fg > zero(F) ? q_e / fluid.h_fg : zero(F)

    return (
        q_c     = q_c,
        q_q     = q_q,
        q_e     = q_e,
        q_total = q_total,
        mdot_e  = mdot_e,
        A_b     = A_b,
        D_bw    = D_bw,
        N_pp    = N_pp,
        f       = f,
        ΔT_sub  = ΔT_sub,
        ΔT_sup  = ΔT_sup,
    )
end


# =============================================================================
# Inverse problem — solve for T_w given prescribed wall heat flux q_w
# =============================================================================

"""
    rpi_solve_wall_temperature(q_w, T_l, T_sat, h_c, fluid, model;
                               T_w_guess=T_sat+5.0, tol=1e-3, max_iter=50)

Newton iteration for the wall temperature satisfying

    rpi_partition(T_w, ...).q_total = q_w

i.e. given a prescribed *heat-flux density* at a heated wall (W/m²),
find the wall temperature that produces that flux through the
combined convective + quenching + evaporative paths. Useful for
heat-flux-controlled boundary conditions.

Returns the converged `T_w` plus the partitioning NamedTuple at that
temperature. Errors if no convergence is reached in `max_iter` steps.
"""
function rpi_solve_wall_temperature(q_w::F, T_l::F, T_sat::F, h_c::F,
                                     fluid::NamedTuple, model::RPI;
                                     T_w_guess::F = T_sat + F(5.0),
                                     tol::F       = F(1.0e-3),
                                     max_iter::Int = 50) where {F<:AbstractFloat}
    T_w = T_w_guess
    res = rpi_partition(T_w, T_l, T_sat, h_c, fluid, model)
    for it in 1:max_iter
        f_val = res.q_total - q_w
        abs(f_val) < tol * max(abs(q_w), one(F)) && return (T_w = T_w, partition = res)

        # Numerical Jacobian: finite difference. Closed-form is messy
        # because A_b clamps and the quenching branch has D_bw / N''
        # interactions; FD is robust and 1-D.
        h_step = max(F(1.0e-3) * abs(T_w - T_sat), F(1.0e-3))
        res2   = rpi_partition(T_w + h_step, T_l, T_sat, h_c, fluid, model)
        df_dTw = (res2.q_total - res.q_total) / h_step
        abs(df_dTw) < F(1.0e-12) && error("RPI: Jacobian near-zero at iter=$it, T_w=$T_w")

        # Damp Newton step to keep T_w above T_l (boiling regime requires
        # superheat or at least non-negative ΔT). Cap the step magnitude.
        dT      = -f_val / df_dTw
        dT      = clamp(dT, -F(50.0), F(50.0))
        T_w_new = max(T_w + dT, T_l + F(1.0e-3))

        T_w = T_w_new
        res = rpi_partition(T_w, T_l, T_sat, h_c, fluid, model)
    end
    error("RPI: wall-temperature solve did not converge in $max_iter iterations " *
          "(last T_w=$T_w, q_total=$(res.q_total), q_w=$q_w)")
end
