export SatProps, sat_props_N2_at

"""
    SatProps

Saturation-line properties at a specific operating pressure. Holds the
scalars needed by Lee / RPI boiling models. Built once at problem setup
from a Helmholtz EOS evaluation (or hardcoded from NIST).

Fields (all SI):
- `pressure` : operating pressure [Pa]
- `T_sat`    : saturation temperature at `pressure` [K]
- `h_fg`     : latent heat of vapourisation [J/kg]
- `rho_l`    : saturated liquid density [kg/m³]
- `rho_v`    : saturated vapour density [kg/m³]
- `cp_l`     : saturated liquid specific heat [J/(kg·K)]
- `mu_l`     : saturated liquid viscosity [Pa·s]
- `k_l`      : saturated liquid thermal conductivity [W/(m·K)]
- `sigma`    : surface tension at `T_sat` [N/m]    (0 if unknown)
"""
Base.@kwdef struct SatProps{F<:AbstractFloat}
    pressure::F
    T_sat::F
    h_fg::F
    rho_l::F
    rho_v::F
    cp_l::F
    mu_l::F
    k_l::F
    sigma::F = zero(F)
end
Adapt.@adapt_structure SatProps

build_property(s::SatProps, mesh) = s


# -----------------------------------------------------------------------------
# Hardcoded NIST values for N2 at common subcritical pressures.
# Source: NIST REFPROP. Pulled directly to avoid runtime Helmholtz evaluation
# for the operating-point constants. If you need a pressure not listed here,
# add an entry.
# -----------------------------------------------------------------------------

const _SAT_N2_TABLE = Dict{Float64, SatProps{Float64}}(
    #            P [Pa]       T_sat   h_fg     rho_l   rho_v   cp_l   mu_l         k_l       sigma
    0.5e6  => SatProps(0.5e6,  91.4,  192.6e3, 763.0,  21.0,  2080.0, 110.0e-6,  0.1198,   6.8e-3),
    1.0e6  => SatProps(1.0e6, 103.7,  161.4e3, 712.0,  41.4,  2200.0, 76.0e-6,   0.1053,   4.0e-3),
    1.5e6  => SatProps(1.5e6, 112.0,  133.8e3, 670.0,  62.5,  2410.0, 59.0e-6,   0.0942,   2.6e-3),
    2.0e6  => SatProps(2.0e6, 118.5,  108.1e3, 631.0,  85.4,  2680.0, 48.0e-6,   0.0853,   1.7e-3),
    2.5e6  => SatProps(2.5e6, 123.9,   80.7e3, 593.0, 113.0,  3140.0, 39.0e-6,   0.0782,   1.0e-3),
    3.0e6  => SatProps(3.0e6, 128.5,   45.0e3, 543.0, 154.0,  4500.0, 30.0e-6,   0.0729,   3.0e-4),
)

"""
    sat_props_N2_at(pressure_Pa) -> SatProps

Look up saturation properties for N2 at one of the tabulated pressures.
Errors with a helpful message if `pressure_Pa` isn't an exact key — at
the moment no interpolation is done because boiling cases use a single
fixed operating P. Add entries to `_SAT_N2_TABLE` as needed.
"""
function sat_props_N2_at(pressure_Pa::Real)
    key = float(pressure_Pa)
    haskey(_SAT_N2_TABLE, key) ||
        error("SatProps: no tabulated N2 saturation data for P=$key Pa.\n" *
              "Available: $(sort(collect(keys(_SAT_N2_TABLE)))).\n" *
              "Add an entry to _SAT_N2_TABLE in sat_props.jl, or pass a manually-" *
              "constructed `SatProps` instead.")
    return _SAT_N2_TABLE[key]
end
