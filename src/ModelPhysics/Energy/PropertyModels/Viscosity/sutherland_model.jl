export mu_sutherland

const SUTHERLAND_CONSTS = Dict(
    :air   => (mu_ref = 1.8e-5,  T_ref = 273.0,  S = 110.4)
)



function get_sutherland_constants(fluid::Symbol)
    params = get(SUTHERLAND_CONSTS, fluid, nothing)
    
    return params
end

function mu_sutherland(T; fluid::Symbol = :air)
    constants = get_sutherland_constants(fluid)
    (; mu_ref, T_ref, S) = constants

    C = mu_ref * (T/T_ref)^1.5
    mu = C * ((T_ref + S)/(T + S))

    return mu
end