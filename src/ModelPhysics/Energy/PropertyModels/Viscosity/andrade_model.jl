export mu_andrade

const ANDRADE_CONSTS = Dict(
    :water => (B = 1.732e-6, C = 1863.0)
)

function get_andrade_constants(fluid::Symbol)
    params = get(ANDRADE_CONSTS, fluid, nothing)

    return params
end

function mu_andrade(T; fluid::Symbol = :water)
    (B, C) = get_andrade_constants(fluid)
    mu = B * exp(C / T)
    return mu
end