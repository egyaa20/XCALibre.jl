export mu_sutherland


struct SutherlandLaw
    mu_ref::Float64  # Pa*s
    T_ref::Float64   # K
    S::Float64       # K
end

get_sutherland_law(::Type{Air}) = SutherlandLaw(1.8e-5, 273.0, 110.4)


function (law::SutherlandLaw)(T::Real)
    C = law.mu_ref * (T / law.T_ref)^1.5
    mu = C * ((law.T_ref + law.S) / (T + law.S))
    return mu
end

function mu_sutherland(T::Real, fluid::AbstractFluid)
    law = get_sutherland_law(typeof(fluid))
    
    return law(T)
end


# --- Step 3: Implement the functor ---
# This makes objects of type SutherlandLaw callable like functions.
# When you call an instance of SutherlandLaw with a temperature, this method runs.
# `(law::SutherlandLaw)` is the special syntax for defining a functor.
function (law::SutherlandLaw)(T::Real)
    # The calculation logic is the same, but it now uses the fields
    # from the struct instance `law`.
    # This couples the data (the constants) with the operation (the calculation).
    C = law.mu_ref * (T / law.T_ref)^1.5
    mu = C * ((law.T_ref + law.S) / (T + law.S))

    return mu
end



function mu_sutherland(T; fluid::Symbol = :air)
    model = get(SUTHERLAND_MODELS, fluid, nothing)
    
    if model === nothing
        error("Fluid '$fluid' not found. Available fluids are: $(keys(SUTHERLAND_MODELS))")
    end

    return model(T)
end





# const SUTHERLAND_CONSTS = Dict(
#     :air   => (mu_ref = 1.8e-5,  T_ref = 273.0,  S = 110.4)
# )



# function get_sutherland_constants(fluid::Symbol)
#     params = get(SUTHERLAND_CONSTS, fluid, nothing)
    
#     return params
# end

# function mu_sutherland(T; fluid::Symbol = :air)
#     constants = get_sutherland_constants(fluid)
#     (; mu_ref, T_ref, S) = constants

#     C = mu_ref * (T/T_ref)^1.5
#     mu = C * ((T_ref + S)/(T + S))

#     return mu
# end