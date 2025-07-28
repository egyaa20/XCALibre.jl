export calculate_surface_tension

###Refer to "Recommended Correlations for the Surface Tension of Common Fluids", 2012

FLUID_PROPERTIES = Dict(
    :hydrogen => (
        T_c = 32.938,
        σ_0 = 0.005314,
        n_0 = 1.060,
        T_min = 13.80, #temperature range before supercritical regime, then sigma=0
        T_max = 31.00  #temperature range before supercritical regime, then sigma=0
    ),
    :nitrogen => (
        T_c = 126.192,
        σ_0 = 0.02898,
        n_0 = 1.246,
        T_min = 64.80, #temperature range before supercritical regime, then sigma=0
        T_max = 120.24 #temperature range before supercritical regime, then sigma=0
    ),
)


function calculate_surface_tension(fluid::Symbol, T::Float64)
    if !haskey(FLUID_PROPERTIES, fluid)
        throw(KeyError("Properties for fluid <$fluid> are not defined."))
    end

    properties = FLUID_PROPERTIES[fluid]

    if !(properties.T_min <= T <= properties.T_max)
        # throw(DomainError(T, "Input temperature is outside the recommended range of $T_min K to $T_max K."))


        ### QUESTION: Just return sigma=0 ? (above critical)
        return 0.0
    end
    
    reduced_temp_term = 1.0 - (T / properties.T_c)
    
    surface_tension = properties.σ_0 * (reduced_temp_term ^ properties.n_0)

    return surface_tension
end




# test_temperature = 20.0
# st_h2 = calculate_surface_tension(:hydrogen, test_temperature)

# println("Calculated Surface Tension: $(round(st_h2, digits=5)) N/m")