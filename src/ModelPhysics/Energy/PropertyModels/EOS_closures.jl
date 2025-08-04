using XCALibre

export EOS_H2, EOS_N2


struct EOS_functions
    fluid_symbol::Symbol
    eos_wrapper::Function
    mu_high_fidelity::Function
    thermal_conductivity::Function
end


function (eos::EOS_functions)(T_input, P_input)
    is_mp, rho0, cv0, cp0, kT0, kT_ref, internal_energy0, 
            enthalpy0, entropy0, latentHeat0, T_sat0 = eos.eos_wrapper(T_input, P_input, 1.0) #alpha
    
    # DEPENDING ON IS_MP, DEFINE ARRAYS....
    if (is_mp == false)
        nu_bar = eos.mu_high_fidelity(T_input, rho0)
        k0 = eos.thermal_conductivity(rho0, T_input, cp0, cv0, kT0, kT_ref, nu_bar)
        surface_tension = calculate_surface_tension(eos.fluid_symbol, T_input)
        
        return is_mp, rho0, cv0, cp0, internal_energy0, enthalpy0, entropy0, 
                    nu_bar, k0, surface_tension, latentHeat0, T_sat
    else
        nu_bar_vals = [0.0, 0.0]
        k0_vals = [0.0, 0.0]

        for i in eachindex(rho0)
            nu_bar_vals[i] = eos.mu_high_fidelity(T_input, rho0[i])
            k0_vals[i] = eos.thermal_conductivity(rho0[i], T_input, cp0[i], cv0[i], kT0[i], kT_ref[i], nu_bar_vals[i])
        end

        surface_tension = calculate_surface_tension(eos.fluid_symbol, T_input)

        return is_mp, rho0, cv0, cp0, internal_energy0, enthalpy0, entropy0, 
                nu_bar_vals, k0_vals, surface_tension, latentHeat0, T_sat
    end
end

const EOS_H2 = EOS_functions(
    :hydrogen,
    XCALibre.ModelPhysics.EOS_wrapper_H2,
    XCALibre.ModelPhysics.mu_high_fidelity_H2,
    XCALibre.ModelPhysics.thermal_conductivity_H2
)

const EOS_N2 = EOS_functions(
    :nitrogen,
    XCALibre.ModelPhysics.EOS_wrapper_N2,
    XCALibre.ModelPhysics.mu_high_fidelity_N2,
    XCALibre.ModelPhysics.thermal_conductivity_N2
)


T_input = 18.803  # Temperature in K
P_input = 0.1e6

is_mp, rho0, cv0, cp0, internal_energy0, enthalpy0, entropy0, nu_bar, k0, surface_tension = EOS_H2(T_input, P_input)
println("rho: $rho0, is_mp: $is_mp, cv: $cv0, mu: $nu_bar, k: $k0")

# Works well, but vapour side thermal conductivity is questionable - 5% error....


# export EOS_H2
# export EOS_N2

# function EOS_H2()
#     rho0, cv0, cp0, kT0, kT_ref, internal_energy0, enthalpy0, entropy0 = XCALibre.ModelPhysics.EOS_wrapper_H2(T_input, P_input)
#     nu_bar = XCALibre.ModelPhysics.mu_high_fidelity_H2(T_input, rho0)
#     k0 = XCALibre.ModelPhysics.thermal_conductivity_H2(rho0, T_input, cp0, cv0, kT0, kT_ref, nu_bar)

#     return rho0, cv0, cp0, internal_energy0, enthalpy0, entropy0, nu_bar, k0
# end


# function EOS_H2()
#     rho0, cv0, cp0, kT0, kT_ref, internal_energy0, enthalpy0, entropy0 = XCALibre.ModelPhysics.EOS_wrapper_N2(T_input, P_input)
#     nu_bar = XCALibre.ModelPhysics.mu_high_fidelity_N2(T_input, rho0)
#     k0 = XCALibre.ModelPhysics.thermal_conductivity_N2(rho0, T_input, cp0, cv0, kT0, kT_ref, nu_bar)

#     return rho0, cv0, cp0, internal_energy0, enthalpy0, entropy0, nu_bar, k0
# end