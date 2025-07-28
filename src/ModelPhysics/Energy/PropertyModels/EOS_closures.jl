
export EOS_H2, EOS_N2


struct EOS_functions
    eos_wrapper::Function
    mu_high_fidelity::Function
    thermal_conductivity::Function
end


function (eos::EOS_functions)(T_input, P_input)
    rho0, cv0, cp0, kT0, kT_ref, internal_energy0, enthalpy0, entropy0 = eos.eos_wrapper(T_input, P_input)
    nu_bar = eos.mu_high_fidelity(T_input, rho0)
    k0 = eos.thermal_conductivity(rho0, T_input, cp0, cv0, kT0, kT_ref, nu_bar)

    return rho0, cv0, cp0, internal_energy0, enthalpy0, entropy0, nu_bar, k0
end

const EOS_H2 = EOS_functions(
    XCALibre.ModelPhysics.EOS_wrapper_H2,
    XCALibre.ModelPhysics.mu_high_fidelity_H2,
    XCALibre.ModelPhysics.thermal_conductivity_H2
)

const EOS_N2 = EOS_functions(
    XCALibre.ModelPhysics.EOS_wrapper_N2,
    XCALibre.ModelPhysics.mu_high_fidelity_N2,
    XCALibre.ModelPhysics.thermal_conductivity_N2
)



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