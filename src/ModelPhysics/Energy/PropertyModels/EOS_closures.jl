using XCALibre


_eos_wrapper(fluid::H2, T, P, alpha) = XCALibre.ModelPhysics.EOS_wrapper_H2(T, P, alpha)
_eos_wrapper(fluid::N2, T, P, alpha) = XCALibre.ModelPhysics.EOS_wrapper_N2(T, P, alpha)

_mu_high_fidelity(fluid::H2, T, rho) = XCALibre.ModelPhysics.mu_high_fidelity_H2(T, rho)
_mu_high_fidelity(fluid::N2, T, rho) = XCALibre.ModelPhysics.mu_high_fidelity_N2(T, rho)

_thermal_conductivity(fluid::H2, args...) = XCALibre.ModelPhysics.thermal_conductivity_H2(args...)
_thermal_conductivity(fluid::N2, args...) = XCALibre.ModelPhysics.thermal_conductivity_N2(args...)


function (eos::HelmholtzEnergy)(T_input, P_input, alpha_input)
    is_mp, rho0, cv0, cp0, kT0, kT_ref, internal_energy0, 
            enthalpy0, entropy0, latentHeat0, T_sat0 = _eos_wrapper(eos.name, T_input, P_input, alpha_input)
        
    # DEPENDING ON IS_MP, DEFINE ARRAYS....
    if !is_mp
        nu_bar = _mu_high_fidelity(eos.name, T_input, rho0)
        k0 = _thermal_conductivity(eos.name, rho0, T_input, cp0, cv0, kT0, kT_ref, nu_bar)
        surface_tension = calculate_surface_tension(eos.name, T_input)
        
        return (is_mp=is_mp, rho=rho0, cv=cv0, cp=cp0, u=internal_energy0, h=enthalpy0, s=entropy0, 
                mu=nu_bar, k=k0, sigma=surface_tension, L_vap=latentHeat0, T_sat=T_sat0)
    else
        nu_bar_vals = [0.0, 0.0]
        k0_vals = [0.0, 0.0]

        for i in eachindex(rho0)
            nu_bar_vals[i] = _mu_high_fidelity(eos.name, T_input, rho0[i])
            k0_vals[i] = _thermal_conductivity(eos.name, rho0[i], T_input, cp0[i], cv0[i], kT0[i], kT_ref[i], nu_bar_vals[i])
        end

        surface_tension = calculate_surface_tension(eos.name, T_input)

        return (is_mp=is_mp, rho=rho0, cv=cv0, cp=cp0, u=internal_energy0, h=enthalpy0, s=entropy0, 
                mu=nu_bar_vals, k=k0_vals, sigma=surface_tension, L_vap=latentHeat0, T_sat=T_sat0)
    end
end



