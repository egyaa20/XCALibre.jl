export mu_high_fidelity

#Refer to "Correlation for the Viscosity of Normal Hydrogen Obtained from Symbolic Regression", 2013

# Constants for Hydrogen (H2)
const M         = 2.01588       # Molar mass
const sigma         = 0.297         # Length scale
const epsilon_div_kb  = 30.41         # Energy scale
const T_c       = 33.145        # Critical temperature
const rho_sc      = 90.5          # Symbolic-regression scaling rho


# Tables of coefficients
const a = ( 0.209630, -0.455274, 0.143602, -0.0335325, 0.00276981 )

const b = (-0.1870, 2.4871, 3.7151, -11.0972, 9.0965, -3.8292, 0.5166)

const c = ( 6.43449673,
            4.56334068e-2,
            0.232797868,
            0.958326120,
            0.127941189,
            0.363576595 )


function mu_0(T::Float64) # so-called 'Zero-density Viscosity'
    Tstar = T / epsilon_div_kb
    ln_Tstar = log(Tstar)
    
    ln_sum = a[1] + a[2]*ln_Tstar + a[3]*ln_Tstar^2 + a[4]*ln_Tstar^3 + a[5]*ln_Tstar^4
    S_starT_star = exp(ln_sum)
    
    return ( 0.021357 * sqrt(M * T) ) / ( (sigma)^2 * S_starT_star )
end


function beta_mu(T::Float64) # 'Second viscosity viral coefficient'
    Tstar = T / epsilon_div_kb
    
    Bstar = 0.0
    for i in eachindex(b)
        Bstar += b[i] * (Tstar^(-1))
    end
    
    return (sigma^3) * Bstar
end


function mu_1(T::Float64) # 'The initial-density coefficient of viscosity'
    return mu_0(T) * beta_mu(T)
end




function mu_high_fidelity(T::Float64, rho::Float64) # End result
    T_r  = T / T_c
    rho_r  = rho / rho_sc

    exp_term = exp( c[2]*T_r + (c[3]/(T_r)) + ((c[4]*rho_r^2)/(c[5]+T_r)) + c[6]*(rho_r^6))

    return mu_0(T) + ( mu_1(T) * rho ) + ( c[1] * (rho_r^2) ) * exp_term
end


