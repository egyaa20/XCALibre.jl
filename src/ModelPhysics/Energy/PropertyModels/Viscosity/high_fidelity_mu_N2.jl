export mu_high_fidelity_N2

### NOTE: this function works also with Oxygen, Air, Argon; different coeffs are required though

###Refer to "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air", 2004

# Constants for Nitrogen (H2)

M  = 28.01348
sigma =  0.3656
epsilon_div_kb =  98.94
T_c   = 126.192
rho_c = 11.1839


# Tables of coefficients
b = ( 0.431, −0.4623, 0.08406, 0.005341, −0.00331 )

N_coeffs = (10.72, 0.03989, 0.001208, -7.402, 4.620)
t_coeffs = (0.1, 0.25, 3.2, 0.9, 0.3)
d_coeffs = (2, 10, 12, 2, 1)
l_coeffs = (0, 1, 1, 2, 3)


function mu_0(T::Float64)
    Tstar = T / epsilon_div_kb
    ln_Tstar = log(Tstar)
    
    ln_sum = b[1] + b[2]*ln_Tstar + b[3]*ln_Tstar^2 + b[4]*ln_Tstar^3 + b[5]*ln_Tstar^4
    omega  = exp(ln_sum)
    
    return ( 0.0266958 * sqrt(M * T) ) / ( (sigma)^2 * omega )
end

function mu_r(δ::Float64, τ::Float64)
    residual_viscosity = 0.0

    for i in eachindex(N_coeffs)
        
        gamma = 0.0
        
        if l_coeffs[i] == 0
            # If the exponent l is zero, gamma is zero
            gamma = 0.0
        else
            # If the exponent l is not zero, gamma is one
            gamma = 1.0
        end

        term = N_coeffs[i] * (τ^t_coeffs[i]) * (δ^d_coeffs[i]) * exp(-gamma * (δ^l_coeffs[i]))
        residual_viscosity += term
    end

    return residual_viscosity
end

function mu_high_fidelity_N2(T::Float64, rho::Float64)
    rho = rho / M

    τ = T_c / T
    δ = rho / rho_c

    return mu_0(T) + mu_r(δ, τ)
end


# println(mu_high_fidelity_N2(78.151, 825.72))

# Great result!