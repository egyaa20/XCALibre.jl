T_c   = 126.192   # K
rho_c = 11.1839   # mol dm^-3
M_N2  = 28.01348  # g mol^-1
R_univ = 8.314472 # J (mol * K)^-1

a_nitro = [
    2.5,
   -12.76952708,
    -0.00784163, 
    -1.934819e-4,
    -1.247742e-5,
     6.678326e-8,
     1.012941,
     26.65788
]

N = [
    # Polynomial Part (k=1 to 6)
    0.924803575275, -0.492448489428, 0.661883336938, -0.192902649201e1, 
    -0.622469309629e-1, 0.349943957581,
    # Exponential Part (k=7 to 32)
    0.564857472498, -0.161720005987e1, -0.481395031883, 0.421150636384,
    -0.161962230825e-1, 0.172100994165, 0.735448924933e-2, 0.16807730535479e-1,
    -0.107626664179e-2, -0.137318088513e-1, 0.635466899859e-3, 0.304432279419e-2,
    -0.4357623366045e-1, -0.723174889316e-1, 0.389644315272e-1, -0.212201363910e-1,
    0.40882298181509e-2, -0.55199017984e-4, -0.462016716479e-1, -0.30031116011e-2,
    0.36882591208e-1, -0.255856846220e-2, 0.896915264558e-2, -0.44151337070350e-2,
    0.133722924858e-2, 0.26483249191957e-3,
    # Gaussian Part (k=33 to 36)
    0.19668894015e2, -0.209115600730e2, 0.1677883066989e-1, 0.2627675665274e4
]

d = [
    1.0, 1.0, 2.0, 2.0, 3.0, 3.0, # Polynomial
    1.0, 1.0, 1.0, 3.0, 3.0, 4.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 1.0, 2.0, 3.0, 
    4.0, 5.0, 8.0, 4.0, 5.0, 5.0, 8.0, 3.0, 5.0, 6.0, 9.0, # Exponential
    1.0, 1.0, 3.0, 2.0  # Gaussian
]

t = [
    0.25, 0.875, 0.5, 0.875, 0.375, 0.75, # Polynomial
    0.5, 0.75, 2.0, 1.25, 3.5, 1.0, 0.5, 3.0, 0.0, 2.75, 0.75, 2.5, 4.0, 6.0,
    6.0, 3.0, 3.0, 6.0, 16.0, 11.0, 15.0, 12.0, 12.0, 7.0, 4.0, 16.0, # Exponential
    0.0, 1.0, 2.0, 3.0 # Gaussian
]

p = [
    0, 0, 0, 0, 0, 0, # Polynomial
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, # Exponential
    0, 0, 0, 0 # Gaussian
]

α = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, # Polynomial & Exponential
    20.0, 20.0, 15.0, 25.0 # Gaussian
]

β = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, # Polynomial & Exponential
    325.0, 325.0, 300.0, 275.0 # Gaussian
]

γ = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, # Polynomial & Exponential
    1.16, 1.16, 1.13, 1.25 # Gaussian
]

D = [
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, # Polynomial & Exponential
    1.0, 1.0, 1.0, 1.0 # Gaussian
]

function alpha_0(δ::Float64, τ::Float64)
    summations = ( a_nitro[4] * (τ^(-1)) ) + ( a_nitro[5] * (τ^(-2)) ) + ( a_nitro[6] * (τ^(-3)) )

    last_term = a_nitro[7] * log(1 - exp(-a_nitro[8]*τ))

    return log(δ) + ( a_nitro[1] * log(τ) ) + a_nitro[2] + ( a_nitro[3] * τ ) + summations + last_term
end


function d_alpha_0_d_delta(δ::Float64, τ::Float64)
    return 1.0 / δ
end

function d2_alpha_0_d_delta2(δ::Float64, τ::Float64)
    return -1.0 / (δ^2)
end


function d_alpha_0_d_tau(δ::Float64, τ::Float64)
    last_term = (a_nitro[7]*a_nitro[8])/(exp(a_nitro[8]*τ)-1)
    return (a_nitro[1]/τ) + a_nitro[3] - (a_nitro[4]/(τ^2)) - ((2*a_nitro[5])/(τ^3)) - ((3*a_nitro[6])/(τ^4)) + last_term
end

function d2_alpha_0_d_tau2(δ::Float64, τ::Float64)
    exponent_term = exp(a_nitro[8]*τ)
    last_term = (a_nitro[7]*(a_nitro[8]^2)*exponent_term)/((exponent_term-1)^2)
    return -(a_nitro[1]/(τ^2)) + ((2*a_nitro[4])/(τ^3)) + ((6*a_nitro[5])/(τ^4)) + ((12*a_nitro[6])/(τ^5)) - last_term
end





### Residual Part
function alpha_r(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:6
        term1 +=  N[i] * (δ^d[i]) * (τ^t[i])
    end


    for i in 7:32
        term2 +=  N[i] * (δ^d[i]) * (τ^t[i]) * exp( -δ^p[i] )
    end


    for i in 33:36
        expo = -α[i]*(δ - D[i])^2 - β[i]*(τ-γ[i])^2

        term3  += N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
    end

    return term1 + term2 + term3
end


function d_alpha_r_d_tau(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:6
        if t[i] != 0
            term1 += N[i] * (δ^d[i]) * t[i] * (τ^(t[i]-1))
        end
    end

    for i in 7:32
        if t[i] != 0
            term2 += N[i] * (δ^d[i]) * t[i] * (τ^(t[i]-1)) * exp(-δ^p[i])
        end
    end

    for i in 33:36
        expo = -α[i]*(δ - D[i])^2 - β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        factor = (t[i] / τ) - 2 * β[i] * (τ - γ[i])
        term3 += term_i * factor
    end

    return term1 + term2 + term3
end

function d_alpha_r_d_delta(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:6
        if d[i] != 0
            term1 += N[i] * d[i] * (δ^(d[i]-1)) * (τ^t[i])
        end
    end

    for i in 7:32
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(-δ^p[i])
        factor = (d[i] / δ) - p[i] * (δ^(p[i]-1))
        term2 += term_i * factor
    end

    for i in 33:36
        expo = -α[i]*(δ - D[i])^2 - β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        factor = (d[i] / δ) - 2 * α[i] * (δ - D[i])
        term3 += term_i * factor
    end

    return term1 + term2 + term3
end


function d2_alpha_r_d_tau2(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:6
        term1 += N[i] * (δ^d[i]) * t[i] * (t[i]-1) * (τ^(t[i]-2))
    end

    for i in 7:32
        term2 += N[i] * (δ^d[i]) * t[i] * (t[i]-1) * (τ^(t[i]-2)) * exp(-δ^p[i])
    end

    for i in 33:36
        expo = -α[i]*(δ - D[i])^2 - β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        
        A_i = (t[i] / τ) - 2 * β[i] * (τ - γ[i])
        A_prime_i = -(t[i] / τ^2) - 2 * β[i]
        
        term3 += term_i * (A_i^2 + A_prime_i)
    end

    return term1 + term2 + term3
end

function d2_alpha_r_d_delta2(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:6
        term1 += N[i] * d[i] * (d[i]-1) * (δ^(d[i]-2)) * (τ^t[i])
    end

    for i in 7:32
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(-δ^p[i])
        B_i = (d[i] / δ) - p[i] * (δ^(p[i]-1))
        B_prime_i = -(d[i] / δ^2) - p[i] * (p[i]-1) * (δ^(p[i]-2))
        term2 += term_i * (B_i^2 + B_prime_i)
    end

    for i in 33:36
        expo = -α[i]*(δ - D[i])^2 - β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        C_i = (d[i] / δ) - 2 * α[i] * (δ - D[i])
        C_prime_i = -(d[i] / δ^2) - 2 * α[i]
        term3 += term_i * (C_i^2 + C_prime_i)
    end

    return term1 + term2 + term3
end


function d2_alpha_r_d_delta_d_tau(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:6
        term1 += N[i] * d[i] * t[i] * (δ^(d[i]-1)) * (τ^(t[i]-1))
    end

    for i in 7:32
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(-δ^p[i])
        B_i = (d[i] / δ) - p[i] * (δ^(p[i]-1))
        term2 += term_i * (t[i] / τ) * B_i
    end

    for i in 33:36
        expo = -α[i]*(δ - D[i])^2 - β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        A_i = (t[i] / τ) - 2 * β[i] * (τ - γ[i])
        C_i = (d[i] / δ) - 2 * α[i] * (δ - D[i])
        term3 += term_i * A_i * C_i
    end

    return term1 + term2 + term3
end

### Lambdas

function lambda_0_01(δ::Float64, τ::Float64)
    return δ * d_alpha_0_d_delta(δ, τ)
end
function lambda_0_02(δ::Float64, τ::Float64)
    return (δ^2) * d2_alpha_0_d_delta2(δ, τ)
end

function lambda_0_10(δ::Float64, τ::Float64)
    return τ * d_alpha_0_d_tau(δ, τ)
end
function lambda_0_20(δ::Float64, τ::Float64)
    return (τ^2) * d2_alpha_0_d_tau2(δ, τ)
end

function lambda_0_11(δ::Float64, τ::Float64)
    return 0.0 # delta * dau * d2alpha_0 / ddelta*dtau
end




function lambda_r_01(δ::Float64, τ::Float64)
    return δ * d_alpha_r_d_delta(δ, τ)
end
function lambda_r_02(δ::Float64, τ::Float64)
    return (δ^2) * d2_alpha_r_d_delta2(δ, τ)
end

function lambda_r_10(δ::Float64, τ::Float64)
    return τ * d_alpha_r_d_tau(δ, τ)
end
function lambda_r_20(δ::Float64, τ::Float64)
    return (τ^2) * d2_alpha_r_d_tau2(δ, τ)
end

function lambda_r_11(δ::Float64, τ::Float64)
    return δ * τ * d2_alpha_r_d_delta_d_tau(δ, τ)
end




function lambda_total_01(δ::Float64, τ::Float64)
    return lambda_0_01(δ, τ) + lambda_r_01(δ, τ)
end
function lambda_total_02(δ::Float64, τ::Float64)
    return lambda_0_02(δ, τ) + lambda_r_02(δ, τ)
end

function lambda_total_10(δ::Float64, τ::Float64)
    return lambda_0_10(δ, τ) + lambda_r_10(δ, τ)
end
function lambda_total_20(δ::Float64, τ::Float64)
    return lambda_0_20(δ, τ) + lambda_r_20(δ, τ)
end

function lambda_total_11(δ::Float64, τ::Float64)
    return lambda_0_11(δ, τ) + lambda_r_11(δ, τ)
end




### CORE CALCULATIONS SECTION

function find_density(T::Float64, P_target::Float64; max_iter=25, tol=1.0e-8) # DOES NOT SUPPORT GPU
    τ = T_c / T
    RT = R_univ * T

    # Initial guess for molar density using the Ideal Gas Law
    rho = P_target / RT

    for it in 1:max_iter
        δ = rho / rho_c

        Z = lambda_total_01(δ, τ)

        # EOS pressure
        p_calc = pressure(T, rho)
        # p_calc = rho * RT * Z

        # Residual for Newton
        f = p_calc - P_target

        if abs(f/P_target) < tol
            # DO NOT!!! Convert molar density (mol/dm^3) to mass density (kg/m^3)
            println("CONVERGED at iteration $it")
            # return rho * M_H2
            return rho
        end

        # Analytical derivative dp/drho
        dp_drho = RT * (1.0 + 2.0 * lambda_r_01(δ, τ) + lambda_r_02(δ, τ))

        # Newton step
        rho_new = rho - f / dp_drho
        if rho_new <= 0 || !isfinite(rho_new)
            error("DIVERGENCE at iteration $it")
        end
        rho = rho_new
    end

    error("Did not converge in $max_iter iterations!")
end


function pressure(T::Float64, rho::Float64)
    τ = T_c / T
    δ = rho / rho_c
    Z = lambda_total_01(δ, τ)
    
    return Z * rho * R_univ * T
end


function c_v(δ::Float64, τ::Float64)
    return ( - R_univ * lambda_total_20(δ, τ) )
end

function c_p(δ::Float64, τ::Float64)
    cv_term = c_v(δ, τ)

    numerator = (1.0 + lambda_r_01(δ, τ) - lambda_r_11(δ, τ))^2
    denominator = (1.0 + 2.0 * lambda_r_01(δ, τ) + lambda_r_02(δ, τ))

    fraction = numerator / denominator

    return cv_term + (R_univ * fraction)
end

function k_T(T::Float64, rho::Float64, δ::Float64, τ::Float64)
    term1 = rho*R_univ*T
    term2 = 1.0 + 2.0*lambda_r_01(δ, τ) + lambda_r_02(δ, τ)

    return 1.0 / (term1*term2)
end





function cv0_calc(T::Float64)
    τ = T_c / T
    
    return - R_univ * lambda_0_20(1.0, τ)
end

function u0_calc(T::Float64)
    τ = T_c / T
    
    return R_univ * T * lambda_0_10(1.0, τ)
end

function h0_calc(T::Float64)
    return u0_calc(T) + R_univ * T
end

function internal_energy_calc(T::Float64, δ::Float64, τ::Float64)
    u_ideal = u0_calc(T)
    
    u_residual = R_univ * T * lambda_r_10(δ, τ)

    return u_ideal + u_residual
end

function enthalpy_calc(T::Float64, δ::Float64, τ::Float64)
    h_ideal = h0_calc(T)
    
    h_residual = R_univ * T * (lambda_r_10(δ, τ) + lambda_r_01(δ, τ))
    return h_ideal + h_residual
end

function entropy_calc(δ::Float64, τ::Float64)
    term1 = lambda_total_10(δ, τ)
    term2 = alpha_0(δ, τ) + alpha_r(δ, τ)
    
    return R_univ * (term1 - term2)
end


function EOS_wrapper_N2(T::Float64, p::Float64)
    rho_mol = find_density(T, p)

    τ = T_c / T
    δ = rho_mol / rho_c

    cv_mol = c_v(δ, τ)
    cp_mol = c_p(δ, τ)
    kT = k_T(T, rho_mol, δ, τ)

    internal_energy_mol = internal_energy_calc(T, δ, τ)
    enthalpy_mol = enthalpy_calc(T, δ, τ)
    entropy_mol = entropy_calc(δ, τ)

    conversion_factor = 1.0 / M_N2

    rho = rho_mol * M_N2

    cv = cv_mol*conversion_factor
    cp = cp_mol*conversion_factor
    internal_energy = internal_energy_mol*conversion_factor
    enthalpy = enthalpy_mol*conversion_factor
    entropy = entropy_mol*conversion_factor



    return rho, cv, cp, kT, internal_energy, enthalpy, entropy
end



T_input = 68.151 # Temperature in K
P_input = 10.0e3 # Pressure in kPa


rho0, cv0, cp0, kT0, internal_energy0, enthalpy0, entropy0 = EOS_wrapper_N2(T_input, P_input)


println("DENSITY:")
println(rho0)
# println("CV:")
# println(cv0)
# println("CP:")
# println(cp0)
# println("KT:")
# println(kT0)
# println("Internal E:")
# println(internal_energy0)
# println("Enthalpy:")
# println(enthalpy0)
# println("Entropy:")
# println(entropy0)