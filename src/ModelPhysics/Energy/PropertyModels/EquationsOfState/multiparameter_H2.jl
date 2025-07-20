T_c = 32.938       # K
rho_c = 15.538       # mol dm^-3
R_univ = 8.314472 # J (mol * K)^-1
M_H2 = 2.01588    # g mol^-1

a_para = [
    2.5,
    -1.4485891134,
    1.884521239,
    4.30256,
    13.0289,
    -47.7365,
    50.0013,
    -18.6261,
    0.993973,
    0.536078
]

k_para = [
    0.0,
    0.0,
    0.0,
    499.0,
    826.5,
    970.8,
    1166.2,
    1341.4,
    5395.0,
    10185.0
]

N = [-7.33375, 0.01, 2.60375, 4.66279, 0.68239,
           -1.47078, 0.135801, -1.05327, 0.328239,
           -0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.11951]

t = [0.6855, 1.0, 1.0, 0.489, 0.774, 1.133, 1.386,
           1.619, 1.162, 3.96, 5.276, 0.99, 6.791, 3.19]

d = [1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1]

p = [0,0,0,0,0,0,0, 1,1, 0,0,0,0,0]

α = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   -1.7437, -0.5516, -0.0634, -2.1341, -1.777
]

β = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   -0.194, -0.2019, -0.0301, -0.2383, -0.3253
]

γ = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.8048, 1.5248, 0.6648, 0.6832, 1.493
]

D = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.5487, 0.1785, 1.28, 0.6319, 1.7104
]


### PREPARATORY FUNCTIONS


function alpha_0(δ::Float64, τ::Float64)
    sum_val = 0.0
    for i in 4:10
        sum_val += a_para[i] * log(1.0 - exp((-k_para[i] * τ)/T_c))
    end
    return log(δ) + ( (a_para[1] - 1) * log(τ) ) + a_para[2] + ( a_para[3] * τ ) + sum_val
end

function d_alpha_0_d_delta(δ::Float64, τ::Float64)
    return 1.0 / δ
end

function d2_alpha_0_d_delta2(δ::Float64, τ::Float64)
    return -1.0 / (δ^2)
end

function d_alpha_0_d_tau(δ::Float64, τ::Float64)
    sum_val = 0.0
    for i in 4:10
        sum_val += (a_para[i] * k_para[i])/(T_c*(exp( (k_para[i]*τ)/T_c ) - 1))
    end
    return ((a_para[1] - 1)/τ) + a_para[3] + sum_val
end

function d2_alpha_0_d_tau2(δ::Float64, τ::Float64)
    sum_val = 0.0
    for i in 4:10
        sum_val += (a_para[i] * ( k_para[i]/T_c )^2) * (exp( (k_para[i]*τ)/T_c )/((exp( (k_para[i]*τ)/T_c ) - 1)^2))
    end
    return -((a_para[1] - 1)/(τ^2)) - sum_val
end


### Residual Part
function alpha_r(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:7
        term1 +=  N[i] * (δ^d[i]) * (τ^t[i])
    end


    for i in 8:9
        term2 +=  N[i] * (δ^d[i]) * (τ^t[i]) * exp( -δ^p[i] )
    end


    for i in 10:14
        expo = α[i]*(δ - D[i])^2 + β[i]*(τ-γ[i])^2

        term3  += N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
    end

    return term1 + term2 + term3
end


function d_alpha_r_d_tau(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:7
        if t[i] != 0
            term1 += N[i] * (δ^d[i]) * t[i] * (τ^(t[i]-1))
        end
    end

    for i in 8:9
        if t[i] != 0
            term2 += N[i] * (δ^d[i]) * t[i] * (τ^(t[i]-1)) * exp(-δ^p[i])
        end
    end

    for i in 10:14
        expo = α[i]*(δ - D[i])^2 + β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        factor = (t[i] / τ) + 2 * β[i] * (τ - γ[i])
        term3 += term_i * factor
    end

    return term1 + term2 + term3
end

function d_alpha_r_d_delta(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:7
        if d[i] != 0
            term1 += N[i] * d[i] * (δ^(d[i]-1)) * (τ^t[i])
        end
    end

    for i in 8:9
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(-δ^p[i])
        factor = (d[i] / δ) - p[i] * (δ^(p[i]-1))
        term2 += term_i * factor
    end

    for i in 10:14
        expo = α[i]*(δ - D[i])^2 + β[i]*(τ-γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        factor = (d[i] / δ) + 2 * α[i] * (δ - D[i])
        term3 += term_i * factor
    end

    return term1 + term2 + term3
end


function d2_alpha_r_d_tau2(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:7
        term1 += N[i] * (δ^d[i]) * t[i] * (t[i]-1) * (τ^(t[i]-2))
    end

    for i in 8:9
        term2 += N[i] * (δ^d[i]) * t[i] * (t[i]-1) * (τ^(t[i]-2)) * exp(-δ^p[i])
    end

    for i in 10:14
        expo = α[i]*(δ - D[i])^2 + β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        
        A_i = (t[i] / τ) + 2 * β[i] * (τ - γ[i])
        A_prime_i = -(t[i] / τ^2) + 2 * β[i]
        
        term3 += term_i * (A_i^2 + A_prime_i)
    end

    return term1 + term2 + term3
end

function d2_alpha_r_d_delta2(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:7
        term1 += N[i] * d[i] * (d[i]-1) * (δ^(d[i]-2)) * (τ^t[i])
    end

    for i in 8:9
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(-δ^p[i])
        B_i = (d[i] / δ) - p[i] * (δ^(p[i]-1))
        B_prime_i = -(d[i] / δ^2) - p[i] * (p[i]-1) * (δ^(p[i]-2))
        term2 += term_i * (B_i^2 + B_prime_i)
    end

    for i in 10:14
        expo = α[i]*(δ - D[i])^2 + β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        C_i = (d[i] / δ) + 2 * α[i] * (δ - D[i])
        C_prime_i = -(d[i] / δ^2) + 2 * α[i]
        term3 += term_i * (C_i^2 + C_prime_i)
    end

    return term1 + term2 + term3
end


function d2_alpha_r_d_delta_d_tau(δ::Float64, τ::Float64)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0

    for i in 1:7
        term1 += N[i] * d[i] * t[i] * (δ^(d[i]-1)) * (τ^(t[i]-1))
    end

    for i in 8:9
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(-δ^p[i])
        B_i = (d[i] / δ) - p[i] * (δ^(p[i]-1))
        term2 += term_i * (t[i] / τ) * B_i
    end

    for i in 10:14
        expo = α[i]*(δ - D[i])^2 + β[i]*(τ - γ[i])^2
        term_i = N[i] * (δ^d[i]) * (τ^t[i]) * exp(expo)
        A_i = (t[i] / τ) + 2 * β[i] * (τ - γ[i])
        C_i = (d[i] / δ) + 2 * α[i] * (δ - D[i])
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
    return 0.0 # delta * tau * d2alpha_0 / ddelta*dtau
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
            # return rho * M_H2
            println("CONVERGED at iteration $it")
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


function EOS_wrapper(T::Float64, p::Float64)
    rho_mol = find_density(T, p)

    τ = T_c / T
    δ = rho_mol / rho_c

    cv_mol = c_v(δ, τ)
    cp_mol = c_p(δ, τ)
    kT = k_T(T, rho_mol, δ, τ)

    internal_energy_mol = internal_energy_calc(T, δ, τ)
    enthalpy_mol = enthalpy_calc(T, δ, τ)
    entropy_mol = entropy_calc(δ, τ)

    conversion_factor = 1.0 / M_H2

    rho = rho_mol * M_H2

    cv = cv_mol*conversion_factor
    cp = cp_mol*conversion_factor
    internal_energy = internal_energy_mol*conversion_factor
    enthalpy = enthalpy_mol*conversion_factor
    entropy = entropy_mol*conversion_factor



    return rho, cv, cp, kT, internal_energy, enthalpy, entropy
end



T_input = 28.80  # Temperature in K
P_input = 10.0e3 # Pressure in kPa


rho0, cv0, cp0, kT0, internal_energy0, enthalpy0, entropy0 = EOS_wrapper(T_input, P_input)

println("DENSITY:")
println(rho0)
println("CV:")
println(cv0)
println("CP:")
println(cp0)
println("KT:")
println(kT0)
println("Internal E:")
println(internal_energy0)
println("Enthalpy:")
println(enthalpy0)
println("Entropy:")
println(entropy0)