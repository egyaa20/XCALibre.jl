export EOS_wrapper_H2
# export constants_EoS_H2

### CONSTANTS AND COEFFICIENTS ###
# T_c = 32.938      # K
# rho_c = 15.538      # mol dm^-3
# R_univ = 8.314472 # J (mol * K)^-1
# M_H2 = 2.01588   # g mol^-1
# T_ref = 49.7175     # K
#
# a_para = [2.5, -1.4485891134, 1.884521239, 4.30256, 13.0289,
#           -47.7365, 50.0013, -18.6261, 0.993973, 0.536078]
#
# k_para = [0.0, 0.0, 0.0, 499.0, 826.5, 970.8, 1166.2, 1341.4, 5395.0, 10185.0]
#
# N = [-7.33375, 0.01, 2.60375, 4.66279, 0.68239, -1.47078, 0.135801,
#      -1.05327, 0.328239, -0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.11951]
#
# t = [0.6855, 1.0, 1.0, 0.489, 0.774, 1.133, 1.386, 1.619, 1.162, 3.96,
#      5.276, 0.99, 6.791, 3.19]
#
# d = [1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1]
#
# p = [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0]
#
# α = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#      -1.7437, -0.5516, -0.0634, -2.1341, -1.777]
#
# β = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#      -0.194, -0.2019, -0.0301, -0.2383, -0.3253]
#
# γ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#      0.8048, 1.5248, 0.6648, 0.6832, 1.493]
#
# D = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#      1.5487, 0.1785, 1.28, 0.6319, 1.7104]

struct constants_EoS_H2
    T_c::Float64
    rho_c::Float64
    R_univ::Float64
    M_H2::Float64
    T_ref::Float64
    a_para::Vector{Float64}
    k_para::Vector{Float64}
    N::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    p::Vector{Int}
    α::Vector{Float64}
    β::Vector{Float64}
    γ::Vector{Float64}
    D::Vector{Float64}
end




### PREPARATORY FUNCTIONS


function alpha_0(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; T_c, a_para, k_para) = constants
    sum_val = 0.0
    for i in 4:10
        sum_val += a_para[i] * log(1.0 - exp((-k_para[i] * τ)/T_c))
    end
    return log(δ) + ( (a_para[1] - 1) * log(τ) ) + a_para[2] + ( a_para[3] * τ ) + sum_val
end

function d_alpha_0_d_delta(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return 1.0 / δ
end

function d2_alpha_0_d_delta2(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return -1.0 / (δ^2)
end

function d_alpha_0_d_tau(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; T_c, a_para, k_para) = constants
    sum_val = 0.0
    for i in 4:10
        sum_val += (a_para[i] * k_para[i])/(T_c*(exp( (k_para[i]*τ)/T_c ) - 1))
    end
    return ((a_para[1] - 1)/τ) + a_para[3] + sum_val
end

function d2_alpha_0_d_tau2(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; T_c, a_para, k_para) = constants
    sum_val = 0.0
    for i in 4:10
        sum_val += (a_para[i] * ( k_para[i]/T_c )^2) * (exp( (k_para[i]*τ)/T_c )/((exp( (k_para[i]*τ)/T_c ) - 1)^2))
    end
    return -((a_para[1] - 1)/(τ^2)) - sum_val
end


### Residual Part
function alpha_r(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; N, d, t, p, α, β, γ, D) = constants
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


function d_alpha_r_d_tau(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; N, d, t, p, α, β, γ, D) = constants
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

function d_alpha_r_d_delta(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; N, d, t, p, α, β, γ, D) = constants
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


function d2_alpha_r_d_tau2(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; N, d, t, p, α, β, γ, D) = constants
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

function d2_alpha_r_d_delta2(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; N, d, t, p, α, β, γ, D) = constants
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


function d2_alpha_r_d_delta_d_tau(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; N, d, t, p, α, β, γ, D) = constants
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

function lambda_0_01(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return δ * d_alpha_0_d_delta(δ, τ, constants)
end
function lambda_0_02(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return (δ^2) * d2_alpha_0_d_delta2(δ, τ, constants)
end

function lambda_0_10(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return τ * d_alpha_0_d_tau(δ, τ, constants)
end
function lambda_0_20(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return (τ^2) * d2_alpha_0_d_tau2(δ, τ, constants)
end

function lambda_0_11(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return 0.0 # delta * tau * d2alpha_0 / ddelta*dtau
end




function lambda_r_01(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return δ * d_alpha_r_d_delta(δ, τ, constants)
end
function lambda_r_02(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return (δ^2) * d2_alpha_r_d_delta2(δ, τ, constants)
end

function lambda_r_10(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return τ * d_alpha_r_d_tau(δ, τ, constants)
end
function lambda_r_20(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return (τ^2) * d2_alpha_r_d_tau2(δ, τ, constants)
end

function lambda_r_11(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return δ * τ * d2_alpha_r_d_delta_d_tau(δ, τ, constants)
end




function lambda_total_01(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return lambda_0_01(δ, τ, constants) + lambda_r_01(δ, τ, constants)
end
function lambda_total_02(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return lambda_0_02(δ, τ, constants) + lambda_r_02(δ, τ, constants)
end

function lambda_total_10(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return lambda_0_10(δ, τ, constants) + lambda_r_10(δ, τ, constants)
end
function lambda_total_20(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return lambda_0_20(δ, τ, constants) + lambda_r_20(δ, τ, constants)
end

function lambda_total_11(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    return lambda_0_11(δ, τ, constants) + lambda_r_11(δ, τ, constants)
end




### CORE CALCULATIONS SECTION

function find_density(T::Float64, P_target::Float64, constants::constants_EoS_H2;
    max_iter=25, tol=1.0e-8) # DOES NOT SUPPORT GPU
    (; T_c, rho_c, R_univ, M_H2, T_ref, a_para, k_para, N, t, d, p, α, β, γ, D) = constants

    τ = T_c / T
    RT = R_univ * T

    # Initial guess for molar density using the Ideal Gas Law
    rho = P_target / RT

    for it in 1:max_iter
        δ = rho / rho_c

        # Z = lambda_total_01(δ, τ, constants)

        # EOS pressure
        p_calc = pressure(T, rho, constants)
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
        dp_drho = RT * (1.0 + 2.0 * lambda_r_01(δ, τ, constants) + lambda_r_02(δ, τ, constants))

        # Newton step
        rho_new = rho - f / dp_drho
        if rho_new <= 0 || !isfinite(rho_new)
            error("DIVERGENCE at iteration $it")
        end
        rho = rho_new
    end

    error("Did not converge in $max_iter iterations!")
end


function pressure(T::Float64, rho::Float64, constants::constants_EoS_H2)
    (; T_c, rho_c, R_univ) = constants
    τ = T_c / T
    δ = rho / rho_c
    Z = 1.0 + lambda_r_01(δ, τ, constants)
    return Z * rho * R_univ * T
end

function c_v(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; R_univ) = constants
    return -R_univ * lambda_total_20(δ, τ, constants)
end

function c_p(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; R_univ) = constants
    cv_term = c_v(δ, τ, constants)
    numerator = (1.0 + lambda_total_01(δ, τ, constants) - lambda_total_11(δ, τ, constants))^2
    denominator = 1.0 + 2.0 * lambda_total_01(δ, τ, constants) + lambda_total_02(δ, τ, constants)
    return cv_term + (R_univ * (numerator / denominator))
end

function k_T(T::Float64, rho::Float64, constants::constants_EoS_H2)
    (; T_c, rho_c, R_univ) = constants
    τ = T_c / T
    δ = rho / rho_c
    term1 = rho * R_univ * T
    term2 = 1.0 + 2.0 * lambda_r_01(δ, τ, constants) + lambda_r_02(δ, τ, constants)
    return 1.0 / (term1 * term2)
end

function cv0_calc(T::Float64, constants::constants_EoS_H2)
    (; T_c, R_univ) = constants
    τ = T_c / T
    return -R_univ * lambda_0_20(1.0, τ, constants)
end

function u0_calc(T::Float64, constants::constants_EoS_H2)
    (; T_c, R_univ) = constants
    τ = T_c / T
    return R_univ * T * lambda_0_10(1.0, τ, constants)
end

function h0_calc(T::Float64, constants::constants_EoS_H2)
    (; R_univ) = constants
    return u0_calc(T, constants) + R_univ * T
end

function internal_energy_calc(T::Float64, δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; R_univ) = constants
    u_ideal = u0_calc(T, constants)
    u_residual = R_univ * T * lambda_r_10(δ, τ, constants)
    return u_ideal + u_residual
end

function enthalpy_calc(T::Float64, δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; R_univ) = constants
    h_ideal = h0_calc(T, constants)
    h_residual = R_univ * T * (lambda_r_10(δ, τ, constants) + lambda_r_01(δ, τ, constants))
    return h_ideal + h_residual
end

function entropy_calc(δ::Float64, τ::Float64, constants::constants_EoS_H2)
    (; R_univ) = constants
    term1 = lambda_total_10(δ, τ, constants)
    term2 = alpha_0(δ, τ, constants) + alpha_r(δ, τ, constants)
    return R_univ * (term1 - term2)
end




function gibbs_free_energy(T::Float64, rho::Float64, constants::constants_EoS_H2)
    (; T_c, rho_c, R_univ) = constants
    
    # if rho <= 0.0
    #     return Inf
    # end

    τ = T_c / T
    δ = rho / rho_c

    alpha_val = alpha_0(δ, τ, constants) + alpha_r(δ, τ, constants)
    
    p_val = pressure(T, rho, constants)
    Z = p_val / (rho * R_univ * T)

    # Calculate dimensional Gibbs free energy
    g = (alpha_val + Z) * R_univ * T
    
    return g
end


function vapour_pressure_ancillary(T::Float64, constants::constants_EoS_H2)
    (; T_c) = constants
    
    # Critical pressure for parahydrogen from the paper (in Pa)
    p_c = 1.2858e6

    # Coefficients and exponents from Table 8 from the paper mentioned above
    N = [-4.87767, 1.03359, 0.82668, -0.129412]
    k = [1.0, 1.5, 2.65, 7.4]
    
    θ = 1.0 - T / T_c
    
    sum_val = 0.0
    for i in eachindex(N)
        sum_val += N[i] * (θ^k[i])
    end
    
    ln_pr = (T_c / T) * sum_val
    
    return exp(ln_pr) * p_c
end



function find_density_advanced(T::Float64, P_target::Float64, rho_guess::Float64, constants::constants_EoS_H2;
    max_iter=100, tol=1.0e-8) # DOES NOT SUPPORT GPU

    (; T_c, rho_c, R_univ, M_H2, T_ref, a_para, k_para, N, t, d, p, α, β, γ, D) = constants

    τ = T_c / T
    RT = R_univ * T

    # Initial guess for density can now account for discontinuity
    rho = rho_guess


    for it in 1:max_iter
        δ = rho / rho_c
        
        p_calc = pressure(T, rho, constants)
        f = p_calc - P_target

        if abs(f / P_target) < tol
            return rho
        end

        # Analytical derivative dp/drho
        dp_drho = RT * (1.0 + 2.0 * lambda_r_01(δ, τ, constants) + lambda_r_02(δ, τ, constants))
        
        # Newton-Raphson step
        step = f / dp_drho
        relaxation_coeff = 1.0
        rho_new = rho - ( relaxation_coeff * step )
        
        # Check for nonphysical results
        if rho_new <= 0 || !isfinite(rho_new)
            return NaN 
        end
        rho = rho_new
    end

    error("Did not converge in $max_iter iterations!")
end




function find_saturation_properties(T::Float64, constants::constants_EoS_H2; max_iter=20, tol=1.0e-7)
    (; T_c, rho_c, R_univ) = constants
    if T >= T_c # Ensure the temperature is in the valid range (below critical)
        error("Temperature must be below the critical temperature for saturation calculation.")
    end

    p_guess = vapour_pressure_ancillary(T, constants) # This function would vary depending on fluid

    # Secant method for iterative solver
    p_current = p_guess
    p_prev = p_guess * 0.999 # Slightly perturb the previous pressure value for first step

    g_diff_current = 0.0
    g_diff_prev = 0.0

    for it in 1:max_iter
        p_ideal_gas = p_current / (R_univ * T)
        p_multiplied =  2.0 * rho_c

        rho_l = find_density_advanced(T, p_current, p_multiplied, constants) # Higher guess for liquid
        rho_v = find_density_advanced(T, p_current, p_ideal_gas, constants) # Ideal gas guess for vapour

        # Check if roots were found
        if isnan(rho_l) || isnan(rho_v)
            error("Failed to find liquid or vapor density root at T=$T K, P=$p_current Pa...")
        end

        # Compute difference in gibbs free energy
        g_l = gibbs_free_energy(T, rho_l, constants)
        g_v = gibbs_free_energy(T, rho_v, constants)
        g_diff_current = g_l - g_v

        # Convergence check based on gibbs energy
        if abs(g_diff_current / g_l) < tol
            println("Saturation solver converged in $it iterations.")
            return (p_current, rho_l * 1000.0, rho_v * 1000.0) # Conversion into appropriate units
        end

        # Secant method step
        if it > 1
            if abs(g_diff_current - g_diff_prev) < 1e-9 # Avoid division by zero
                p_next = p_current * 1.001
            else
                p_next = p_current - g_diff_current * (p_current - p_prev) / (g_diff_current - g_diff_prev)
            end
            
            # Update values for next itration
            p_prev = p_current
            g_diff_prev = g_diff_current
            p_current = p_next
        else
            
            p_prev = p_current
            g_diff_prev = g_diff_current
            
            p_current = p_guess * 1.001
        end

        if p_current <= 0 # Check for -ve temps
            error("Solver stepped to a -ve pressure.")
        end
    end

    error("Saturation solver did not converge in $max_iter iterations.")
end


function EOS_wrapper_H2(T::Float64, pressure::Float64)
    constants = constants_EoS_H2(
        32.938, # T_c
        15.538e3, # rho_c, multiplied by e3 for convenience
        8.314472, # R_univ
        2.01588e-3, # M_H2, multiplied by e-3 for convenience
        49.7175, # T_ref
        
        [2.5, -1.4485891134, 1.884521239, 4.30256, 13.0289,
        -47.7365, 50.0013, -18.6261, 0.993973, 0.536078], # a_para

        [0.0, 0.0, 0.0, 499.0, 826.5, 970.8, 1166.2, 1341.4, 5395.0, 10185.0], # k_para

        [-7.33375, 0.01, 2.60375, 4.66279, 0.68239, -1.47078, 0.135801,
        -1.05327, 0.328239, -0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.11951], # N

        [0.6855, 1.0, 1.0, 0.489, 0.774, 1.133, 1.386, 1.619, 1.162, 3.96,
        5.276, 0.99, 6.791, 3.19], # t

        [1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1], # d

        [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], # p
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

        -1.7437, -0.5516, -0.0634, -2.1341, -1.777], # α
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

        -0.194, -0.2019, -0.0301, -0.2383, -0.3253], # β
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.8048, 1.5248, 0.6648, 0.6832, 1.493], # γ

        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.5487, 0.1785, 1.28, 0.6319, 1.7104] # D
    )

    (; T_c, rho_c, R_univ, M_H2, T_ref, a_para, k_para, N, t, d, p, α, β, γ, D) = constants


    rho_mol = 0.0

    # Firstly, account for supercritical fluid state
    if T >= constants.T_c # Maybe add fancier check "within tolerance"
        rho_guess = pressure / (constants.R_univ * T) # Ideal gas guess
        rho_mol = find_density_advanced(T, pressure, rho_guess, constants)

    else # Else it is not supercritical
        (Psat_Pa, rho_l_sat, rho_v_sat) = find_saturation_properties(T, constants)

        if abs(pressure - Psat_Pa) < 1e-3 # TWO PHASE REGION
            error("TWO PHASE REGION AT T=$T !")

        elseif pressure < Psat_Pa # VAPOUR REGION
            println("VAPOUR")
            rho_guess = pressure / (constants.R_univ * T) # Ideal gas guess
            rho_mol = find_density_advanced(T, pressure, rho_guess, constants)

        elseif pressure > Psat_Pa # LIQUID REGION
            println("LIQUID")
            rho_guess = 2.0 * constants.rho_c # Higher density guess
            rho_mol = find_density_advanced(T, pressure, rho_guess, constants)
        end

        if isnan(rho_mol)
            error("Failed to find a valid density for the given single-phase state.")
        end
    end

    # rho_mol = find_density(T, pressure, constants) ## OLD WAY

    τ = T_c / T
    δ = rho_mol / rho_c

    cv_mol = c_v(δ, τ, constants)
    cp_mol = c_p(δ, τ, constants)
    kT = k_T(T, rho_mol, constants)
    kT_ref = k_T(T_ref, rho_mol, constants)

    internal_energy_mol = internal_energy_calc(T, δ, τ, constants)
    enthalpy_mol = enthalpy_calc(T, δ, τ, constants)
    entropy_mol = entropy_calc(δ, τ, constants)

    conversion_factor = 1.0 / M_H2

    rho = rho_mol * M_H2

    cv = cv_mol*conversion_factor
    cp = cp_mol*conversion_factor
    internal_energy = internal_energy_mol*conversion_factor
    enthalpy = enthalpy_mol*conversion_factor
    entropy = entropy_mol*conversion_factor

    return rho, cv, cp, kT, kT_ref, internal_energy, enthalpy, entropy
end


# T_input = 35.0  # Temperature in K
# P_input = 10.0e6


# rho0, cv0, cp0, kT0, kT_ref, internal_energy0, enthalpy0, entropy0 = EOS_wrapper_H2(T_input, P_input)

# println("Rho: $rho0")

# T_input = 28.80  # Temperature in K
# P_input = 10.0e3 # Pressure in kPa


# rho0, cv0, cp0, kT0, kT_ref, internal_energy0, enthalpy0, entropy0 = EOS_wrapper(T_input, P_input)

# println("DENSITY:")
# println(rho0)
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