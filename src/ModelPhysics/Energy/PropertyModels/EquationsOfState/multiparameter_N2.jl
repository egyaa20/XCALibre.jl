export EOS_wrapper_N2

### CONSTANTS AND COEFFICIENTS ###
# T_c   = 126.192   # K
# rho_c = 11.1839   # mol dm^-3
# M_N2  = 28.01348  # g mol^-1
# R_univ = 8.314472 # J (mol * K)^-1
# 
# a_nitro = [
#     2.5,
#    -12.76952708,
#     -0.00784163, 
#     -1.934819e-4,
#     -1.247742e-5,
#      6.678326e-8,
#      1.012941,
#      26.65788
# ]
# 
# N = [
#     # Polynomial Part (k=1 to 6)
#     0.924803575275, -0.492448489428, 0.661883336938, -0.192902649201e1, 
#     -0.622469309629e-1, 0.349943957581,
#     # Exponential Part (k=7 to 32)
#     0.564857472498, -0.161720005987e1, -0.481395031883, 0.421150636384,
#     -0.161962230825e-1, 0.172100994165, 0.735448924933e-2, 0.16807730535479e-1,
#     -0.107626664179e-2, -0.137318088513e-1, 0.635466899859e-3, 0.304432279419e-2,
#     -0.4357623366045e-1, -0.723174889316e-1, 0.389644315272e-1, -0.212201363910e-1,
#     0.40882298181509e-2, -0.55199017984e-4, -0.462016716479e-1, -0.30031116011e-2,
#     0.36882591208e-1, -0.255856846220e-2, 0.896915264558e-2, -0.44151337070350e-2,
#     0.133722924858e-2, 0.26483249191957e-3,
#     # Gaussian Part (k=33 to 36)
#     0.19668894015e2, -0.209115600730e2, 0.1677883066989e-1, 0.2627675665274e4
# ] 






# d = [
#     1.0, 1.0, 2.0, 2.0, 3.0, 3.0, # Polynomial
#     1.0, 1.0, 1.0, 3.0, 3.0, 4.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 1.0, 2.0, 3.0, 
#     4.0, 5.0, 8.0, 4.0, 5.0, 5.0, 8.0, 3.0, 5.0, 6.0, 9.0, # Exponential
#     1.0, 1.0, 3.0, 2.0  # Gaussian
# ]
# 
# t = [
#     0.25, 0.875, 0.5, 0.875, 0.375, 0.75, # Polynomial
#     0.5, 0.75, 2.0, 1.25, 3.5, 1.0, 0.5, 3.0, 0.0, 2.75, 0.75, 2.5, 4.0, 6.0,
#     6.0, 3.0, 3.0, 6.0, 16.0, 11.0, 15.0, 12.0, 12.0, 7.0, 4.0, 16.0, # Exponential
#     0.0, 1.0, 2.0, 3.0 # Gaussian
# ]
# 
# p = [
#     0, 0, 0, 0, 0, 0, # Polynomial
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, # Exponential
#     0, 0, 0, 0 # Gaussian
# ]
# 
# α = [
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#     0.0, 0.0, # Polynomial & Exponential
#     20.0, 20.0, 15.0, 25.0 # Gaussian
# ]
# 
# β = [
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#     0.0, 0.0, # Polynomial & Exponential
#     325.0, 325.0, 300.0, 275.0 # Gaussian
# ]
# 
# γ = [
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#     0.0, 0.0, # Polynomial & Exponential
#     1.16, 1.16, 1.13, 1.25 # Gaussian
# ]
# 
# D = [
#     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
#     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
#     1.0, 1.0, # Polynomial & Exponential
#     1.0, 1.0, 1.0, 1.0 # Gaussian
# ]

struct constants_EoS_N2
    T_c::AbstractFloat
    rho_c::AbstractFloat
    M_N2::AbstractFloat
    R_univ::AbstractFloat
    T_ref::AbstractFloat
    a_nitro::Vector{AbstractFloat}
    N::Vector{AbstractFloat}
    d::Vector{AbstractFloat}
    t::Vector{AbstractFloat}
    p::Vector{Int}
    α::Vector{AbstractFloat}
    β::Vector{AbstractFloat}
    γ::Vector{AbstractFloat}
    D::Vector{AbstractFloat}
end

function alpha_0(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; a_nitro) = constants
    summations = ( a_nitro[4] * (τ^(-1)) ) + ( a_nitro[5] * (τ^(-2)) ) + ( a_nitro[6] * (τ^(-3)) )

    last_term = a_nitro[7] * log(1 - exp(-a_nitro[8]*τ))

    return log(δ) + ( a_nitro[1] * log(τ) ) + a_nitro[2] + ( a_nitro[3] * τ ) + summations + last_term
end


function d_alpha_0_d_delta(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return 1.0 / δ
end

function d2_alpha_0_d_delta2(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return -1.0 / (δ^2)
end


function d_alpha_0_d_tau(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; a_nitro) = constants
    last_term = (a_nitro[7]*a_nitro[8])/(exp(a_nitro[8]*τ)-1)
    return (a_nitro[1]/τ) + a_nitro[3] - (a_nitro[4]/(τ^2)) - ((2*a_nitro[5])/(τ^3)) - ((3*a_nitro[6])/(τ^4)) + last_term
end

function d2_alpha_0_d_tau2(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; a_nitro) = constants
    exponent_term = exp(a_nitro[8]*τ)
    last_term = (a_nitro[7]*(a_nitro[8]^2)*exponent_term)/((exponent_term-1)^2)
    return -(a_nitro[1]/(τ^2)) + ((2*a_nitro[4])/(τ^3)) + ((6*a_nitro[5])/(τ^4)) + ((12*a_nitro[6])/(τ^5)) - last_term
end





### Residual Part
function alpha_r(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; N, d, t, p, α, β, γ, D) = constants
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


function d_alpha_r_d_tau(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; N, d, t, p, α, β, γ, D) = constants
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

function d_alpha_r_d_delta(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; N, d, t, p, α, β, γ, D) = constants
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


function d2_alpha_r_d_tau2(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; N, d, t, p, α, β, γ, D) = constants
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

function d2_alpha_r_d_delta2(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; N, d, t, p, α, β, γ, D) = constants
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


function d2_alpha_r_d_delta_d_tau(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; N, d, t, p, α, β, γ, D) = constants
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

function lambda_0_01(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return δ * d_alpha_0_d_delta(δ, τ, constants)
end
function lambda_0_02(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return (δ^2) * d2_alpha_0_d_delta2(δ, τ, constants)
end

function lambda_0_10(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return τ * d_alpha_0_d_tau(δ, τ, constants)
end
function lambda_0_20(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return (τ^2) * d2_alpha_0_d_tau2(δ, τ, constants)
end

function lambda_0_11(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return 0.0 # delta * dau * d2alpha_0 / ddelta*dtau
end




function lambda_r_01(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return δ * d_alpha_r_d_delta(δ, τ, constants)
end
function lambda_r_02(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return (δ^2) * d2_alpha_r_d_delta2(δ, τ, constants)
end

function lambda_r_10(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return τ * d_alpha_r_d_tau(δ, τ, constants)
end
function lambda_r_20(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return (τ^2) * d2_alpha_r_d_tau2(δ, τ, constants)
end

function lambda_r_11(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return δ * τ * d2_alpha_r_d_delta_d_tau(δ, τ, constants)
end




function lambda_total_01(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return lambda_0_01(δ, τ, constants) + lambda_r_01(δ, τ, constants)
end
function lambda_total_02(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return lambda_0_02(δ, τ, constants) + lambda_r_02(δ, τ, constants)
end

function lambda_total_10(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return lambda_0_10(δ, τ, constants) + lambda_r_10(δ, τ, constants)
end
function lambda_total_20(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return lambda_0_20(δ, τ, constants) + lambda_r_20(δ, τ, constants)
end

function lambda_total_11(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    return lambda_0_11(δ, τ, constants) + lambda_r_11(δ, τ, constants)
end




### CORE CALCULATIONS SECTION

function find_density(T::AbstractFloat, P_target::AbstractFloat, constants::constants_EoS_N2;
    max_iter=25, tol=1.0e-8) # DOES NOT SUPPORT GPU
    (; T_c, rho_c, R_univ) = constants
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


function pressure(T::AbstractFloat, rho::AbstractFloat, constants::constants_EoS_N2)
    (; T_c, rho_c, R_univ) = constants
    τ = T_c / T
    δ = rho / rho_c
    
    Z = 1.0 + lambda_r_01(δ, τ, constants)
    return Z * rho * R_univ * T



    # τ = T_c / T
    # δ = rho / rho_c
    # Z = lambda_total_01(δ, τ, constants)
    
    # return Z * rho * R_univ * T
end


function c_v(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; R_univ) = constants
    return ( - R_univ * lambda_total_20(δ, τ, constants) )
end

function c_p(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; R_univ) = constants
    cv_term = c_v(δ, τ, constants)

    numerator = (1.0 + lambda_r_01(δ, τ, constants) - lambda_r_11(δ, τ, constants))^2
    denominator = (1.0 + 2.0 * lambda_r_01(δ, τ, constants) + lambda_r_02(δ, τ, constants))

    fraction = numerator / denominator

    return cv_term + (R_univ * fraction)
end

function k_T(T::AbstractFloat, rho::AbstractFloat, constants::constants_EoS_N2)
    (; T_c, rho_c, R_univ) = constants
    τ = T_c / T
    δ = rho / rho_c
    term1 = rho * R_univ * T
    term2 = 1.0 + 2.0 * lambda_r_01(δ, τ, constants) + lambda_r_02(δ, τ, constants)
    return 1.0 / (term1 * term2)
end




function cv0_calc(T::AbstractFloat, constants::constants_EoS_N2)
    (; T_c, R_univ) = constants
    τ = T_c / T
    
    return - R_univ * lambda_0_20(1.0, τ, constants)
end

function u0_calc(T::AbstractFloat, constants::constants_EoS_N2)
    (; T_c, R_univ) = constants
    τ = T_c / T
    
    return R_univ * T * lambda_0_10(1.0, τ, constants)
end

function h0_calc(T::AbstractFloat, constants::constants_EoS_N2)
    (; R_univ) = constants
    return u0_calc(T, constants) + R_univ * T
end

function internal_energy_calc(T::AbstractFloat, δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; R_univ) = constants
    u_ideal = u0_calc(T, constants)
    
    u_residual = R_univ * T * lambda_r_10(δ, τ, constants)

    return u_ideal + u_residual
end

function enthalpy_calc(T::AbstractFloat, δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; R_univ) = constants
    h_ideal = h0_calc(T, constants)
    
    h_residual = R_univ * T * (lambda_r_10(δ, τ, constants) + lambda_r_01(δ, τ, constants))
    return h_ideal + h_residual
end

function entropy_calc(δ::AbstractFloat, τ::AbstractFloat, constants::constants_EoS_N2)
    (; R_univ) = constants
    term1 = lambda_total_10(δ, τ, constants)
    term2 = alpha_0(δ, τ, constants) + alpha_r(δ, τ, constants)
    
    return R_univ * (term1 - term2)
end





function gibbs_free_energy(T::AbstractFloat, rho::AbstractFloat, constants::constants_EoS_N2)
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

function vapour_pressure_ancillary(T::AbstractFloat, constants::constants_EoS_N2)
    (; T_c) = constants
    
    # Critical pressure for Nitrogen from the paper (in Pa)
    p_c = 3.3958e6 

    # Coefficients and exponents from Eqn. 4
    N = [-6.12445284, 1.26327220, -0.765910082, -1.77570564]
    k = [1.0, 1.5, 2.5, 5.0]
    
    θ = 1.0 - T / T_c
    
    sum_val = 0.0
    for i in eachindex(N)
        sum_val += N[i] * (θ^k[i])
    end
    
    ln_pr = (T_c / T) * sum_val
    
    return exp(ln_pr) * p_c
end


function dPsat_dT(T::AbstractFloat, p_pair::AbstractFloat, constants::constants_EoS_N2)
    (; T_c) = constants
    
    # Critical pressure for Nitrogen from the paper (in Pa)
    p_c = 3.3958e6 

    # Coefficients and exponents from Eqn. 4
    N = [-6.12445284, 1.26327220, -0.765910082, -1.77570564]
    k = [1.0, 1.5, 2.5, 5.0]
    
    θ = 1.0 - T / T_c
    
    S_T = 0.0
    sum_dS_dT_term = 0.0
    for i in eachindex(N)
        S_T += N[i] * (θ^k[i])
        sum_dS_dT_term += N[i] * k[i] * (θ^(k[i] - 1.0))
    end

    derivative_term = ( (-T_c / (T^2.0)) * S_T ) - ( (1.0 / T) * sum_dS_dT_term )
    
    return p_pair * derivative_term # apply chain rule
end


function find_saturation_temperature(P_target::AbstractFloat, constants::constants_EoS_N2; max_iter=20, tol=1.0e-7)
    # Newton-Raphson method, similar to density

    (; T_c) = constants
    p_c = 3.3958e6 # Critical pressure in Pa (NITROGEN)
    T_triple = 63.151 # Triple point temperature in K (PARAHYDROGEN)
    p_triple = 12.528e3 # Triple point pressure in Pa (PARAHYDROGEN)

    # Check if there is such thing as saturation temperature
    if P_target > p_c
        println("Pressure ($P_target Pa) is above the critical pressure ($p_c Pa). Saturation temperature is not defined!!!")
        return 0.0
    end
    if P_target < p_triple
         println("Pressure ($P_target Pa) is below the triple point pressure ($p_triple Pa). Saturation temperature is not defined11!")
        return 0.0
    end



    # Initial guess for temp using a simplified inversion of the ancillary equation
    N1 = -6.12445284 # First coefficient from the ancillary equation
    T_guess = T_c / (1.0 + log(P_target / p_c) / N1)
    
    T = clamp(T_guess, T_triple, T_c - 1e-6) # Make sure guess is physical


    

    for it in 1:max_iter
        # Calculate residual
        P_calc = vapour_pressure_ancillary(T, constants)
        f = P_calc - P_target

        # Check for convergence
        if abs(f / P_target) < tol
            return T
        end

        # Calculate derivative for Newton-Raphson step
        dP_dT = dPsat_dT(T, P_calc, constants)
        
        if abs(dP_dT) < 1e-9 # Avoid division by zero
             error("T_SAT : DERIVATIVE IS ALMOST ZERO")
        end

        # Newton-Raphson step
        T_new = T - f / dP_dT
        
        T = clamp(T_new, T_triple, T_c - 1e-6)
    end
    
    error("T_saturation solver did not converge in $max_iter iterations.")
end



function find_density_advanced(T::AbstractFloat, P_target::AbstractFloat, rho_guess::AbstractFloat, constants::constants_EoS_N2;
    max_iter=100, tol=1.0e-8) # DOES NOT SUPPORT GPU
    (; T_c, rho_c, R_univ, M_N2, T_ref, a_nitro, N, t, d, p, α, β, γ, D) = constants

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



function find_saturation_properties(T::AbstractFloat, pressure::AbstractFloat, constants::constants_EoS_N2; max_iter=20, tol=1.0e-7)
    (; T_c, rho_c, R_univ) = constants
    if T >= T_c # Ensure the temperature is in the valid range (below critical)
        error("Temperature must be below the critical temperature for saturation calculation.")
    end

    p_guess = vapour_pressure_ancillary(T, constants) # This function would vary depending on fluid

    # Secant method for iterative solver
    p_sat = p_guess
    p_prev = p_guess * 0.999 # Slightly perturb the previous pressure value for first step

    g_diff_current = 0.0
    g_diff_prev = 0.0

    for it in 1:max_iter
        p_ideal_gas = p_sat / (R_univ * T)
        p_multiplied =  5.0 * rho_c

        rho_l = find_density_advanced(T, p_sat, p_multiplied, constants) # Higher guess for liquid
        rho_v = find_density_advanced(T, p_sat, p_ideal_gas, constants) # Ideal gas guess for vapour

        # Check if roots were found
        if isnan(rho_l) || isnan(rho_v)
            error("Failed to find liquid or vapor density root at T=$T K, P=$p_sat Pa...")
        end

        # Compute difference in gibbs free energy
        g_l = gibbs_free_energy(T, rho_l, constants)
        g_v = gibbs_free_energy(T, rho_v, constants)
        g_diff_current = g_l - g_v

        # Convergence check based on gibbs energy
        if abs(g_diff_current / g_l) < tol
            # println("Saturation solver converged in $it iterations.")


            T_sat = 0.0
            T_sat = find_saturation_temperature(pressure, constants::constants_EoS_N2)

            return (p_sat, T_sat, rho_l, rho_v) # Conversion into appropriate units
        end

        # Secant method step
        if it > 1
            if abs(g_diff_current - g_diff_prev) < 1e-9 # Avoid division by zero
                p_next = p_sat * 1.001
            else
                p_next = p_sat - g_diff_current * (p_sat - p_prev) / (g_diff_current - g_diff_prev)
            end
            
            # Update values for next itration
            p_prev = p_sat
            g_diff_prev = g_diff_current
            p_sat = p_next
        else
            
            p_prev = p_sat
            g_diff_prev = g_diff_current
            
            p_sat = p_guess * 1.001
        end

        if p_sat <= 0 # Check for -ve temps
            error("Solver stepped to a -ve pressure.")
        end
    end

    error("Saturation solver did not converge in $max_iter iterations.")
end




function EOS_wrapper_N2(T::AbstractFloat, pressure::AbstractFloat, alpha::AbstractFloat, c_τ_evap::AbstractFloat=0.0, c_τ_cond::AbstractFloat=0.0)
    constants = constants_EoS_N2(
        126.192,  # T_c
        11.1839e3,  # rho_c, multiplied by e3 for convenience
        28.01348e-3, # M_N2, multiplied by e-3 for convenience
        8.314472, # R_univ
        252.384, # T_ref = T_c * 2
        [
            2.5, -12.76952708, -0.00784163, -1.934819e-4, -1.247742e-5,
            6.678326e-8, 1.012941, 26.65788
        ], # a_nitro
        [
            # Polynomial Part (k=1 to 6)
            0.924803575275, -0.492448489428, 0.661883336938, -0.192902649201e1, 
            -0.622469309629e-1, 0.349943957581,
            # Exponential Part (k=7 to 32)
            0.564857472498, -0.161720005987e1, -0.481395031883, 0.421150636384,
            -0.161962230825e-1, 0.172100994165, 0.735448924933e-2, 0.168077305479e-1,
            -0.107626664179e-2, -0.137318088513e-1, 0.635466899859e-3, 0.304432279419e-2,
            -0.4357623366045e-1, -0.723174889316e-1, 0.389644315272e-1, -0.212201363910e-1,
            0.40882298181509e-2, -0.551990017984e-4, -0.462016716479e-1, -0.300311716011e-2,
            0.368825891208e-1, -0.255856846220e-2, 0.896915264558e-2, -0.441513370350e-2,
            0.133722924858e-2, 0.264832491957e-3,
            # Gaussian Part (k=33 to 36)
            0.19668894015e2, -0.209115600730e2, 0.1677883066989e-1, 0.2627675665274e4
        ], # N
        [
            1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 1.0, 1.0, 1.0, 3.0, 3.0, 4.0, 6.0, 
            6.0, 7.0, 7.0, 8.0, 8.0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 4.0, 5.0, 
            5.0, 8.0, 3.0, 5.0, 6.0, 9.0, 1.0, 1.0, 3.0, 2.0
        ], # d
        [
            0.25, 0.875, 0.5, 0.875, 0.375, 0.75, 0.5, 0.75, 2.0, 1.25, 3.5, 
            1.0, 0.5, 3.0, 0.0, 2.75, 0.75, 2.5, 4.0, 6.0, 6.0, 3.0, 3.0, 6.0, 
            16.0, 11.0, 15.0, 12.0, 12.0, 7.0, 4.0, 16.0, 0.0, 1.0, 2.0, 3.0
        ], # t
        [
            0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
            2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 0, 0, 0, 0
        ], # p
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 20.0, 15.0, 25.0
        ], # α
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 325.0, 325.0, 300.0, 275.0
        ], # β
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.16, 1.16, 1.13, 1.25
        ], # γ
        [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ]  # D
    )
    
    (; T_c, rho_c, R_univ, M_N2, T_ref) = constants

    pressure_tol = 1e-4
    alpha_tol = 1e-8

    rho_mol = 0.0
    rho_list = [0.0, 0.0]

    latentHeat = 0.0
    T_sat = 0.0

    m_lv = 0.0
    m_vl = 0.0
    
    
    # Firstly, account for supercritical fluid state
    if T >= T_c # Maybe add fancier check "within tolerance"
        rho_guess = pressure / (R_univ * T) # Ideal gas guess
        rho_mol = find_density_advanced(T, pressure, rho_guess, constants)

    else # Else it is not supercritical
        (P_sat, T_sat, rho_l_sat, rho_v_sat) = find_saturation_properties(T, pressure, constants)

        # println("rho_l : $rho_l_sat, rho_v: $rho_v_sat")

        if ( abs(pressure - P_sat) / P_sat ) < pressure_tol # TWO PHASE REGION, pressure matched saturation line
            # error("TWO PHASE REGION AT T=$T !")
            # println("WE ARE HERE")
            # isMultiphase = true

            rho_mol_liquid = find_density_advanced(T, pressure, rho_l_sat, constants)
            rho_mol_vapour = find_density_advanced(T, pressure, rho_v_sat, constants)

            # println("rho_mol_LIQUID: $rho_mol_liquid, rho_mol_VAPOUR: $rho_mol_vapour")

            rho_list[1] = rho_mol_liquid
            rho_list[2] = rho_mol_vapour

            # compute everything for both phases and use alpha to find mixtures?


        # elseif (alpha <= (1.0 - alpha_tol) && alpha >= (0.0 + alpha_tol)) # TWO PHASE REGION, based on alpha
        #     # need to also account for mass flow rate!!!
        #     isMultiphase = true

        #     rho_mol_liquid = find_density_advanced(T, pressure, rho_l_sat, constants)
        #     rho_mol_vapour = find_density_advanced(T, pressure, rho_v_sat, constants)

        #     rho_list[1] = rho_mol_liquid
        #     rho_list[2] = rho_mol_vapour

        elseif pressure < P_sat # VAPOUR REGION
            # println("VAPOUR")
            rho_guess = pressure / (R_univ * T) # Ideal gas guess
            rho_guess_ = 5.0 * rho_c # Higher density guess
            rho_mol = find_density_advanced(T, pressure, rho_guess, constants)
            rho_liquid = find_density_advanced(T_sat, pressure, rho_guess_, constants)

            rho_list[1] = rho_liquid
            rho_list[2] = rho_mol

            m_lv = c_τ_evap * alpha * rho_liquid * ( (T - T_sat)/T_sat ) # ASSUME alpha=1 is liquid

        elseif pressure > P_sat # LIQUID REGION
            # println("LIQUID")
            rho_guess = 5.0 * rho_c # Higher density guess
            rho_guess_ = pressure / (R_univ * T) # Ideal gas guess
            rho_mol = find_density_advanced(T, pressure, rho_guess, constants)
            rho_vapour = find_density_advanced(T_sat, pressure, rho_guess_, constants)

            rho_list[1] = rho_mol
            rho_list[2] = rho_vapour

            m_vl = c_τ_cond * (1.0-alpha) * rho_vapour * ( (T_sat - T)/T_sat )# ASSUME alpha=0 is liquid vapour
        end

        if isnan(rho_mol)
            error("Failed to find a valid density for the given single-phase state.")
        end
    end

    if T >= T_c
        # DO NOT CARE ABOUT THIS FOR NOW. SHOULD BE SINGLE PHASE THOUGH!
        # rho_val, cp_val, cv_val, kT_val, kT_ref_val, 
        # internal_energy_val, enthalpy_val, entropy_val = params_computation(rho_mol, T, constants)

        # return rho_val, cp_val, cv_val, kT_val, kT_ref_val, 
        #         internal_energy_val, enthalpy_val, entropy_val, latentHeat, T_sat
    else
        rho_vals = [0.0, 0.0]
        cv_vals = [0.0, 0.0]
        cp_vals = [0.0, 0.0]
        kT_vals = [0.0, 0.0]
        kT_ref_vals = [0.0, 0.0]
        internal_energy_vals = [0.0, 0.0]
        enthalpy_vals = [0.0, 0.0]
        entropy_vals = [0.0, 0.0]

        for i in eachindex(rho_list)
            rho_vals[i], cp_vals[i], cv_vals[i], kT_vals[i], kT_ref_vals[i], 
            internal_energy_vals[i], enthalpy_vals[i], entropy_vals[i] = params_computation(rho_list[i], T, constants)
        end

        latentHeat = enthalpy_vals[2] - enthalpy_vals[1] #enthalpy_V - enthalpy_L
        

        return rho_vals, cp_vals, cv_vals, kT_vals, 
                kT_ref_vals, internal_energy_vals, enthalpy_vals, entropy_vals, latentHeat, T_sat, m_lv, m_vl

    end
end

function params_computation(rho_mol::AbstractFloat, T::AbstractFloat, constants::constants_EoS_N2)
    (; T_c, rho_c, R_univ, M_N2, T_ref) = constants

    τ = T_c / T
    δ = rho_mol / rho_c

    cv_mol = c_v(δ, τ, constants)
    cp_mol = c_p(δ, τ, constants)
    kT = k_T(T, rho_mol, constants)
    kT_ref = k_T(T_ref, rho_mol, constants)

    internal_energy_mol = internal_energy_calc(T, δ, τ, constants)
    enthalpy_mol = enthalpy_calc(T, δ, τ, constants)
    entropy_mol = entropy_calc(δ, τ, constants)

    conversion_factor = 1.0 / (M_N2*1.0e3)

    rho = rho_mol * M_N2

    cv = cv_mol*conversion_factor
    cp = cp_mol*conversion_factor
    internal_energy = internal_energy_mol*conversion_factor
    enthalpy = enthalpy_mol*conversion_factor
    entropy = entropy_mol*conversion_factor

    return rho, cv, cp, kT, kT_ref, internal_energy, enthalpy, entropy
end



function u_into_T(u_target::AbstractFloat, P_target::AbstractFloat, T_guess::AbstractFloat, constants::constants_EoS_N2;
                    max_iter=20, tol=1.0e-7)

    (; T_c, rho_c) = constants
    # u_target is u_n+1
    # p_target is p_n+1
    # T_guess is T_n

    T = T_guess # Start with an initial guess
    
    for it in 1:max_iter

        rho_current = find_density(T, P_target, constants) # A MORE ADVANCED APPROACH IS REQUIRED!

        if isnan(rho_current)
            error("Density solver failed during T inversion.")
        end

        τ = T_c / T
        δ = rho_current / rho_c #mol value!

        u_calc = internal_energy_calc(T, δ, τ, constants) #WARNING : need to check units!!!

        # Calculate the residual
        f = u_calc - u_target
        if abs(f / u_target) < tol
            return T # Converged
        end

        # Calculate the derivative (cv)
        cv_val = c_v(T, rho_current, constants)
        if abs(cv_val) < 1e-9
            error("Derivative (cv) is near zero.")
        end

        # Newton's step
        T = T - f / cv_val
    end

    error("Temperature inversion failed to converge.")
end



# T_input = 73.151 # Temperature in K
# P_input = 0.1e6 # Pressure in kPa


# rho0, cv0, cp0, kT0, internal_energy0, enthalpy0, entropy0 = EOS_wrapper_N2(T_input, P_input)


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