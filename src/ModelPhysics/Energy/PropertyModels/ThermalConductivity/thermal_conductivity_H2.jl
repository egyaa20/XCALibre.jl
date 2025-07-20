
###Refer to "Correlation of the Thermal Conductivity of Normal and Parahydrogen
###                              from the Triple Point to 1000 K and up to 100MPa", 2011

using Printf

T_c   = 32.938
rho_c = 31.323

A1 = [-1.24500e0,  3.10212e2, -3.31004e2,  2.46016e2,
            -6.57810e1,  1.08260e1, -5.19659e-1,  1.43979e-2]
A2 = [ 1.42304e4, -1.93922e4,  1.58379e4, -4.81812e3,
             7.28639e2, -3.57365e1,  1.00000e0]

B1 = [ 2.65975e-2, -1.33826e-3,  1.30219e-2, -5.67678e-3, -9.23380e-5]
B2 = [-1.21727e-3,  3.66663e-3,  3.88715e-3, -9.21055e-3,  4.00723e-3]


C1, C2, C3 = 3.57e-4, -2.46e-2, 0.2


R_D   = 1.01
ν     = 0.63
γ_crit = 1.2415
ξ₀    = 1.5e-10
Γ₀    = 0.172
qD    = 1 / 5.0e-10
T_ref = 1.5 * T_c

k_B   = 1.380649e-23

function lambda0(T::Float64)
    T_r = T / T_c
    numerator = sum(A1[i+1] * (T_r)^i for i in 0:7)
    denominator = sum(A2[i+1] *(T_r)^i for i in 0:6)
    return numerator / denominator
end


function delta_lambda(rho::Float64, T::Float64)
    T_r = T / T_c
    rho_r = rho / rho_c

    term_sum = 0.0
    for i in 1:5
        term_sum += (B1[i] + (B2[i] * T_r) ) * (rho_r)^i
    end
    return term_sum
end


function delta_lambda_c_empirical(rho::Float64, T::Float64) #The easy version
    delta_T_c   = (T / T_c) - 1
    delta_rho_c   = (rho / rho_c) - 1
    denominator  = C2 + abs(delta_T_c)

    if denominator <= 0.0
        return 0.0
    else
        return (C1 / denominator) * exp(-(C3 * delta_rho_c)^2)
    end
end

function delta_lambda_c(rho::Float64, T::Float64) #The tricky one

end
