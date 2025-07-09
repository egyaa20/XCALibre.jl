# Note : values are for parahydrogen


# A 1,i table for hydrogen
# A 2,i table for hydrogen

# B 1,i table for hydrogen
# B 2,i table for hydrogen


# Hydrogen constants: T_c, rho_c, T_ref

# General constants: qD, Г, 

# Known values that are required: Cp, Cv, T, rho

T = 50
rho = 0.1

cp = 500
cv = 500


# 


T_c = 33
rho_c = 31.3
T_ref = (3/2) * T_c


q_bar_D = 2.0e9
q_bar_D_inversed = 5.0e-10
R_D = 1.01
nu = 0.63 #the one that looks like 'v'
gamma = 1.2415
GAMMA = 0.052
xi_0 = 1.5e-10

k_b = 1.380649e-23


# Initial calculations





# Compute dilute-gas term, k0

function k0(T) #not adapted for kernel
    local_fraction = T / T_c
    num = zero(local_fraction)
    den = zero(local_fraction)
    for i in 0:6
        num += A1[i+1] * local_fraction^i
    end
    for i in 0:3
        den += A2[i+1] * local_fraction^i
    end
    return num / den
end

# Compute the excess term, delta_k

function delta_k(rho, T) #not adapted for kernel
    delta = 0.0
    for i in 1:5
        delta += (B1[i] + B2[i] * (T / T_c)) * (rho / rho_c)^i
    end
    return delta
end

# Compute the critical enhancement term, delta_k_c

function delta_k_c(rho, T, Cp, Cv, mu, d_rho_dp)
    # correlation length xi
    #....
    return rho * Cp * R_D * k_B * T / (6π * mu * xi) * X_minus_X0 
end

# Compute k using three terms

function thermal_conductivity(rho, T; Cp, Cv, mu, d_rho_dp)
    return k0(T) + delta_k(rho, T) + delta_k_c(rho, T, Cp, Cv, mu, d_rho_dp)
end