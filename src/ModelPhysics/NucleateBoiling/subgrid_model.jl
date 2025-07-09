export boiling_flux


g = 9.81

function boiling_flux(T_w::Float64, T_sat::Float64, sigma::Float64, rho_l::Float64, 
    rho_v::Float64, c_pl::Float64, k_l::Float64, h_lv::Float64, q_con::Float64, alpha_v::Float64)
    
    Npp = nucleate_site_density(T_w, T_sat)

    D_dep = bubble_departure_diameter(sigma, rho_l, rho_v)

    f = bubble_departure_frequency(rho_l, rho_v, D_dep)

    t_wait = t_w(f)

    Kq = K_q(2.0, D_dep, Npp) # First arg is F_A = 2

    h_q = h_q(f, t_wait, rho_l, c_pl, k_l, Kq)

    qq = q_q(h_q, T_w, T_sat)
    qe = q_e(Npp, f, D_dep, rho_v, h_lv)

    q_subgrid = q_sgnb(alpha_v, qq, qe)
    q_total = q_t(q_con, alpha_v, qq, qe)

    return q_subgrid, q_total
end

function nucleate_site_density(T_w::Float64, T_sat::Float64)
    m = 185.0
    n = 1.805

    return (m * abs(T_w - T_sat))^n
end

function bubble_departure_diameter(theta_c::Float64=0.772, sigma::Float64, rho_l::Float64, rho_v::Float64)
    # Theta_c = contact angle, we define it as 0.772 rad here
    D_cal = 1.5126e-2

    term1 = sqrt(sigma/(g*(rho_l - rho_v)))
    term2 = ((rho_l - rho_v)/rho_v)^0.9
    return D_cal * theta_c * term1 * term2
end

function bubble_departure_frequency(rho_l::Float64, rho_v::Float64, D_dep::Float64)
    return sqrt(4*g*(rho_l - rho_v)/(3 * D_dep * rho_l))
end

function h_q(f::Float64, t_w::Float64, rho_l::Float64, c_pl::Float64, k_l::Float64, K_q::Float64)
    #c_pl is the constant pressure heat capacity of the liquid
    return ( 2 * K_q * f ) * sqrt((rho_l * c_pl * k_l * t_w) / pi)
end

function q_e(Npp::Float64, f::Float64, D_dep::Float64, rho_v::Float64, h_lv::Float64) #qe''
    return Npp * f * ((pi*D_dep^3)/6) * rho_v * h_lv # how to find h_lv??? - refer to NIST REFPROP
end

function q_q(h_q::Float64, T_w::Float64, T_q::Float64) #qq''
    return h_q * (T_w - T_q) # T_q is usually assumed to equal T_sat
end

function K_q(F_A::Float64, D_dep::Float64, Npp::Float64) #F_A = 2 from article
    return ( F_A * pi * (D_dep^2) * Npp ) / 4
end

function t_w(f::Float64) #waiting time
    c_w = 0.8 #waiting coefficient
    return c_w / f
end


function q_t(q_con::Float64, alpha_v::Float64, q_q::Float64, q_e::Float64) # return this I guess??
    return q_con + (1 - alpha_v)*(q_q + q_e) #q_con calculated from CHT
end


function q_sgnb(alpha_v::Float64, q_q::Float64, q_e::Float64)
    return (1 - alpha_v)*(q_q + q_e)
end


# CHT coupling functions are also needed... e.g. T_w,eff