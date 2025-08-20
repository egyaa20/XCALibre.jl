export multiphase!


function multiphase!(
    model, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2)


    residuals = setup_multiphase_solvers(
        MULTIPHASE, model, config; 
        output=output,
        pref=pref,
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops
        )
        
    return residuals
end


function setup_multiphase_solvers(
    solver_variant, model, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    ) 

    (; solvers, schemes, runtime, hardware, boundaries) = config

    @info "Extracting configuration and input fields..."

    (; U, p, Uf, pf) = model.momentum
    (; rho, alpha, alphaf) = model.fluid

    phases = model.fluid.phases

    mesh = model.domain

    @info "Pre-allocating fields..."
    
    TF = _get_float(mesh)
    time = zero(TF) # assuming time=0


    @info "Computing Fluid Properties..."


    ∇alpha = Grad{schemes.alpha.gradient}(alpha)  
    n_f = FaceVectorField(mesh)



    ∇p = Grad{schemes.p.gradient}(p)
    grad!(∇p, pf, p, boundaries.p, time, config)
    limit_gradient!(schemes.p.limiter, ∇p, p, config)

    phif = FaceScalarField(mesh)
    phiC = FaceScalarField(mesh)
    U_c_f = FaceScalarField(mesh)

    rhof = FaceScalarField(mesh)
    mdotf = FaceScalarField(mesh)
    nueff = FaceScalarField(mesh)

    alpha_source = ScalarField(mesh)


    interpolate!(alphaf, alpha, config)

    flux!(phif, Uf, config)



    # Umprev = VectorField(mesh)
    # initialise!(Umprev, [0.0, 0.0, 3.0])
    # a_test = construct_acceleration_field(U, Umprev, model.fluid.physics_properties, 0.5, config.hardware.workgroup, mesh)



    rho_fractions = [ScalarField(mesh) for _ in phases]
    nu_fractions  = [ScalarField(mesh) for _ in phases]
    nuf_fractions  = [FaceScalarField(mesh) for _ in phases]

    
    # DUMMY FIELD

    T_temp = ScalarField(mesh)
    initialise!(T_temp, 70.0)
    
    # DUMMY FIELD

    
    if typeof(model.fluid.phases[1].eos) <: HelmholtzEnergy
        helmholtz_energy!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, model, config)
    else
        update_properties!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, phases, model, config)
    end


    # DRIFT VELOCITY LOOP (check if driftVelocity is activated in sources)

    U_prev = VectorField(mesh)
    U_prev = U

    ∇U = Grad{schemes.U.gradient}(U)
    grad!(∇U, Uf, U, boundaries.U, time, config)
    limit_gradient!(schemes.U.limiter, ∇U, U, config)

    a_field = construct_acceleration_field(U, ∇U, U_prev, dt)
    isInitialisation = true

    v_pq = VectorField(mesh)
    v_pq_prev = VectorField(mesh)

    v_dr_p = VectorField(mesh)
    v_dr_q = VectorField(mesh)
    v_p = VectorField(mesh)
    v_q = VectorField(mesh)

    d_p = ScalarField(mesh)
    initialise!(d_p, model.fluid.physics_properties.driftVelocity.d_p)


    # LOOK UP U FIELD AND RETURN MAX VALUE IN EACH AXIS THEN USE IT IN HERE
    # LOOK UP VELOCITY BCs AND RETURN MAX VALUE
    
    find_Vpq!(rho_fractions[2], rho_fractions[1], rho, nu_fractions[2], nu_fractions[1], d_p, a_field, v_pq_prev, isInitialisation, SVector(1.0e-9, 1.0e-9, 1.0e-9), SVector(10.0, 10.0, 10.0), v_pq, mesh)

    update_velocities!(U, v_pq, v_dr_p, v_dr_q, v_p, v_q, rho, rho_fractions[2], rho_fractions[1], alpha, mesh) # liquid <> p ; vapour <> q


    # DRIFT VELOCITY LOOP



    flux!(mdotf, Uf, rhof, config)


    grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)

    interpolate!(n_f, ∇alpha.result, config)


    println(∇alpha.result) # vector field 
    println(n_f) #face vector field
    # correct_boundaries!(n_f, ∇alpha.result, boundaries.alpha, time, config) # PROBLEMATIC LINE!



    cAlpha = 1.0

    initialise!(U_c_f, 1.0)
    initialise!(phiC, 1.0)
    # @. U_c_f.values *= n_f.values
    @. U_c_f.values *= cAlpha

    flux!(U_c_f, Uf, config)

    alphaf_minus = FaceScalarField(mesh)
    initialise!(alphaf_minus, 1.0)
    @. alphaf_minus.values -= alphaf.values

    @. phiC.values *= U_c_f.values
    @. phiC.values *= alphaf.values
    @. phiC.values *= alphaf_minus.values

    div!(alpha_source, phiC, config)

    # update_nueff!(nueff, nu, model.turbulence, config) #need to work with nu




    psi = ScalarField(mesh)
    initialise!(psi, 1.0)

    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)

    divHv = ScalarField(mesh)
    div!(divHv, mdotf, config)


    @info "Defining models..."

    zero_field = ScalarField(mesh)

    
    # Alpha_p... ? rho_p ?
    alpha_eqn = ( # Volume Fraction Transport Eqn
        Time{schemes.alpha.time}(alpha)
        + Divergence{schemes.alpha.divergence}(phif, alpha) #phif is U * n * alpha 
        # divergence multiplied by rho_l ????
        ==
        - Source(zero_field) #zero_field
        # Lee Source
        # Subgrid source
    ) → ScalarEquation(alpha, boundaries.alpha)

    # println(phif.values)

    props = model.fluid.physics_properties

    momentum_rhs = - Source(∇p.result)
    momentum_rhs = momentum_eqn_sources(props, momentum_rhs, rho, mesh)



    U_eqn = ( # Momentum Eqn
        Time{schemes.U.time}(rho, U) # rho * U
        + Divergence{schemes.U.divergence}(mdotf, U)  #mdotf is rho*U*n
        - Laplacian{schemes.U.laplacian}(nueff, U) # nu_mixture; nu needs update!!!
        == momentum_rhs
        # - Source(∇p.result) # -∇p
        # + Source(rhoG) # (0, -9.81 * rho, 0)
        # + Source(Su) # nucleate boiling subgrid model
        # - Source(slipVelocity) # slip velocity term
    ) → VectorEquation(U, boundaries.U)

    
    g = -9.81
    prgh = ScalarField(mesh)
    ∇prgh = Grad{schemes.p.gradient}(prgh)

    p_eqn = (
       Time{schemes.p.time}(psi, p)
      - Laplacian{schemes.p.laplacian}(rDf, prgh)
      ==
      - Source(divHv)
    ) → ScalarEquation(p, boundaries.p)

    # @. prgh.values = p.values - rho.values * g # * yCoords
    # grad!(∇prgh, pf, prgh, boundaries.p, time, config)
    # limit_gradient!(schemes.p.limiter, ∇prgh, prgh, config)


    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."
     
    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))


    
    if typeof(model.energy) <: MultiphaseEnergy
        @info "Initialising energy model..."
        energyModel = initialise(model.energy, model, mdotf, rho, p_eqn, config) #change
    end

    @info "Initialising turbulence model..."
    turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        model, turbulenceModel, alpha_eqn, ∇p, U_eqn, p_eqn, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        time, rDf, nueff, phif, prgh, ∇prgh, mdotf, divHv
        )

    return residuals    
end # end function






# It might make sense to just dispatch HelmholtzEnergy case
function MULTIPHASE(
    model, turbulenceModel, alpha_eqn, ∇p, U_eqn, p_eqn, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2,
    time, rDf, nueff, phif, prgh, ∇prgh, mdotf, divHv
    )

    inner_loops = 2
    

    g = 9.81

    (; U, p, Uf, pf) = model.momentum
    # (; nu) = model.fluid
    (; rho, rhof, alpha, alphaf) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries) = config
    (; iterations, write_interval, dt) = runtime
    (; backend) = hardware

    phases = model.fluid.phases

    rho_fractions = [ScalarField(mesh) for _ in phases] #need to pass them in args
    nu_fractions  = [ScalarField(mesh) for _ in phases]
    nuf_fractions = [FaceScalarField(mesh) for _ in phases]


    Hv = VectorField(mesh)


    rD = ScalarField(mesh)
    outputWriter = initialise_writer(output, model.domain)
    
    TF = _get_float(mesh)

    n_cells = length(mesh.cells)
    prev = KernelAbstractions.zeros(backend, TF, n_cells) 

    # Pre-allocate vectors to hold residuals 
    R_ux = ones(TF, iterations)
    R_uy = ones(TF, iterations)
    R_uz = ones(TF, iterations)
    R_p = ones(TF, iterations)
    R_alpha = ones(TF, iterations)
    cellsCourant =adapt(backend, zeros(TF, length(mesh.cells)))
    

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting PISO loops..."

    progress = Progress(iterations; dt=1.0, showspeed=true)

    @time for iteration ∈ 1:iterations
        time = iteration *dt
        
        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config; time=time)

        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config)
        interpolate!(rDf, rD, config)
        remove_pressure_source!(U_eqn, ∇p, config)
        
        rp = 0.0
        # println(inner_loops)
        for i ∈ 1:inner_loops
            H!(Hv, U, U_eqn, config)
            
            # Interpolate faces
            interpolate!(Uf, Hv, config) # Careful: reusing Uf for interpolation
            correct_boundaries!(Uf, Hv, boundaries.U, time, config)
            # div!(divHv, Uf, config)

            # new approach
            flux!(mdotf, Uf, config)
            div!(divHv, mdotf, config)
            
            # Pressure calculations (previous implementation)
            @. prev = p.values
            rp = solve_equation!(p_eqn, p, boundaries.p, solvers.p, config; ref=pref, time=time)
            if i == inner_loops
                explicit_relaxation!(p, prev, 1.0, config)
            else
                explicit_relaxation!(p, prev, solvers.p.relax, config)
            end

            grad!(∇p, pf, p, boundaries.p, time, config) 
            limit_gradient!(schemes.p.limiter, ∇p, p, config)

            # new approach
            interpolate!(Uf, U, config) # velocity from momentum equation

            correct_boundaries!(Uf, U, boundaries.U, time, config)
            flux!(mdotf, Uf, config)
            correct_mass_flux(mdotf, p, rDf, config)
            correct_velocity!(U, Hv, ∇p, rD, config)
        end

        
    if !(typeof(model.turbulence) <: Laminar) #isa
        println("TURBULENCE ON")
        turbulence!(turbulenceModel, model, S, prev, time, config) 
    end


    update_nueff!(nueff, nueff, model.turbulence, config)
    

    flux!(phif, Uf, config)
    
    ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config; time=time)

    @. alpha.values = clamp(alpha.values, 0.0, 1.0)

    interpolate!(alphaf, alpha, config)
    # flux!(phif, Uf, config)


    # DUMMY FIELD

    T_temp = ScalarField(mesh)
    initialise!(T_temp, 70.0)

    # DUMMY FIELD

    #update_thermo_properties!(phase1 phase2 etc)
    
    if typeof(model.fluid.phases[1].eos) <: HelmholtzEnergy #HelmholtzEnergyEOS
        helmholtz_energy!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, model, config)
    else
        update_properties!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, phases, model, config)
    end


    if typeof(model.energy) <: MultiphaseEnergy
        println("ENERGY ON")
        energy!(energyModel, model, prev, mdotf, rho, mueff, time, config)
    end

    maxCourant = max_courant_number!(cellsCourant, model, config)

    R_ux[iteration] = rx
    R_uy[iteration] = ry
    R_uz[iteration] = rz
    R_p[iteration] = rp
    R_alpha[iteration] = ralpha

    ProgressMeter.next!(
        progress, showvalues = [
            (:time, iteration*runtime.dt),
            (:Courant, maxCourant),
            (:Ux, R_ux[iteration]),
            (:Uy, R_uy[iteration]),
            (:Uz, R_uz[iteration]),
            (:p, R_p[iteration]),
            (:alpha, R_alpha[iteration]),
            turbulenceModel.state.residuals...
            ]
        )

    if iteration%write_interval + signbit(write_interval) == 0
        
        save_output(model, outputWriter, iteration, time, config)
    end

    end # end for loop

    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p, alpha=R_alpha)
end





function momentum_eqn_sources(props, rhs, rho, mesh)
    if ( hasproperty(props, :gravity) ) && ( props.gravity isa Gravity )
        g = props.gravity.g

        x0, y0, z0 = g[1], g[2], g[3]
        # x = ConstantScalar(x0)
        # y = ConstantScalar(y0)
        # z = ConstantScalar(z0)
        # CAN VECTOR FIELD BE DISPATCHED TO INCLUDE CONSTANT SCALARS ???

        x = ScalarField(mesh)
        y = ScalarField(mesh)
        z = ScalarField(mesh)
        initialise!(x, x0)
        initialise!(y, y0)
        initialise!(z, z0)

        @. x.values = x.values * rho.values
        @. y.values = y.values * rho.values
        @. z.values = z.values * rho.values


        rhoG = VectorField(x, y, z, mesh)

        @info "Adding Gravity term..."
        rhs -= Source(rhoG)
    end

    return rhs
end

function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0.values * alpha_field.values) + (property_1.values * (1.0 - alpha_field.values))
    nothing
end




function update_properties!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, phases, model, config)
    (; p) = model.momentum
    (; rho, rhof, alpha, alphaf) = model.fluid

    for phase in eachindex(phases)
        current_phase = phases[phase]
        rho_ = rho_fractions[phase]
        nu_ = nu_fractions[phase]
        nuf_ = nuf_fractions[phase]

        if current_phase.eos isa ConstEos
            initialise!(rho_, current_phase.eos.rho)
        else current_phase.eos isa PerfectGas
            R = current_phase.eos.R
            perfectGas!(rho_, 273.0, R, p)
        end

        
        if current_phase.mu isa ConstMu
            mu0 = current_phase.mu.mu
            @. nu_.values = mu0 / rho_.values
        elseif current_phase.mu isa Sutherland # CHECK FOR ENERGY MODEL!!! T FIELD IS REQUIRED
            mu_ref = current_phase.mu.mu_ref
            S = current_phase.mu.S
            sutherland!(nu_, T_temp, mu_ref, S) # mu is computed
            @. nu_.values *= 1.0/rho_.values # conversion
        elseif current_phase.mu isa Andrade # CHECK FOR ENERGY MODEL!!! T FIELD IS REQUIRED
            B = current_phase.mu.B
            C = current_phase.mu.C
            andrade!(nu_, T_temp, B, C) # mu is computed
            @. nu_.values *= 1.0/rho_.values # conversion
        end

        interpolate!(nuf_, nu_, config)
    end
    

    blend_properties!(rho, alpha, rho_fractions[1], rho_fractions[2]) # 1 is alpha=0,   2 is alpha=1
    blend_properties!(nueff, alphaf, nuf_fractions[1], nuf_fractions[2]) # weird mu <-> nu conversion


    interpolate!(rhof, rho, config)
end


function perfectGas!(density_field, T_ref, R, pressure_field)
    @. density_field.values = (pressure_field.values) / (R * T_ref) # CAREFUL WITH p=0 initialisation
    nothing
end

function sutherland!(mu_field, T_field, mu_ref, S)
    T_ref = 273.0
    @. mu_field.values = (mu_ref * (T_field.values/T_ref)^1.5) * ((T_ref + S)/(T_field.values + S))
    nothing
end

function andrade!(mu_field, T_field, B, C) # maybe do it in Phase(mu=(...)) instead?
    @. mu_field.values = B * exp(C / T_field.values)
    nothing
end



function construct_acceleration_field(U_m, ∇U_m, U_m_prev, dt) #, workgroup
    grad_U = ∇U_m.result
    a = VectorField(mesh)

    g = props.DriftVelocity.gravity
    # g = [0.0, -9.81, 0.0] #test

    x0, y0, z0 = g[1], g[2], g[3]

    x = ScalarField(mesh)
    y = ScalarField(mesh)
    z = ScalarField(mesh)
    initialise!(x, x0)
    initialise!(y, y0)
    initialise!(z, z0)
    G = VectorField(x, y, z, mesh)

    dUdt = VectorField(mesh)
    
    for i in eachindex(U_m)
        dUdt[i] = (U_m[i] - U_m_prev[i]) / dt
        a[i] = G[i] - (grad_U[i] * U_m[i]) - dUdt[i]
    end
    
    return a
end


function get_E!(E, rho_p, d_p, mu_p, rho, a_field)
    rho_p = rho_p.values
    d_p = d_p.values
    mu_p = mu_p.values
    rho = rho.values

    for i in eachindex(E)
        E[i] = (rho_p[i]*((d_p[i])^2))/(18*mu_p[i]) * ((rho_p[i]-rho[i])/rho_p[i]) * a_field[i]
    end

    return nothing
end

function compute_RE!(RE, rho_q, v_pq, d_p, mu_q)
    rho_q = rho_q.values
    d_p = d_p.values
    mu_q = mu_q.values
    RE = RE.values
    
    for i in eachindex(RE)
        RE[i] = (rho_q[i] * d_p[i] * norm(v_pq[i]))/mu_q[i]
    end

    return nothing
end

function find_Vpq!(rho_q, rho_p, rho, mu_q, mu_p, d_p, a_field, v_pq_prev, isInitialisation, v_low, v_high, v_pq, mesh)
    @. mu_q.values = mu_q.values * rho_q.values #conversion again..... NU -> MU
    @. mu_p.values = mu_p.values * rho_p.values #conversion again.....
    
    E = VectorField(mesh)
    expr = ScalarField(mesh)
    v_pq_trial = VectorField(mesh)
    RE_trial = ScalarField(mesh)

    get_E!(E, rho_p, d_p, mu_p, rho, a_field)

    @. expr.values = mu_q.values/(0.0183 * rho_q.values * d_p.values)
    
    for i in eachindex(v_pq_trial)
        v_pq_trial[i] = sqrt.(E[i]*expr.values[i])
    end

    compute_RE!(RE_trial, rho_q, v_pq_trial, d_p, mu_q)

    K = ScalarField(mesh)
    @. K.values = 0.15 * ((rho_q.values * d_p.values)/mu_q.values)^0.687

    for i in eachindex(RE_trial)
        if RE_trial[i] > 1000
            v_pq[i] = v_pq_trial[i]
        else
            if isInitialisation
                v_pq[i] = Vpq_bisect(K[i], E[i], v_low, v_high)
            else
                v_pq[i] = Vpq_newton_raphson(E[i], K[i], v_pq_prev[i])
            end
        end
    end
end



function Vpq_newton_raphson(E, K, v; max_iter=20, tol=1.0e-7)
    for it in 1:max_iter # E IS A VECTOR FIELD
        v_pow_0_687 = v .^ 0.687
        f = v .+ K .* v_pow_0_687 .* v .- E
        # f = v + K*(v^1.687) - E

        f_prime = 1.0 .+ 1.687 .* K .* v_pow_0_687
        # f_prime = 1 + 1.687 * K * (v^1.687)

        v_new = v .- f ./ f_prime # division by zero in one direction is possible
        # v_new = v - (f/f_prime)

        err = abs.((v_new .- v) ./ (v_new .+ eps())) # Add eps to prevent division by zero
        # err = abs((v_new - v) / v_new)

        if maximum(err) < tol
            v_pq = v_new
            println(v_pq)

            return v_pq
        end

        v = v_new
    end
end


function bisect_f(v, K, E)
    return v .+ K .* (v .^ 1.687) .- E
end

function Vpq_bisect(K, E, v_low, v_high; max_iter=50, tol=1.0e-7) #v_mid=SVector(0.0, 0.0, 0.0)
    for it in 1:max_iter
        v_mid = (v_low .+ v_high) ./ 2.0
        
        f_mid = bisect_f(v_mid, K, E)
        f_low = bisect_f(v_low, K, E)

        vectorised_check = sign.(f_mid) .== sign.(f_low)
        v_low = ifelse.(vectorised_check, v_mid, v_low)
        v_high = ifelse.(vectorised_check, v_high, v_mid)

        # if sign.(f_mid) == sign.(f_low)
        #     v_low = v_mid
        # else
        #     v_high = v_mid
        # end

        err = abs.(v_high .- v_low)

        if maximum(err) < tol
            v_pq = v_mid

            return v_pq
        end
    end
end

function update_velocities!(U_m, v_pq, v_dr_p, v_dr_q, v_p, v_q, rho, rho_q, rho_p, alpha, mesh)
    alpha = alpha.values
    rho = rho.values
    rho_q = rho_q.values
    rho_p = rho_p.values

    C_p = ScalarField(mesh)
    C_q = ScalarField(mesh)

    @. C_p.values = (alpha * rho_p) / rho          #liquid
    @. C_q.values = ((1.0-alpha) * rho_q) / rho    #vapour

    
    for i in eachindex(U_m)
        v_dr_p[i] = (1.0 + C_p.values[i]) * v_pq[i]
        v_dr_q[i] = -(1.0 + C_q.values[i]) * v_pq[i]

        v_p[i] = v_dr_p[i] - U_m[i]
        v_q[i] = v_dr_q[i] - U_m[i]
    end
end






function hydrogen_viscosity!() # maybe a functor?

end

function nitrogen_viscosity!() # maybe a functor?

end

function peng_robinson!(T, p, T_c, p_c, ω, M)
    # using Polynomials
    # ASK HUMBERTO !!!

    # T_c = 33.145
    # p_c = 1.2964e6
    R_univ = 8.314
    M = M * 1.0e-3
    # ω = -0.216 # Acentric factor for H2

    T_r = T / T_c

    a_H2 = 0.45724 * (((R_univ^2)*(T_c^2))/p_c)
    b_H2 = 0.07780 * ((R_univ*T_c)/p_c)

    κ = 0.37464 + 1.54226*ω - 0.26992*(ω^2)
    α = (1 + κ*(1 - sqrt(T_r)))^2
    
    A = ( a_H2 * α * p ) / ( (R_univ^2) * (T^2) )
    B = ( b_H2 * p)  / (R_univ * T)
    
    coeffs = [
        -(A*B - B^2 - B^3),   # constant term (c0)
         A - 2B - 3B^2,       # Z term (c1)
        -(1 - B),             # Z^2 term (c2)
         1.0                  # Z^3 term
    ]
    
    poly = Polynomial(coeffs)
    Z_roots = roots(poly)

    # Extract only real roots
    Z_real = [z for z in Z_roots if isreal(z)]
    Z_real = real.(Z_real)

    Z_liq = minimum(Z_real)   # liquid root
    Z_vap = maximum(Z_real)   # vapour root

    # Convert to molar volumes
    Vm_liq = Z_liq * R_univ * T / p
    Vm_vap = Z_vap * R_univ * T / p

    # Convert to densities [kg/m3]
    rho_liq = M / Vm_liq
    rho_vap = M / Vm_vap

    return Z_liq, Z_vap # liq = vap for single phase so its handled by alpha I suppose....
end


function helmholtz_energy!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, model, config)
    (; p) = model.momentum
    (; rho, rhof, alpha, alphaf) = model.fluid
    HelmholtzEOS = model.fluid.phases[1].eos

    # Horrible solution, replace later.....
    for cell in eachindex(rho)
        _, rho_temp, _, _, _, _, _, mu_temp, _, _, _, _, _, _ = HelmholtzEOS(T_temp[cell], p[cell], alpha[cell])
        
        for i in eachindex(rho_fractions)
            rho_fractions[i][cell] = rho_temp[i]
            nu_fractions[i][cell]  = mu_temp[i]
        end
    end

    @. nu_fractions[1].values *= nu_fractions[1].values/rho_fractions[1].values # conversion
    @. nu_fractions[2].values *= nu_fractions[2].values/rho_fractions[2].values # conversion
    interpolate!(nuf_fractions[1], nu_fractions[1], config)
    interpolate!(nuf_fractions[2], nu_fractions[2], config)
    
    blend_properties!(rho, alpha, rho_fractions[1], rho_fractions[2]) # 1 is alpha=0,   2 is alpha=1
    blend_properties!(nueff, alphaf, nuf_fractions[1], nuf_fractions[2]) # weird mu <-> nu conversion

    interpolate!(rhof, rho, config)
end
