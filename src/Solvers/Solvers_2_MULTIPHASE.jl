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
    



    rho_fractions = [ScalarField(mesh) for _ in phases]
    nu_fractions  = [ScalarField(mesh) for _ in phases]
    nuf_fractions  = [FaceScalarField(mesh) for _ in phases]

    
    # DUMMY FIELD

    T_temp = ScalarField(mesh)
    initialise!(T_temp, 100.0)
    
    # DUMMY FIELD

    
    if typeof(model.fluid.phases[1].eos) <: HelmholtzEnergy
        model.fluid.phases[1].eos(50.0, 1.0e5, 1.0)
        # helmholtz_properties!()
    else
        update_properties!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, phases, model, config)
    end



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
    initialise!(T_temp, 100.0)
    # DUMMY FIELD


    
    if typeof(model.fluid.phases[1].eos) <: HelmholtzEnergy
        model.fluid.phases[1].eos(50.0, 1.0e5, 1.0)
        # helmholtz_properties!()
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
    @. property_field.values = (property_0.values * alpha_field.values) + (property_1.values * (1 - alpha_field.values)) #cant do 1.0 (OR CAN I?)
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
    p_ref = 101_325.0 # for conversion to absolute pressure, avoid rho=0
    @. density_field.values = (pressure_field.values+p_ref) / (R * T_ref) # NOT SURE ABOUT THIS
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

function construct_acceleration_field(U_m, U_m_prev, props, dt) #, workgroup
    a = VectorField(mesh)
    g = props.DriftVelocity.gravity

    x0, y0, z0 = g[1], g[2], g[3]

    x = ScalarField(mesh)
    y = ScalarField(mesh)
    z = ScalarField(mesh)
    initialise!(x, x0)
    initialise!(y, y0)
    initialise!(z, z0)
    G = VectorField(x, y, z, mesh)

    dUdt = VectorField(mesh)
    @. dUdt[i].values = ( U_m[i] - U_m_prev[i] ) / dt

    @. a[i] = g[i] - (U[i] * ∇U.result[i]) - dUdt[i]
    
    # AK.foreachindex(dUdt, min_elems=workgroup, block_size=workgroup) do i 
    #     dUdt[i] = ( U_m[i] - U_m_prev[i] ) / dt
    # end
    return a
end

function construct_drag_field(rho_q, mu_q, props) # here q = vapour (from kassemi paper)
    if ( hasproperty(props, :drag) ) && ( props.drag isa Drag_SchillerNaumann )
        d_p = props.DriftVelocity.d_p
        # Re = compute
    end
end