export multiphase!

export blend_properties!


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
    (; alpha, alphaf, rho, rhof) = model.fluid

    phases = model.fluid.phases
    props = model.fluid.physics_properties

    backend = hardware.backend
    workgroup = hardware.workgroup

    mesh = model.domain
    isInit = true
    

    @info "Pre-allocating fields..."
    
    TF = _get_float(mesh)
    time = zero(TF) # assuming time=0


    ∇p = Grad{schemes.p.gradient}(p)
    grad!(∇p, pf, p, boundaries.p, time, config)
    limit_gradient!(schemes.p.limiter, ∇p, p, config)

    ∇U = Grad{schemes.U.gradient}(U)
    grad!(∇U, Uf, U, boundaries.U, time, config)
    limit_gradient!(schemes.U.limiter, ∇U, U, config)

    phif = FaceScalarField(mesh)
    mdotf = FaceScalarField(mesh)
    nueff = FaceScalarField(mesh)

    interpolate!(alphaf, alpha, config)

    rho_l = phases[1].rho
    rhof_l = FaceScalarField(mesh)
    interpolate!(rhof_l, rho_l, config)


    flux!(phif, Uf, rhof_l, config)


    @info "Computing Fluid Properties..."

    phase_eos = [phases[1].eosModel, phases[2].eosModel]

    if typeof(model.energy) <: Nothing # Isothermal
        T = ScalarField(mesh)
        initialise!(T, 273.0) # THIS PROBABLY NEEDS TO BE DEFINED BY USER! Redesign Isothermal Energy ?
    end

    update_phase_thermodynamics!(phase_eos[1], phases[1], nueff, T, model, config) # For 3+ phases can just wrap in a for loop
    update_phase_thermodynamics!(phase_eos[2], phases[2], nueff, T, model, config)

    blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)

    for property in props
        update_extra_physics!(property, ∇U, model, config, isInit, 1.0, mesh) # currently computes v_pq and updates velocities in drift model
    end

    flux!(mdotf, Uf, rhof, config)

    psi = ScalarField(mesh)
    initialise!(psi, 1.0)

    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)

    divHv = ScalarField(mesh)
    div!(divHv, mdotf, config)


    @info "Defining models..."

    zero_field = ConstantScalar(0.0)


    alpha_rhs = - Source(zero_field)

    alpha_eqn = ( # Volume Fraction Transport Eqn
        Time{schemes.alpha.time}(rho_l, alpha)
        + Divergence{schemes.alpha.divergence}(phif, alpha)
        == alpha_rhs
    ) → ScalarEquation(alpha, boundaries.alpha)

    update_alpha_sources!(props, alpha_rhs, alpha_eqn, alpha, rho, mesh, config, phases, isInit)

    # sub Source <: Src

    # alpha_rhs = - Source(zero_field) - Source(1) + Source(2) - Source(4)
    

   momentum_rhs = - Source(∇p.result)

   U_eqn = ( # Momentum Eqn
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(mdotf, U)
        - Laplacian{schemes.U.laplacian}(nueff, U) #turbulence VS laminar ?
        == momentum_rhs
    ) → VectorEquation(U, boundaries.U)

    momentum_rhs = update_momentum_sources!(props, momentum_rhs, U_eqn, alpha, rho, mesh, config, phases, isInit)

    println("CHECKING ON SOURCES")
    println(momentum_rhs)



    p_eqn = (
       Time{schemes.p.time}(psi, p) # take psi (compressibility) from Helmholtz EoS if possible?
      - Laplacian{schemes.p.laplacian}(rDf, p) # something to do with momentum - figure this out!
      ==
      - Source(divHv)
    ) → ScalarEquation(p, boundaries.p)



    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."
     
    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))


    @info "Initialising turbulence model..."
    turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)

    if typeof(model.energy) <: MultiphaseEnergy
        @info "Initialising energy model..."

        m_qp = props.leeModel.m_qp
        m_pq = props.leeModel.m_pq
        latentHeat = props.leeModel.latentHeat

        nut = model.turbulence.nut

        # k_eff? dt ?
        energyModel = initialise(model.energy, model, nut, latentHeat, m_pq, m_qp, ∇p, config, mesh, 1.0) #ASSUME LEE MODEL IS TURNED ON
    end

    residuals  = solver_variant(
        model, turbulenceModel, alpha_eqn, ∇p, ∇U, U_eqn, p_eqn, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        time, rDf, nueff, phif, prgh, ∇prgh, mdotf, divHv
        )

    return residuals    
end # end function





function MULTIPHASE(
    model, turbulenceModel, alpha_eqn, ∇p, ∇U, U_eqn, p_eqn, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2,
    time, rDf, nueff, phif, prgh, ∇prgh, mdotf, divHv
    )

    # inner_loops = 2 DOUBLE CHECK IF THIS IS FINE
    
    (; U, p, Uf, pf) = model.momentum
    (; alpha, alphaf, rho, rhof) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries) = config
    (; iterations, write_interval, dt) = runtime
    (; backend) = hardware

    phases = model.fluid.phases
    props = model.fluid.physics_properties

    backend = hardware.backend
    workgroup = hardware.workgroup
    
    isInit = false
    phase_eos = [phases[1].eosModel, phases[2].eosModel]

    mesh = model.domain

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
    
    #REWORK THIS BIT OF CODE!!!
    if typeof(model.energy) <: Nothing # Isothermal
        T = ScalarField(mesh)
        initialise!(T, 273.0) # THIS PROBABLY NEEDS TO BE DEFINED BY USER! Redesign Isothermal Energy ?
    end

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting Multiphase Solver..."

    progress = Progress(iterations; dt=1.0, showspeed=true)

    @time for iteration ∈ 1:iterations
        time = iteration *dt

        # BEFORE PISO LOOPS

        update_phase_thermodynamics!(phase_eos[1], phases[1], nueff, T, model, config) # For 3+ phases can just wrap in a for loop
        update_phase_thermodynamics!(phase_eos[2], phases[2], nueff, T, model, config)

        blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)

        for property in props
            update_extra_physics!(property, ∇U, model, config, isInit, 1.0, mesh) # currently computes v_pq and updates velocities in drift model
        end

        
        update_momentum_sources!(props, 0.0, U_eqn, alpha, rho, mesh, config, phases, isInit)
        update_alpha_sources!(props, 0.0, alpha_eqn, alpha, rho, mesh, config, phases, isInit)
        
        # BEFORE PISO LOOPS

        @info "Starting PISO loops..."
        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config; time=time)

        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config)
        interpolate!(rDf, rD, config)
        remove_pressure_source!(U_eqn, ∇p, config)
        
        rp = 0.0
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

        # T_temp = ScalarField(mesh)
        # initialise!(T_temp, 70.0)

        # DUMMY FIELD

        #update_thermo_properties!(phase1 phase2 etc)
        
        # if typeof(model.fluid.phases[1].eos) <: HelmholtzEnergy #HelmholtzEnergyEOS
        #     helmholtz_energy!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, model, config)
        # else
        #     update_properties!(rho_fractions, nu_fractions, nuf_fractions, nueff, T_temp, phases, model, config)
        # end


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








# function update_momentum_sources!(props, rhs, alpha, rho, mesh, config, phases)
#     if ( hasproperty(props, :gravity) ) && ( props.gravity isa Gravity )
#         g = props.gravity.g
#         rho_p = phases[1].rho
#         rho_q = phases[2].rho

#         # rho = ScalarField(mesh)
#         # blend_properties!(rho, alpha, rho_p, rho_q)

#         x0, y0, z0 = g[1], g[2], g[3]
        
#         x = ScalarField(mesh)
#         y = ScalarField(mesh)
#         z = ScalarField(mesh)
#         initialise!(x, x0)
#         initialise!(y, y0)
#         initialise!(z, z0)

#         @. x.values = x.values * rho.values
#         @. y.values = y.values * rho.values
#         @. z.values = z.values * rho.values

#         rhoG = VectorField(x, y, z, mesh)

#         @info "Adding Gravity term in Momentum Equation..."
#         rhs -= Source(rhoG)
#     end

#     if ( hasproperty(props, :driftVelocity) ) && ( props.driftVelocity isa DriftVelocity )
#         schemes = config.schemes
#         backend = config.hardware.backend
#         workgroup = config.hardware.workgroup

#         (vdr_p, vdr_q) = (props.driftVelocity.v_dr_p, props.driftVelocity.v_dr_q)

#         rho_p = phases[1].rho
#         rho_q = phases[2].rho

#         alpha = alpha.values
#         rho_p = rho_p.values
#         rho_q = rho_q.values
#         # prob would make sense to put drift velocities inside DriftVelocity() object

#         slipVelocity = VectorField(mesh)
#         slipVelocityTensor = TensorField(mesh)

#         ndrange = length(slipVelocityTensor)
#         kernel! = _momentum_driftVelocity!(_setup(backend, workgroup, ndrange)...)
#         kernel!(slipVelocityTensor, vdr_p, vdr_q, alpha, rho_p, rho_q)

#         div_tensor!(slipVelocityTensor, slipVelocity, mesh, config)

#         @info "Adding Slip Velocity term in Momentum Equation..."

#         rhs -= Source(slipVelocity)
#     end

#     # return rhs
# end
@kernel inbounds=true function _momentum_driftVelocity!(slipVelocityTensor, vdr_p, vdr_q, alpha, rho_p, rho_q)
    i = @index(Global)

    cross_p = vdr_p[i] * vdr_p[i]'
    cross_q = vdr_q[i] * vdr_q[i]'
    slipVelocityTensor[i] = (alpha[i] * rho_p[i] * cross_p) + ((1.0-alpha[i]) * rho_q[i] * cross_q) #must be cross products
end




function update_extra_physics!(property, ∇U, model, config, isInit, dt, mesh)
    return nothing
end
function update_extra_physics!(property::DriftVelocityState, ∇U, model, config, isInit, dt, mesh)
    (; U) = model.momentum
    (; alpha) = model.fluid
    (; hardware) = config

    phases = model.fluid.phases
    props = model.fluid.physics_properties

    backend = hardware.backend
    workgroup = hardware.workgroup

    U_prev      = property.U_prev
    v_pq        = property.v_pq
    v_pq_prev   =  property.v_pq_prev

    if isInit
        U_prev = U # maybe perturb e.g. * 0.9
    else
        U_prev = U
    end
    
    a_field = construct_acceleration_field(U, ∇U, U_prev, dt, mesh, props, backend, workgroup)
    slip_velocity!(alpha, a_field, v_pq, v_pq_prev, mesh, model, phases, config, isInit)

end



struct InitialisationMode end
struct UpdateMode end

function compute_alpha_source(model::DriftVelocityState, alpha, rho, phases, config, mesh)

    schemes = config.schemes
    backend = config.hardware.backend
    workgroup = config.hardware.workgroup

    vdr_p = model.v_dr_p
    rho_p = phases[1].rho
    rho_q = phases[2].rho

    # rho = ScalarField(mesh)
    # blend_properties!(rho, alpha, rho_p, rho_q)

    alpha = alpha.values
    rho_p = rho_p.values

    div_term = ScalarField(mesh)

    inside_divergence = VectorField(mesh)

    ndrange = length(inside_divergence)
    kernel! = _alpha_driftVelocity!(_setup(backend, workgroup, ndrange)...)
    kernel!(inside_divergence, vdr_p, alpha, rho_p)
    
    div_interpolated = FaceVectorField(mesh)
    interpolate!(div_interpolated, inside_divergence, config)
    div!(div_term, div_interpolated, config)

    return [div_term, 1.0] # Source, sign
end

function compute_alpha_source(model::LeeModelState, alpha, rho, phases, config, mesh)
    return nothing
    # Special case - do nothing
end


function process_alpha_source!(::InitialisationMode, model::LeeModelState, alpha_eqn, rhs, source, pointerID) ## SPECIAL CASE FOR LEE MODEL
    
    m_qp = model.m_qp
    m_pq = model.m_pq
    
    @info "Adding Lee Model in Alpha Equation..."
    
    rhs += Source(m_qp)
    rhs -= Source(m_pq)
end
function process_alpha_source!(::InitialisationMode, model, alpha_eqn, rhs, source, pointerID) ## DURING INITIALISATION
    @info "Adding source for $(nameof(typeof(model))) in Alpha Equation..."

    source_field    = source[1]
    source_sign     = source[2]

    if source_sign > 0
        rhs += Source(source_field)
    else
        rhs -= Source(source_field)
    end

    model.alpha_source_pointer[] = pointerID
end

function process_alpha_source!(::UpdateMode, model, alpha_eqn, rhs, source_object, pointerID) ## DURING SOLVING
    source_field = source_object[1]

    pointer = model.alpha_source_pointer[]
    temp_source = get_source(alpha_eqn, pointer)
    temp_source = source_field # Check if this is not a Source() object
end

function update_alpha_sources!(props, rhs, eqn, alpha, rho, mesh, config, phases, is_init::Bool)
    
    if is_init
        mode = InitialisationMode()
    else
        mode = UpdateMode()
    end

    filtered_models = (m for m in props if hasproperty(m, :alpha_source_pointer)) # sort properties that are only relevant to momentum eqn.

    for (i, model) in enumerate(filtered_models)
        source_object = compute_alpha_source(model, alpha, rho, phases, config, mesh)
        process_alpha_source!(mode, model, eqn, rhs, source_object, i+1) #i+1 because first source is always defined
    end
end




function compute_momentum_source(model::DriftVelocityState, alpha, rho, phases, config, mesh)
    schemes = config.schemes
    backend = config.hardware.backend
    workgroup = config.hardware.workgroup

    (vdr_p, vdr_q) = (model.v_dr_p, model.v_dr_q)

    rho_p = phases[1].rho
    rho_q = phases[2].rho

    alpha = alpha.values
    rho_p = rho_p.values
    rho_q = rho_q.values
    # prob would make sense to put drift velocities inside DriftVelocity() object

    slipVelocity = VectorField(mesh)
    slipVelocityTensor = TensorField(mesh)

    ndrange = length(slipVelocityTensor)
    kernel! = _momentum_driftVelocity!(_setup(backend, workgroup, ndrange)...)
    kernel!(slipVelocityTensor, vdr_p, vdr_q, alpha, rho_p, rho_q)

    div_tensor!(slipVelocityTensor, slipVelocity, mesh, config)

    return [slipVelocity, -1.0] # Source, sign
end

function compute_momentum_source(model::GravityState, alpha, rho, phases, config, mesh)
    g = model.g
    rho_p = phases[1].rho
    rho_q = phases[2].rho

    # rho = ScalarField(mesh)
    # blend_properties!(rho, alpha, rho_p, rho_q)

    x0, y0, z0 = g[1], g[2], g[3]
    
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

    return [rhoG, 1.0] # Source, sign
end

function process_momentum_source(::InitialisationMode, model, config, U_eqn, rhs, source, pointerID) ## DURING INITIALISATION
    @info "Adding source for $(nameof(typeof(model))) in Momentum Equation..."

    source_field    = source[1]
    source_sign     = source[2]

    println(rhs)    

    if source_sign > 0
        rhs += Source(source_field)
    else
        rhs -= Source(source_field)
    end
    # push_to_RHS!(source_field, source_sign, rhs, config)
    
    println("AFTER ADDING")
    println(rhs)
    model.momentum_source_pointer[] = pointerID

    return rhs
end
function process_momentum_source(::UpdateMode, model, config, U_eqn, rhs, source_object, pointerID) ## DURING SOLVING
    source_field = source_object[1]

    pointer = model.momentum_source_pointer[]
    println("Pointer from the model: ", pointer)
    println("Pointer ID: ", pointerID)
    temp_source = get_source(U_eqn, pointer)
    temp_source = source_field # Check if this is not a Source() object

    return nothing
end

function update_momentum_sources!(props, rhs, eqn, alpha, rho, mesh, config, phases, is_init::Bool)
    
    if is_init
        mode = InitialisationMode()
    else
        mode = UpdateMode()
    end

    filtered_models = (m for m in props if hasproperty(m, :momentum_source_pointer)) # sort properties that are only relevant to momentum eqn.

    for (i, model) in enumerate(filtered_models)
        println("INSIDE THE LOOP, MODE : $mode")
        source_object = compute_momentum_source(model, alpha, rho, phases, config, mesh)
        rhs = process_momentum_source(mode, model, config, eqn, rhs, source_object, i+1) #i+1 because first source is always defined
    end

    return rhs
end



function push_to_RHS!(source_field::ScalarField, source_sign, RHS, config)
    @. source_field.values = source_field.values * source_sign

    push!(RHS, Source(source_field))
end
function push_to_RHS!(source_field::VectorField, source_sign, RHS, config)
    backend = config.hardware.backend
    workgroup = config.hardware.workgroup

    ndrange = length(source_field)
    kernel! = _push_to_RHS!(_setup(backend, workgroup, ndrange)...)
    kernel!(source_field, source_sign)
    
    push!(RHS, Source(source_field))
end
@kernel inbounds=true function _push_to_RHS!(source_field, source_sign)
    i = @index(Global)

    source_field[i] = source_field[i] * source_sign
end


# function update_alpha_sources!(props, rhs, alpha, rho, mesh, config, phases)
#     if ( hasproperty(props, :driftVelocity) ) && ( props.driftVelocity isa DriftVelocity )
#         schemes = config.schemes
#         backend = config.hardware.backend
#         workgroup = config.hardware.workgroup

#         vdr_p = props.driftVelocity.v_dr_p
#         rho_p = phases[1].rho
#         rho_q = phases[2].rho

#         # rho = ScalarField(mesh)
#         # blend_properties!(rho, alpha, rho_p, rho_q)

#         alpha = alpha.values
#         rho_p = rho_p.values

#         div_term = ScalarField(mesh)

#         inside_divergence = VectorField(mesh)

#         ndrange = length(inside_divergence)
#         kernel! = _alpha_driftVelocity!(_setup(backend, workgroup, ndrange)...)
#         kernel!(inside_divergence, vdr_p, alpha, rho_p)
        
#         div_interpolated = FaceVectorField(mesh)
#         interpolate!(div_interpolated, inside_divergence, config)
#         div!(div_term, div_interpolated, config)

#         # EITHER DO THIS
#         @info "Adding Slip Velocity term in Alpha Equation..."
        
#         rhs -= Source(div_term)

#         #OR JUST GET THE FIELD AND UPDATE IT
#     end
#     if ( hasproperty(props, :leeModel) ) && ( props.leeModel isa LeeModel )
#         m_qp = props.leeModel.m_qp
#         m_pq = props.leeModel.m_pq
        
#         @info "Adding Lee Model in Alpha Equation..."
        
#         rhs += Source(m_qp)
#         rhs -= Source(m_pq)
#     end
    
#     # return rhs
# end
@kernel inbounds=true function _alpha_driftVelocity!(inside_divergence, vdr_p, alpha, rho_p)
    i = @index(Global)

    inside_divergence[i] = vdr_p[i] * alpha[i] * rho_p[i]
end











function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0.values * alpha_field.values) + (property_1.values * (1.0 - alpha_field.values))
    nothing
end




function div_tensor!(tensor_field, output_vector, mesh, config)
    resX = ScalarField(mesh)
    resY = ScalarField(mesh)
    resZ = ScalarField(mesh)

    r1 = VectorField(mesh)
    r2 = VectorField(mesh)
    r3 = VectorField(mesh)

    r1f = FaceVectorField(mesh)
    r2f = FaceVectorField(mesh)
    r3f = FaceVectorField(mesh)

    r1.x.values .= tensor_field.xx.values
    r1.y.values .= tensor_field.xy.values
    r1.z.values .= tensor_field.xz.values

    r2.x.values .= tensor_field.yx.values
    r2.y.values .= tensor_field.yy.values
    r2.z.values .= tensor_field.yz.values

    r3.x.values .= tensor_field.zx.values
    r3.y.values .= tensor_field.zy.values
    r3.z.values .= tensor_field.zz.values
    
    interpolate!(r1f, r1, config)
    interpolate!(r2f, r2, config)
    interpolate!(r3f, r3, config)

    div!(resX, r1f, config)
    div!(resY, r2f, config)
    div!(resZ, r3f, config)

    # println(output_vector)
    # println(output_vector)
    output_vector.x.values .= resX.values
    output_vector.y.values .= resY.values
    output_vector.z.values .= resZ.values

    # nothing
end


function slip_velocity!(alpha, a_field, v_pq, v_pq_prev, mesh, model, phases, config, isInitialisation)
    (; U) = model.momentum

    rho = model.fluid.rho
    props = model.fluid.physics_properties

    d_p = ScalarField(mesh)
    initialise!(d_p, props.driftVelocity.d_p)

    # max_U = max_field_value(U)
    max_U = max_BC_value(config.boundaries.U)


    find_Vpq!(alpha, rho, d_p, a_field, v_pq_prev, isInitialisation, max_U, v_pq, mesh, phases, config)
    
    update_velocities!(U, props, alpha, mesh, phases, config, rho) # liquid <> p ; vapour <> q
end


function construct_acceleration_field(U_m, ∇U_m, U_m_prev, dt, mesh, props, backend, workgroup) #, workgroup
    grad_U = ∇U_m.result
    a = VectorField(mesh)

    g = props.driftVelocity.gravity.g

    x0, y0, z0 = g[1], g[2], g[3]

    x = ScalarField(mesh)
    y = ScalarField(mesh)
    z = ScalarField(mesh)
    initialise!(x, x0)
    initialise!(y, y0)
    initialise!(z, z0)
    G = VectorField(x, y, z, mesh)

    dUdt = VectorField(mesh)


    ndrange = length(U_m)
    kernel! = _construct_acceleration_field!(_setup(backend, workgroup, ndrange)...)
    kernel!(U_m, dUdt, U_m_prev, dt, a, G, grad_U)
    
    return a
end
@kernel inbounds=true function _construct_acceleration_field!(U_m, dUdt, U_m_prev, dt, a, G, grad_U)
    i = @index(Global)

    dUdt[i] = (U_m[i] - U_m_prev[i]) / dt
    a[i] = G[i] - (grad_U[i] * U_m[i]) - dUdt[i]
end








function get_E!(E, alpha, rho, rho_p, rho_q, d_p, mu_p, a_field, backend, workgroup, mesh)
    # rho = ScalarField(mesh)
    # blend_properties!(rho, alpha, rho_p, rho_q)

    ndrange = length(E)
    kernel! = _get_E!(_setup(backend, workgroup, ndrange)...)
    kernel!(E, rho_p, d_p, mu_p, rho, a_field)


    return nothing
end
@kernel inbounds=true function _get_E!(E, rho_p, d_p, mu_p, rho, a_field)
    i = @index(Global)
 
    E[i] = (rho_p[i]*((d_p[i])^2))/(18*mu_p[i]) * ((rho_p[i]-rho[i])/rho_p[i]) * a_field[i]
end


function compute_RE!(RE, rho_q, v_pq, d_p, mu_q, backend, workgroup)
    ndrange = length(RE)
    kernel! = _compute_RE!(_setup(backend, workgroup, ndrange)...)
    kernel!(RE, rho_q, v_pq, d_p, mu_q)
    

    return nothing
end
@kernel inbounds=true function _compute_RE!(RE, rho_q, v_pq, d_p, mu_q)
    i = @index(Global)
 
    RE[i] = (rho_q[i] * d_p[i] * norm(v_pq[i]))/mu_q[i]
end


function find_Vpq!(alpha, rho, d_p, a_field, v_pq_prev, isInitialisation, v_high, v_pq, mesh, phases, config)
    v_low = SVector(1.0e-9, 1.0e-9, 1.0e-9)
    backend = config.hardware.backend
    workgroup = config.hardware.workgroup

    rho_p = phases[1].rho
    rho_q = phases[2].rho

    mu_p = phases[1].mu
    mu_q = phases[2].mu


    E = VectorField(mesh)
    expr = ScalarField(mesh)
    v_pq_trial = VectorField(mesh)
    RE_trial = ScalarField(mesh)

    get_E!(E, alpha, rho, rho_p, rho_q, d_p, mu_p, a_field, backend, workgroup, mesh)

    @. expr.values = (mu_q.values*rho_q.values)/(0.0183 * rho_q.values * d_p.values)
    
    for i in eachindex(v_pq_trial)
        v_pq_trial[i] = sqrt.(E[i]*expr.values[i])
    end

    compute_RE!(RE_trial, rho_q, v_pq_trial, d_p, mu_q, backend, workgroup)

    K = ScalarField(mesh)
    @. K.values = 0.15 * ((rho_q.values * d_p.values)/mu_q.values)^0.687

    
    ndrange = length(RE_trial)
    kernel! = _find_Vpq!(_setup(backend, workgroup, ndrange)...)
    kernel!(RE_trial, v_pq, v_pq_trial, K, E, v_low, v_high, v_pq_prev, isInitialisation)
end
@kernel inbounds=true function _find_Vpq!(RE_trial, v_pq, v_pq_trial, K, E, v_low, v_high, v_pq_prev, isInitialisation)
    i = @index(Global)
 
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
            # println(v_pq)

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

function update_velocities!(U, props, alpha, mesh, phases, config, rho)
    v_pq = props.driftVelocity.v_pq
    
    v_dr_p  = props.driftVelocity.v_dr_p
    v_dr_q  = props.driftVelocity.v_dr_q
    v_p     = props.driftVelocity.v_p
    v_q     = props.driftVelocity.v_q

    backend = config.hardware.backend
    workgroup = config.hardware.workgroup
    rho_p = phases[1].rho
    rho_q = phases[2].rho

    # rho = ScalarField(mesh)
    # blend_properties!(rho, alpha, rho_p, rho_q)

    alpha = alpha.values
    rho_q = rho_q.values
    rho_p = rho_p.values
    rho = rho.values

    C_p = ScalarField(mesh)
    C_q = ScalarField(mesh)

    @. C_p.values = (alpha * rho_p) / rho          #liquid
    @. C_q.values = ((1.0-alpha) * rho_q) / rho    #vapour

    ndrange = length(U)
    kernel! = _update_velocities!(_setup(backend, workgroup, ndrange)...)
    kernel!(U, v_dr_p, v_dr_q, v_p, v_q, C_p, C_q, v_pq)
end
@kernel inbounds=true function _update_velocities!(U, v_dr_p, v_dr_q, v_p, v_q, C_p, C_q, v_pq)
    i = @index(Global)

    v_dr_p[i] = (1.0 + C_p[i]) * v_pq[i]
    v_dr_q[i] = -(1.0 + C_q[i]) * v_pq[i]

    v_p[i] = v_dr_p[i] - U[i]
    v_q[i] = v_dr_q[i] - U[i]
end





function peng_robinson!(T, p, T_c, p_c, ω, M)
    # using Polynomials

    # T_c = 33.145 # H2
    # p_c = 1.2964e6 # H2
    R_univ = 8.314
    M = M * 1.0e-3
    # ω = -0.216 # H2

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


    # Rho 
    # Cp 
    # Beta
    # m_qp, m_pq, latentHeat (Enthalpy?)

    return rho_liq, rho_vap # liq = vap for single phase so its handled by alpha I suppose....
end



function update_phase_thermodynamics!(EoS::Union{ConstEos, PerfectGas}, phase, nueff, T, model, config)
    phase.eosModel(phase, model)
    phase.viscosityModel(phase, T)
end



#SPECIALISE TO AVOID IF

    # model::Physics{T,F,SO,M,Tu,E,D,BI}, config; 
    # output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    # ) where{T,F<:Multiphase,SO,M,Tu,E,D,BI} =  #T<:Transient ???


function update_phase_thermodynamics!(EoS::HelmholtzEnergy, phase, nueff, T, model, config)
    if phase == model.fluid.phases[1] # IF IT IS SECOND HELMHOLTZ CALL WE JUST IGNORE IT
        (; p) = model.momentum
        alpha = model.fluid.alpha

        HelmholtzModel = phase.eosModel
        rho_field = phase.rho

        backend = config.hardware.backend
        workgroup = config.hardware.workgroup

        lee_model_state = if hasproperty(model.fluid.physics_properties, :leeModel)
            fluid.physics_properties.leeModel
        else
            nothing
        end
        
        ndrange = length(rho_field)
        kernel! = _helmholtz_kernel!(_setup(backend, workgroup, ndrange)...)
        kernel!(model.fluid, HelmholtzModel, T, p, alpha, lee_model_state)
    end
end
@kernel inbounds=true function _helmholtz_kernel!(fluid, HelmholtzModel, T, p, alpha, lee_model_state)
    i = @index(Global)

    _, rho_temp, _, cp_temp, _, _, _, beta_temp, mu_temp, k_temp, _, latentHeat_temp, _, m_lv_temp, m_vl_temp = HelmholtzModel(T[i], p[i], alpha[i])
            
    fluid.phases[1].rho[i] = rho_temp[1]
    fluid.phases[2].rho[i] = rho_temp[2]
            
    fluid.phases[1].mu[i] = mu_temp[1]
    fluid.phases[2].mu[i] = mu_temp[2]
            
    fluid.phases[1].cp[i] = cp_temp[1]
    fluid.phases[2].cp[i] = cp_temp[2]

    fluid.phases[1].beta[i] = beta_temp[1]
    fluid.phases[2].beta[i] = beta_temp[2]

    fluid.phases[1].k[i] = k_temp[1]
    fluid.phases[2].k[i] = k_temp[2]

    _update_lee_model_fields!(lee_model_state, i, latentHeat_temp, m_lv_temp, m_vl_temp) #only works if lee model is there
end

# @inline should optimise it??
@inline function _update_lee_model_fields!(::Nothing, i, latentHeat, m_lv, m_vl) 
    return nothing
end

@inline function _update_lee_model_fields!(lee_model_state::LeeModelState, i, latentHeat, m_lv, m_vl)
    lee_model_state.latentHeat[i] = latentHeat
    lee_model_state.m_qp[i] = m_lv
    lee_model_state.m_pq[i] = m_vl
    
    return nothing
end




function max_BC_value(BC) # only for vector fields r
    # example: BC = config.boundaries.U
    output = SVector(0.0, 0.0, 0.0)

    for i in eachindex(BC)
        current_mag = norm(BC[i].value)
        output_mag = norm(output)

        if current_mag > output_mag
            output = BC[i].value
        end
    end

    return output
end


function max_field_value(field::ScalarField)
    return maximum(field.values)
end

function max_field_value(field::VectorField)
    max_x = maximum(field.x.values)
    max_y = maximum(field.y.values)
    max_z = maximum(field.z.values)

    return SVector(max_x, max_y, max_z)
end