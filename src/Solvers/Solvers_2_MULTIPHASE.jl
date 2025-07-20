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




# Temporary definitions, taken from openFoam damBreak just to quickly validate the results

#0 = air, 1 = water

rho_0 = 1
rho_1 = 1000

#nu is used in OF:

nu_0 = 1.48e-5
nu_1 = 1.0e-6

# mu_0 = 
# mu_1 = 

sigma = 0.07

# Temporary definitions



function setup_multiphase_solvers(
    solver_variant, model, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    ) 

    (; solvers, schemes, runtime, hardware, boundaries) = config

    @info "Extracting configuration and input fields..."

    (; U, p, Uf, pf) = model.momentum
    (; rho, alpha, alphaf) = model.fluid
    mesh = model.domain

    @info "Pre-allocating fields..."
    
    # ∇p = Grad{schemes.p.gradient}(p)
    # rhorDf = FaceScalarField(mesh)
    # initialise!(rhorDf, 1.0)
    # mueff = FaceScalarField(mesh)
    # mueffgradUt = VectorField(mesh)
    # divHv = ScalarField(mesh)
    # ddtrho = ScalarField(mesh)
    # psidpdt = ScalarField(mesh)
    # divmdotf = ScalarField(mesh)
    # psi = ScalarField(mesh)
    
    TF = _get_float(mesh)
    time = zero(TF) # assuming time=0



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


    x = ScalarField(mesh)
    y = ScalarField(mesh)
    z = ScalarField(mesh)

    initialise!(y, -9.81)

    # rho_water = ScalarField(mesh)
    # rho_air = ScalarField(mesh)

    # initialise!(rho_air, rho_0)
    # initialise!(rho_water, rho_1)
    interpolate!(alphaf, alpha, config)

    flux!(phif, Uf, config)
    # flux!(mdotf, Uf, config)
    # flux!(phiC, Uc, config)

    # flux!(phif, Uf, config)  
    flux!(phif, Uf, alphaf, config)              # (ASK HUMBERTO ABOUT THIS)
    # @. phif.values = phif.values * alphaF.values (ASK HUMBERTO ABOUT THIS)

    @. rho.values = (rho_1 * alpha.values) + (rho_0 * (1 - alpha.values)) # this is rho_m
    @. nueff.values = (nu_1 * alphaf.values) + (nu_0 * (1 - alphaf.values))
    @. y.values *= rho.values 

    interpolate!(rhof, rho, config)
    flux!(mdotf, Uf, rhof, config)

    # g = VectorField(mesh, (0.0, -9.81, 0.0))
    rhoG = VectorField(x, y, z, mesh)

    interpolate!(n_f, ∇alpha.result, config)
    correct_boundaries!(n_f, ∇alpha, boundaries.alpha, time, config)
    # @. n_f.values = n_f.values / norm(n_f.values) # normalised grad alpha
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



    
    # Trying to construct Navier Stokes:

    ### VOF TRANSPORT EQN.:
    
    

    # Assume that all terms but sources need to go on the left side, modify signs accordingly?

    #what is alpha_i?

    println(schemes.alpha.divergence)
    
    alpha_eqn = ( # Volume Fraction Transport Eqn
        Time{schemes.alpha.time}(alpha)
        + Divergence{schemes.alpha.divergence}(phif, alpha) #phif is U * n * alpha
        ==
        - Source(zero_field) #zero_field
    ) → ScalarEquation(alpha, boundaries.alpha)


    U_eqn = ( # Momentum Eqn
        Time{schemes.U.time}(rho, U) # rho * U
        + Divergence{schemes.U.divergence}(mdotf, U)  #mdotf is rho*U*n
        - Laplacian{schemes.U.laplacian}(nueff, U) # nu_mixture
        == 
        - Source(∇p.result) # -∇p
        # + Source(rhoG) # (0, -9.81 * rho, 0)
        # + Source(Mk) # leave it for now
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

    @. prgh.values = p.values - rho.values * g # * yCoords
    grad!(∇prgh, pf, prgh, boundaries.p, time, config)
    limit_gradient!(schemes.p.limiter, ∇prgh, prgh, config)


    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."
     
    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))

    # @info "Initialising energy model..."
    # energyModel = initialise(model.energy, model, mdotf, rho, p_eqn, config)

    @info "Initialising turbulence model..."
    turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        model, turbulenceModel, alpha_eqn, ∇p, U_eqn, p_eqn, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        time, rDf, nueff, phif, prgh, ∇prgh
        )

    return residuals    
end # end function







function MULTIPHASE(
    model, turbulenceModel, alpha_eqn, ∇p, U_eqn, p_eqn, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2,
    time, rDf, nueff, phif, prgh, ∇prgh
    )

    g = 9.81

    # T = 200 # dummy placeholder
    # p = 10e6 # dummy placeholder

    ### THERMOPHYSICAL SECTION
    # Update fluid properties (rho, mu, cp, k, sigma, ?) for all fractions (e.g. liquid and vapour)
    # Dummy values for now:

    # supercritical vs multiphase????
    # rho = 50
    # mu = 7e-6
    # cp = 15000
    # cv = 6750
    # k = 0.12208
    # sigma = 0.003

    # Get mixture values (blend) e.g. rho_m


    ### VOLUME FRACTION SECTION
    # Solve volume fraction transport eqn.
    # Any other VoF-related steps?????


    ### SUBGRID NUCLEATE BOILING SECTION
    # Nucleate boiling subgrid handling => get source terms


    ### TURBULENCE SECTION
    # Solve additional turbulence eqns.


    ### MOMENTUM-CONTINUITY SECTION
    # PISO loop where momentum eqn. is solved and corrected w.r.t. continuity

    
    # Solve energy eqn.


    # Next iteration!








# [DONE] Laplace validate + look up whats happening (analytical solution comparison)
# [DONE!] EoS multiparam
# Try to replicate cases for multiphase
# OpenFoam look into code





    ###############################################################


    (; U, p, Uf, pf) = model.momentum
    # (; nu) = model.fluid
    (; rho, rhof, alpha, alphaf) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries) = config
    (; iterations, write_interval, dt) = runtime
    (; backend) = hardware

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

        # println("alpha eqn:")
        # println(alpha_eqn)
        # println("alpha:")
        # println(alpha)
        # println("alpha boundaries:")
        # println(boundaries.alpha)
        # println("alpha solvers:")
        # println(solvers.alpha)
        


        interpolate!(alphaf, alpha, config)
        flux!(phif, Uf, alphaf, config)
        @. rho.values = (rho_1 * alpha.values) + (rho_0 * (1 - alpha.values)) # this is rho_m
        @. nueff.values = (nu_1 * alphaf.values) + (nu_0 * (1 - alphaf.values))

        interpolate!(rhof, rho, config)


        @. prgh.values = p.values - rho.values * g # * yCoords
        grad!(∇prgh, pf, prgh, boundaries.p, time, config)
        limit_gradient!(schemes.p.limiter, ∇prgh, prgh, config)


            
        # println("ralpha")
        # ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config; time=time)

            # rp = solve_equation!(p_eqn, p, boundaries.p, solvers.p, config; ref=pref, time=time)
        # break
        # println(ralpha)
        # break
        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config; time=time)
        # break
        # break
        # break
          
        # println(rx)
        # println(ry)
        # println(rz)


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

            # nonorthogonal correction (experimental)
            for i ∈ 1:ncorrectors
                discretise!(p_eqn, p, config)       
                apply_boundary_conditions!(p_eqn, boundaries.p, nothing, time, config)
                setReference!(p_eqn, pref, 1, config)
                nonorthogonal_face_correction(p_eqn, ∇p, rDf, config)
                update_preconditioner!(p_eqn.preconditioner, p.mesh, config)
                rp = solve_system!(p_eqn, solvers.p, p, nothing, config)

                if i == ncorrectors
                    explicit_relaxation!(p, prev, 1.0, config)
                else
                    explicit_relaxation!(p, prev, solvers.p.relax, config)
                end
                grad!(∇p, pf, p, boundaries.p, time, config) 
                limit_gradient!(schemes.p.limiter, ∇p, p, config)
            end

            # old approach - keep for now!
            # correct_velocity!(U, Hv, ∇p, rD, config)
            # interpolate!(Uf, U, config)
            # correct_boundaries!(Uf, U, boundaries.U, time, config)
            # flux!(mdotf, Uf, config) # old approach

            # new approach
            interpolate!(Uf, U, config) # velocity from momentum equation
            correct_boundaries!(Uf, U, boundaries.U, time, config)
            flux!(mdotf, Uf, config)
            correct_mass_flux(mdotf, p, rDf, config)
            correct_velocity!(U, Hv, ∇p, rD, config)

        end # corrector loop end
        
        # correct_mass_flux(mdotf, p, rDf, config) # new approach

    # turbulence!(turbulenceModel, model, S, prev, time, config) 

    println("ralpha")
    ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config; time=time)

    update_nueff!(nueff, nueff, model.turbulence, config)

    # if typeof(mesh) <: Mesh3
    #     residual!(R_uz, U_eqn, U.z, iteration, zdir, config)
    # end
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
    # println(ralpha)

    if iteration%write_interval + signbit(write_interval) == 0
        # save_output(model, outputWriter, time, config)
        save_output(model, outputWriter, iteration, time, config)
    end

    end # end for loop

    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p, alpha=R_alpha)
end