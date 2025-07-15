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
    (; rho) = model.fluid
    mesh = model.domain

    @info "Pre-allocating fields..."
    
    ∇p = Grad{schemes.p.gradient}(p)
    mdotf = FaceScalarField(mesh)
    rhorDf = FaceScalarField(mesh)
    initialise!(rhorDf, 1.0)
    mueff = FaceScalarField(mesh)
    mueffgradUt = VectorField(mesh)
    divHv = ScalarField(mesh)
    ddtrho = ScalarField(mesh)
    psidpdt = ScalarField(mesh)
    divmdotf = ScalarField(mesh)
    psi = ScalarField(mesh)


    @info "Defining models..."



    rho_m = 0.0

    for k in n
        rho_m += rho_k 
    end



    # water portion = alpha_water
    # air portion = alpha_air




    
    # Trying to construct Navier Stokes:

    ### VOF TRANSPORT EQN.:

    phif = mdotf  # ρ * U * n  (already used for momentum), where U is averaged by mass
    
    

    # Assume that all terms but sources need to go on the left side, modify signs accordingly?

    #what is alpha_i?

    alpha_eqn = ( # Volume Fraction Transport Eqn
        Time{schemes.alpha.time}(1.0, alphai)
        + Divergence{schemes.alpha.div}(phif, alphai)          
        ==
        - Source(zero_field) 
    ) → ScalarEquation(alpha, boundaries.alpha)


    U_eqn = ( # Momentum Eqn
        Time{schemes.U.time}(alphaRho, U) #alpha * rho
        + Divergence{schemes.U.divergence}(mdotf, U)  #mdotf = 
        - Laplacian{schemes.U.laplacian}(mueff, U) #mueff = 
        == 
        - Source(∇p.result) #
        + Source(rhoAlphaG) #
        + Source(Mk) #
    ) → VectorEquation(U, boundaries.U)


    # P eqn too?






    # rho eqn doesn't work at the moment.
    # rho_eqn = (
    #     Time{schemes.rho.time}(rho) 
    #     == 
    #     -Source(divmdotf)
    # ) → ScalarEquation(mesh)

    # U_eqn = (
    #     Time{schemes.U.time}(rho, U)
    #     + Divergence{schemes.U.divergence}(mdotf, U) 
    #     - Laplacian{schemes.U.laplacian}(mueff, U) 
    #     == 
    #     - Source(∇p.result)
    #     + Source(mueffgradUt)
    # ) → VectorEquation(U, boundaries.U)

    # if typeof(model.fluid) <: WeaklyCompressible
        
    #     p_eqn = (
    #         Time{schemes.p.time}(psi, p)  # correction(fvm::ddt(p)) means d(p)/d(t) - d(pold)/d(t)
    #         - Laplacian{schemes.p.laplacian}(rhorDf, p)
    #         ==
    #         - Source(divHv)
    #     ) → ScalarEquation(p, boundaries.p)

    # elseif typeof(model.fluid) <: Compressible

    #     pconv = FaceScalarField(mesh)

    #     p_eqn = (
    #         Time{schemes.p.time}(psi, p)
    #         - Laplacian{schemes.p.laplacian}(rhorDf, p) 
    #         + Divergence{schemes.p.divergence}(pconv, p)
    #         ==
    #         -Source(divHv)
    #         -Source(ddtrho) # capture correction part of dPdT and explicit drhodt
    #     ) → ScalarEquation(p, boundaries.p)
    # end

    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p.preconditioner, p_eqn)

    @info "Pre-allocating solvers..."
     
    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p.solver, _b(p_eqn))

    @info "Initialising energy model..."
    energyModel = initialise(model.energy, model, mdotf, rho, p_eqn, config)

    @info "Initialising turbulence model..."
    turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        model, turbulenceModel, energyModel, ∇p, U_eqn, p_eqn, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops)

    return residuals    
end # end function







function MULTIPHASE(
    model, turbulenceModel, ∇p, U_eqn, p_eqn, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    )

    T = 200 # dummy placeholder
    p = 10e6 # dummy placeholder

    ### THERMOPHYSICAL SECTION
    # Update fluid properties (rho, mu, cp, k, sigma, ?) for all fractions (e.g. liquid and vapour)
    # Dummy values for now:

    # supercritical vs multiphase????
    rho = 50
    mu = 7e-6
    cp = 15000
    cv = 6750
    k = 0.12208
    sigma = 0.003

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

    





    # Extract model variables and configuration
    # (; U, p, Uf, pf) = model.momentum
    # (; nu) = model.fluid
    # mesh = model.domain
    # (; solvers, schemes, runtime, hardware, boundaries) = config
    # (; iterations, write_interval, dt) = runtime
    # (; backend) = hardware
    
    # mdotf = get_flux(U_eqn, 2)
    # nueff = get_flux(U_eqn, 3)
    # rDf = get_flux(p_eqn, 1)
    # divHv = get_source(p_eqn, 1)

    # outputWriter = initialise_writer(output, model.domain)
    
    # @info "Allocating working memory..."

    # # Define aux fields 
    # gradU = Grad{schemes.U.gradient}(U)
    # gradUT = T(gradU)
    # Uf = FaceVectorField(mesh)
    # S = StrainRate(gradU, gradUT, U, Uf)

    # n_cells = length(mesh.cells)
    # pf = FaceScalarField(mesh)
    # Hv = VectorField(mesh)
    # rD = ScalarField(mesh)

    # # Pre-allocate auxiliary variables
    # TF = _get_float(mesh)
    # TI = _get_int(mesh)
    # # prev = zeros(TF, n_cells)
    # # prev = _convert_array!(prev, backend) 
    # prev = KernelAbstractions.zeros(backend, TF, n_cells)

    # # Pre-allocate vectors to hold residuals 
    # R_ux = ones(TF, iterations)
    # R_uy = ones(TF, iterations)
    # R_uz = ones(TF, iterations)
    # R_p = ones(TF, iterations)
    # cellsCourant =adapt(backend, zeros(TF, length(mesh.cells)))
    
    # # Initial calculations
    # time = zero(TF) # assuming time=0
    # interpolate!(Uf, U, config)   
    # correct_boundaries!(Uf, U, boundaries.U, time, config)
    # flux!(mdotf, Uf, config)
    # grad!(∇p, pf, p, boundaries.p, time, config)
    # limit_gradient!(schemes.p.limiter, ∇p, p, config)

    # update_nueff!(nueff, nu, model.turbulence, config)

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

    turbulence!(turbulenceModel, model, S, prev, time, config) 
    update_nueff!(nueff, nu, model.turbulence, config)

    # if typeof(mesh) <: Mesh3
    #     residual!(R_uz, U_eqn, U.z, iteration, zdir, config)
    # end
    maxCourant = max_courant_number!(cellsCourant, model, config)

    R_ux[iteration] = rx
    R_uy[iteration] = ry
    R_uz[iteration] = rz
    R_p[iteration] = rp

    ProgressMeter.next!(
        progress, showvalues = [
            (:time, iteration*runtime.dt),
            (:Courant, maxCourant),
            (:Ux, R_ux[iteration]),
            (:Uy, R_uy[iteration]),
            (:Uz, R_uz[iteration]),
            (:p, R_p[iteration]),
            turbulenceModel.state.residuals...
            ]
        )

    if iteration%write_interval + signbit(write_interval) == 0
        # save_output(model, outputWriter, time, config)
        save_output(model, outputWriter, iteration, config)
    end

    end # end for loop

    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p)
end