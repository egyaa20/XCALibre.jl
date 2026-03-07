
export multiphase!

"""
    multiphase!(model, config; 
        output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0)

Multiphase solver for immiscible fluids. Solves coupled momentum, phase fraction (transport), and dynamic pressure equations.
Uses a "p_rgh" pressure formulation to handle gravity and hydrostatic pressure stability.

# Input arguments

- `model` reference to a `Physics` model defined by the user.
- `config` Configuration structure defined by the user with solvers, schemes, runtime and hardware structures configuration details.
- `output` select the format used for simulation results from `VTK()` or `OpenFOAM` (default = `VTK()`)
- `pref` Reference pressure value for cases that do not have a pressure defining BC. Incompressible solvers only (default = `nothing`)
- `ncorrectors` number of non-orthogonality correction loops (default = `0`)
- `inner_loops` number to inner loops used in transient solver based on PISO algorithm (default = `0`)

# Output

- `Ux` Vector of x-velocity residuals for each iteration.
- `Uy` Vector of y-velocity residuals for each iteration.
- `Uz` Vector of y-velocity residuals for each iteration.
- `p` Vector of pressure residuals for each iteration.
- `alpha` Vector of phase fraction residuals for each iteration.
"""
function multiphase!(
    model, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    )

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
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    ) 

    (; solvers, schemes, runtime, hardware, boundaries) = config

    @info "Extracting configuration and input fields..."

    (; U, p, Uf, pf) = model.momentum
    (; alpha, alphaf, rho, rhof, nu, nuf, p_rgh, p_rghf) = model.fluid

    phases = model.fluid.phases
    props = model.fluid.physics_properties

    backend = hardware.backend
    workgroup = hardware.workgroup
    mesh = model.domain

    @info "Pre-allocating fields..."

    TF = _get_float(mesh)
    time = zero(TF) # assuming time=0
    
    ∇p = Grad{schemes.p_rgh.gradient}(p)
    
    ∇p_rgh = Grad{schemes.p_rgh.gradient}(p_rgh)
    grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)
    
    ∇alpha = Grad{schemes.alpha.gradient}(alpha)
    grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
    limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

    mdotf = FaceScalarField(mesh)
    rhoPhi = FaceScalarField(mesh)
    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)
    nueff = FaceScalarField(mesh)
    divHv = ScalarField(mesh)
    
    phi_g = VectorField(mesh)
    phi_gf = FaceScalarField(mesh)
    divCompressionFlux = ScalarField(mesh)
    compressionFlux = FaceScalarField(mesh)

    interpolate_upwind!(alphaf, alpha, mdotf, config)

    # Needs to be defined before energyModel
    p_eqn = (
        - Laplacian{schemes.p.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)


    @info "Computing Fluid Properties..."


    phase_eos = [phases[1].density, phases[2].density]
    T_field = model.energy.T

    update_phase_thermodynamics!(phase_eos[1], Val(1), nueff, T_field, model, config)
    update_phase_thermodynamics!(phase_eos[2], Val(2), nueff, T_field, model, config)

    blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)
    blend_properties!(nu, alpha, phases[1].nu, phases[2].nu)

    blend_properties_at_faces!(rhof, alphaf, phases[1].density.rho, phases[2].density.rho)
    blend_properties_at_faces!(nuf, alphaf, phases[2].mu.mu / phases[1].density.rho, phases[2].mu.mu / phases[2].density.rho)

    gh = model.fluid.physics_properties.gravity.gh
    ghf = model.fluid.physics_properties.gravity.ghf
    g = model.fluid.physics_properties.gravity.g

    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)

    ∇rho = Grad{schemes.p_rgh.gradient}(rho)
    grad!(∇rho, rhof, rho, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

    ∇alpha = Grad{schemes.alpha.gradient}(alpha)
    grad!(∇alpha, alphaf, alpha, time, config)
    limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(rhoPhi, U)
        - Laplacian{schemes.U.laplacian}(nueff, U) 
        ==
        - Source(∇p_rgh.result)
    ) → VectorEquation(U, boundaries.U)

    alpha_eqn = (
        Time{schemes.alpha.time}(alpha)
        + Divergence{schemes.alpha.divergence}(mdotf, alpha)
        == 
        - Source(ConstantScalar(0.0))
    ) → ScalarEquation(alpha, boundaries.alpha)

    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p_rgh.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."
     
    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p_rgh.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))

    @info "Initialising turbulence model..."
    turbulenceModel, config = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        model, turbulenceModel, ∇p, ∇p_rgh, ∇rho, ∇alpha, U_eqn, p_eqn, alpha_eqn, mdotf, rhoPhi, gh, ghf, phi_g, phi_gf, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops)

    return residuals
end # end function



function MULTIPHASE(
    model, turbulenceModel, ∇p, ∇p_rgh, ∇rho, ∇alpha, U_eqn, p_eqn, alpha_eqn, mdotf, rhoPhi, gh, ghf, phi_g, phi_gf, config;

    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    )
    
    (; U, p, Uf, pf) = model.momentum
    (; nu, nuf, rho, rhof, alpha, alphaf, p_rgh, p_rghf) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries, postprocess) = config
    (; iterations, write_interval, dt) = runtime
    (; backend) = hardware
    
    if typeof(runtime.adaptive) <: Nothing
        maxAlphaCo = 0.75
    else
        (; maxCo, maxAlphaCo, maxGrow, minShrink) = runtime.adaptive
    end

    dt_cpu = zeros(_get_float(mesh), 1)
    copyto!(dt_cpu, config.runtime.dt)

    phases = model.fluid.phases
    
    rho1 = phases[1].rho
    rho2 = phases[2].rho

    rho1f = FaceScalarField(mesh)
    rho2f = FaceScalarField(mesh)

    postprocess = convert_time_to_iterations(postprocess,model,dt_cpu[1],iterations)

    # nueff = get_flux(U_eqn, 3)
    nueff = FaceScalarField(mesh) #Apparently doing this makes the difference in hydrostatic column test and lets U converge
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)

    @info "Allocating working memory..."

    # Define aux fields 
    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    Uf = FaceVectorField(mesh)
    S = StrainRate(gradU, gradUT, U, Uf)
    ∇p_rghf_deconstructed = FaceScalarField(mesh)
    ∇p_rghf_reconstructed = VectorField(mesh)
    rho_prev = ScalarField(mesh)
    alpha_prev = ScalarField(mesh)

    alphaf_upwind = FaceScalarField(mesh)
    alphaf_HO = FaceScalarField(mesh)

    ∇alphaf_upwind = FaceVectorField(mesh)
    ∇alphaf_HO = FaceVectorField(mesh)
    
    F_final = FaceScalarField(mesh)

    sigma = 71.1e-3
    # sigma = 0.0
    
    nhatf_prep = FaceVectorField(mesh)
    kappa = ScalarField(mesh)
    kappaf = FaceScalarField(mesh)

    n_cells = length(mesh.cells)
    Hv = VectorField(mesh)
    rD = ScalarField(mesh)

    # Pre-allocate auxiliary variables
    TF = _get_float(mesh)
    TI = _get_int(mesh)
    prev = KernelAbstractions.zeros(backend, TF, n_cells)

    # Pre-allocate vectors to hold residuals 
    R_ux = ones(TF, iterations)
    R_uy = ones(TF, iterations)
    R_uz = ones(TF, iterations)
    R_p = ones(TF, iterations)
    R_alpha = ones(TF, iterations)
    cellsCourant = KernelAbstractions.zeros(backend, TF, n_cells)
    cellsAlphaCourant = KernelAbstractions.zeros(backend, TF, n_cells)

    
    # Initial calculations
    time = zero(TF) # assuming time=0
    interpolate!(Uf, U, config)   
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    
    flux!(mdotf, Uf, config)
    @. rhoPhi.values = mdotf.values * rhof.values

    phase_eos = [phases[1].density, phases[2].density]
    T_field = model.energy.T

    update_nueff!(nueff, nuf, model.turbulence, config)

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting multiphase loops..."

    progress = Progress(iterations; dt=1.0, showspeed=true)

    sub_cycle = false


    @time for iteration ∈ 1:iterations
        copyto!(dt_cpu, config.runtime.dt)
        time += dt_cpu[1]
        @. alpha_prev.values = alpha.values
        grad!(∇alpha, alphaf, alpha, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

        interpolate!(alphaf_HO, alpha, config)
        interpolate!(∇alphaf_HO, ∇alpha.result, config)
        interpolate_upwind!(alphaf_upwind, alpha, mdotf, config)
        interpolate_upwind!(∇alphaf_upwind, ∇alpha.result, mdotf, config)

        alpha_explicit!(alpha_prev, alpha, alphaf, mdotf, rho, dt_cpu[1], config, alphaf_upwind, alphaf_HO, ∇alphaf_upwind, ∇alphaf_HO, F_final)

        @. alpha.values = clamp(alpha.values, 0.0, 1.0)
        @. alphaf.values = clamp(alphaf.values, 0.0, 1.0)
        
        # if sub_cycle
        # #     courant = max_courant_number!(cellsCourant, model, config)
        # #     alphaCourant = max_alpha_courant_number!(cellsAlphaCourant, alpha, mdotf, model, config, dt_cpu[1])

        # #     n_sub = ceil(Int, alphaCourant / maxAlphaCo)
        # #     n_sub = clamp(n_sub, 1, 10) # Ensure at least 1 step and max of 10 steps

        # #     println("AlphaCo: $(alphaCourant); maxAlphaCo: $(maxAlphaCo); Co: $courant;  Number of cycles: ($n_sub)")
        # #     sub_dt = dt_cpu[1] / n_sub
        # #     current_time_sum = 0.0

        # #     for sub_step in 1:n_sub
        # #         # Correction for the final step to kill floating point error in dt division
        # #         if sub_step == n_sub
        # #             real_sub_dt = dt_cpu[1] - current_time_sum
        # #         else
        # #             real_sub_dt = sub_dt
        # #         end

        # #         t_sub = time - dt_cpu[1] + current_time_sum + real_sub_dt

        # #         @. alpha_prev.values = alpha.values
        # #         grad!(∇alpha, alphaf, alpha, time, config)
        # #         limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

        # #         interpolate!(alphaf_HO, alpha, config)
        # #         interpolate!(∇alphaf_HO, ∇alpha.result, config)
        # #         interpolate_upwind!(alphaf_upwind, alpha, mdotf, config)
        # #         interpolate_upwind!(∇alphaf_upwind, ∇alpha.result, mdotf, config)
                
        # #         alpha_explicit!(alpha_prev, alpha, alphaf, mdotf, rho, dt_cpu[1], config, alphaf_upwind, alphaf_HO, ∇alphaf_upwind, ∇alphaf_HO, F_final)

        # #         current_time_sum += real_sub_dt
        # #     end
        # # else
        # #     @. alpha_prev.values = alpha.values
        # #     grad!(∇alpha, alphaf, alpha, time, config)
        # #     limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

        # #     interpolate!(alphaf_HO, alpha, config)
        # #     interpolate!(∇alphaf_HO, ∇alpha.result, config)
        # #     interpolate_upwind!(alphaf_upwind, alpha, mdotf, config)
        # #     interpolate_upwind!(∇alphaf_upwind, ∇alpha.result, mdotf, config)

        # #     alpha_explicit!(alpha_prev, alpha, alphaf, mdotf, rho, dt_cpu[1], config, alphaf_upwind, alphaf_HO, ∇alphaf_upwind, ∇alphaf_HO, F_final)

        # #     # @. alpha.values = clamp(alpha.values, 0.0, 1.0)
        # #     # @. alphaf.values = clamp(alphaf.values, 0.0, 1.0)
        # # end

        grad!(∇alpha, alphaf, alpha, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)
        interpolate!(∇alphaf_HO, ∇alpha.result, config)
        interpolate_upwind!(∇alphaf_upwind, ∇alpha.result, mdotf, config)

        # ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config, rho_prev; time=time)
        # interpolate_upwind!(alphaf, alpha, mdotf, config)

        @. rho_prev.values = rho.values

        blend_properties!(rho, alpha, rho1, rho2)
        blend_properties!(nu, alpha, phases[1].nu, phases[2].nu)

        blend_properties_at_faces!(rhof, alphaf, phases[1].density.rho, phases[2].density.rho)
        blend_properties_at_faces!(nuf, alphaf, phases[1].mu.mu / phases[1].density.rho, phases[2].mu.mu / phases[2].density.rho)

        interpolate!(rho1f, rho1, config)
        interpolate!(rho2f, rho2, config)

        grad!(∇rho, rhof, rho, time, config)
        limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

        @. rhoPhi.values = F_final.values * (rho1f.values - rho2f.values) + mdotf.values * rho2f.values
        
        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config, rho_prev; time=time)

        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config)
        interpolate!(rDf, rD, config)

        remove_pressure_source!(U_eqn, ∇p_rgh, config)
        
        rp = 0.0
        for i ∈ 1:inner_loops
            H!(Hv, U, U_eqn, config)
            
            interpolate!(Uf, Hv, config)
            correct_boundaries!(Uf, Hv, boundaries.U, time, config)

            flux!(mdotf, Uf, config)

            phi_gf!(phi_gf, rho, ghf, rDf, model, config)

            nhat_prep!(nhatf_prep, alpha, ∇alphaf_HO, config) #must be HO
            div!(kappa, nhatf_prep, config)
            interpolate!(kappaf, kappa, config)
            surface_tension_flux!(rDf, sigma, kappaf, alpha, ∇alphaf_HO, phi_gf, config)
            
            
            reconstruct_operation!(phi_g, phi_gf, config)
            @. mdotf.values += phi_gf.values

            div!(divHv, mdotf, config)
            
            @. prev = p_rgh.values
            rp = solve_equation!(p_eqn, p_rgh, boundaries.p_rgh, solvers.p_rgh, config, rho_prev; ref=0.0, time=time)

            grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config) 
            limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)
            
            # correct_mass_flux1(mdotf, p_rgh, rDf, config)
            correct_mass_flux(mdotf, p_eqn, config)
            
            pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)
            reconstruct_operation!(∇p_rghf_reconstructed, ∇p_rghf_deconstructed, config) # reconstruction for velocity correction
            correct_velocity_rgh!(U, Hv, ∇p_rghf_reconstructed, rD, phi_g, config)
        end # corrector loop end
        
    # @. rhoPhi.values = mdotf.values * rhof.values
    @. p.values = p_rgh.values + (rho.values * gh.values)

    turbulence!(turbulenceModel, model, S, prev, time, config)
    update_nueff!(nueff, nuf, model.turbulence, config)

    courant = max_courant_number!(cellsCourant, model, config)
    alphaCourant = max_alpha_courant_number!(cellsAlphaCourant, alpha, mdotf, model, config, dt_cpu[1])
    
    update_dt!(config.runtime, courant, alphaCourant)
    # update_dt!(config.runtime, courant)

    R_ux[iteration] = rx
    R_uy[iteration] = ry
    R_uz[iteration] = rz
    R_p[iteration] = rp
    # R_alpha[iteration] = ralpha

    ProgressMeter.next!(
        progress, showvalues = [
            (:dt, dt_cpu[1]),
            (:time, time),
            (:Courant, courant),
            (:AlphaCourant, alphaCourant),
            (:Ux, R_ux[iteration]),
            (:Uy, R_uy[iteration]),
            (:Uz, R_uz[iteration]),
            (:p_rgh, R_p[iteration]),
            # (:alpha, R_alpha[iteration]),
            turbulenceModel.state.residuals...
            ]
        )

    runtime_postprocessing!(postprocess,iteration,iterations)
    
    if iteration%write_interval + signbit(write_interval) == 0
        save_output(model, outputWriter, iteration, time, config)
        save_postprocessing(postprocess,iteration,time,mesh,outputWriter,config.boundaries)
    end

    end # end for loop
    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p)
end


function correct_mass_flux1(mdotf, p, rDf, config)
    # sngrad = FaceScalarField(mesh)
    (; faces, cells, boundary_cellsID) = mdotf.mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_ifaces # length(n_ifaces) was a BUG! should be n_ifaces only!!!!
    kernel! = _correct_mass_flux1(_setup(backend, workgroup, ndrange)...)
    kernel!(mdotf, p, rDf, faces, cells, n_bfaces)
    # KernelAbstractions.synchronize(backend)
end

@kernel function _correct_mass_flux1(mdotf, p, rDf, faces, cells, n_bfaces)
    i = @index(Global)
    fID = i + n_bfaces

    @inbounds begin 
        face = faces[fID]
        (; area, normal, ownerCells, delta) = face 
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        p1 = p[cID1]
        p2 = p[cID2]
        face_grad = area*(p2 - p1)/delta # best option so far!
        mdotf[fID] -= face_grad*rDf[fID]
    end
end



function nhat_prep!(nhatf_prep, alpha, ∇alphaf, config)
    (; hardware) = config
    (; faces, cells, boundary_cellsID) = alpha.mesh
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(nhatf_prep.x) # ?
    kernel! = _nhat_prep!(_setup(backend, workgroup, ndrange)...)
    kernel!(nhatf_prep, alpha, faces, ∇alphaf)
end
@kernel inbounds=true function _nhat_prep!(nhatf_prep, alpha, faces, ∇alphaf_)
    i = @index(Global)
    face = faces[i]
    (; area, normal, ownerCells, delta) = face

    nhatf_prep[i] = ∇alphaf_[i]/(norm(∇alphaf_[i])+eps())

    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    alpha1 = alpha[cID1]
    alpha2 = alpha[cID2]

    # ∇alphaf = normal * ((alpha2 - alpha1)/delta)
    # nhatf_prep[i] = ∇alphaf/(norm(∇alphaf)+eps())
end


function surface_tension_flux!(rDf, sigma, kappaf, alpha, ∇alphaf_HO, phi_gf, config)
    (; hardware) = config
    (; faces, cells, boundary_cellsID) = phi_gf.mesh
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(phi_gf)
    kernel! = _surface_tension_flux!(_setup(backend, workgroup, ndrange)...)
    kernel!(rDf, sigma, kappaf, alpha, ∇alphaf_HO, phi_gf, faces)
end
@kernel inbounds=true function _surface_tension_flux!(rDf, sigma, kappaf, alpha, ∇alphaf_HO, phi_gf, faces)
    i = @index(Global)
    
    mesh = kappaf.mesh
    face = mesh.faces[i]
    (; area, normal, ownerCells, delta) = face
    Sf = area*normal

    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    alpha1 = alpha[cID1]
    alpha2 = alpha[cID2]

    ∇alphaf = normal * ((alpha2 - alpha1)/delta)
    

    phi_gf[i] -= sigma * kappaf[i] * (∇alphaf ⋅ Sf) * rDf[i]

    
        # face_grad = area*(rho2 - rho1)/delta

        # phi_gf[fID] = -ghf[fID] * face_grad * rDf[fID]
end


function alpha_explicit!(alpha_prev, alpha, alphaf, mdotf, rho, dt, config, alphaf_upwind, alphaf_HO, ∇alphaf_upwind, ∇alphaf_HO, F_final)
    (; hardware) = config
    mesh = rho.mesh
    backend = hardware.backend
    workgroup = hardware.workgroup
    cells = rho.mesh.cells

    nbfaces = length(mesh.boundary_cellsID)


    divergence_result = ScalarField(mesh)
    alpha_up = ScalarField(mesh)
    alpha_corr = ScalarField(mesh)
    lambda = ScalarField(mesh)

    lambdaf = FaceScalarField(mesh)
    F_corr = FaceScalarField(mesh)
    F_upwind = FaceScalarField(mesh)

    F_comp = FaceScalarField(mesh)

    F_compression_upwind = FaceScalarField(mesh)
    F_compression_HO = FaceScalarField(mesh)

    ndrange = length(F_corr)
    kernel! = _compute_F_compression(_setup(backend, workgroup, ndrange)...)
    kernel!(F_compression_upwind, F_compression_HO, ∇alphaf_upwind, ∇alphaf_HO, alphaf_upwind, alphaf_HO, mdotf, nbfaces)

    ndrange = length(F_corr)
    kernel! = _compute_alpha_flux(_setup(backend, workgroup, ndrange)...)
    kernel!(mdotf, alphaf_upwind, alphaf_HO, F_corr, F_upwind, F_compression_upwind, F_compression_HO)

    ndrange = length(divergence_result)
    kernel! = _alpha_divergence(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf_upwind, F_upwind)
    ndrange = nbfaces
    kernel! = _alpha_divergence_boundaries(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf_upwind, F_upwind)
    @. alpha_up.values = alpha_prev.values - (dt) * divergence_result.values

    ndrange = length(divergence_result)
    kernel! = _alpha_divergence(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf_upwind, F_corr)
    ndrange = nbfaces
    kernel! = _alpha_divergence_boundaries(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf_upwind, F_corr)
    @. alpha_corr.values = - (dt) * divergence_result.values

    ndrange = length(lambda)
    kernel! = _lambda_calc(_setup(backend, workgroup, ndrange)...)
    kernel!(lambda, alpha_up, alpha_corr)

    ndrange = length(lambdaf)
    kernel! = _lambdaf_calc(_setup(backend, workgroup, ndrange)...)
    kernel!(lambda, lambdaf, F_final, F_upwind, F_corr, nbfaces)

    # ndrange = length(lambdaf)
    # kernel! = _test_flux(_setup(backend, workgroup, ndrange)...)
    # kernel!(n_bfaces, mdotf, alphaf_upwind, F_final)

    ndrange = length(divergence_result)
    kernel! = _alpha_divergence(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf, F_final)
    ndrange = nbfaces
    kernel! = _alpha_divergence_boundaries(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf_upwind, F_final)

    @. alpha.values = alpha_prev.values - (dt) * divergence_result.values
    @. alphaf.values = alphaf_upwind.values + lambdaf.values * (alphaf_HO.values - alphaf_upwind.values)

    

    ndrange = length(divergence_result)
    kernel! = _alpha_divergence(_setup(backend, workgroup, ndrange)...)
    kernel!(divergence_result, mdotf, alphaf, F_compression_HO)
    ndrange = nbfaces
end
@kernel inbounds=true function _test_flux(n_bfaces, mdotf, alphaf_upwind, F_final)
    fID = @index(Global)
    i = fID + n_bfaces
    (; area, normal, ownerCells) = faces[i]

    mesh = mdotf.mesh
    owner_P = F_final[ownerCells[1]]
    owner_N = F_final[ownerCells[2]]
end
@kernel inbounds=true function _compute_alpha_flux(mdotf, alphaf_upwind, alphaf_HO, F_corr, F_upwind, F_compression_upwind, F_compression_HO)
    i = @index(Global)

    mesh = mdotf.mesh

    F_upwind[i] = (mdotf[i] * alphaf_upwind[i]) + F_compression_upwind[i]
    F_HO = (mdotf[i] * alphaf_HO[i]) + F_compression_HO[i]
    F_corr[i] = F_HO - F_upwind[i]
end
@kernel inbounds=true function _compute_F_compression(F_compression_upwind, F_compression_HO, ∇alphaf_upwind, ∇alphaf_HO, alphaf_upwind, alphaf_HO, mdotf, n_bfaces)
    i = @index(Global)

    mesh = mdotf.mesh
    faces = mesh.faces

    (; area, normal) = faces[i]
    TF = eltype(F_compression_upwind)

    # if i <= n_bfaces
    #     F_compression_upwind[i] = zero(TF)
    #     F_compression_HO[i] = zero(TF)
    # else

    C_alpha = 1.0
    Sf = normal * area

    phic_upwind = C_alpha * abs(mdotf[i]) / (area + eps(TF))
    phic_HO = C_alpha * abs(mdotf[i]) / (area + eps(TF))

    n_hat_upwind = ∇alphaf_upwind[i] / (norm(∇alphaf_upwind[i]) + 1e-8)
    n_hat_HO = ∇alphaf_HO[i] / (norm(∇alphaf_HO[i]) + 1e-8)

    phi_r_upwind = phic_upwind * (n_hat_upwind ⋅ Sf)
    phi_r_HO = phic_HO * (n_hat_upwind ⋅ Sf)  # use HO normal for HO flux

    F_compression_upwind[i] = phi_r_upwind * alphaf_upwind[i] * (1.0 - alphaf_upwind[i])
    F_compression_HO[i] = phi_r_HO * alphaf_HO[i] * (1.0 - alphaf_HO[i])
    # end
end
@kernel inbounds=true function _alpha_divergence(divergence_result, mdotf, alphaf, flux)
    i = @index(Global)

    mesh = mdotf.mesh
    cells = mesh.cells
    cell_faces = mesh.cell_faces
    cell_nsign = mesh.cell_nsign

    T = eltype(divergence_result.values)
    sum = zero(T)
    volume = cells[i].volume

    fr = cells[i].faces_range
    @inbounds for k in fr
        pointer = cell_faces[k]
        sum += flux[pointer] * cell_nsign[k]
    end

    divergence_result[i] = sum / volume
end
@kernel inbounds=true function _alpha_divergence_boundaries(divergence_result, mdotf, alphaf, flux)
    i = @index(Global)

    mesh = mdotf.mesh
    cells = mesh.cells
    faces = mesh.faces
    cID = faces[i].ownerCells[1]
    volume = cells[cID].volume

    divergence_result[cID] += flux[i] / volume # must be here! otherwise it goes bad!
end
@kernel inbounds=true function _lambda_calc(lambda, alpha_up, alpha_corr)
    i = @index(Global)

    # we will do ifelse just for testing

    increase_term = (1.0 - alpha_up[i])/(alpha_corr[i]+eps())
    decrease_term = alpha_up[i]/(-alpha_corr[i]+eps())

    if abs(alpha_corr[i]) < 1.0e-8
        lambda[i] = 1.0
    elseif alpha_corr[i] > 0.0
        lambda[i] = min(1, increase_term)
    elseif alpha_corr[i] < 0.0
        lambda[i] = min(1, decrease_term)
    end

    lambda[i] = clamp(lambda[i], 0.0, 1.0)
end
@kernel inbounds=true function _lambdaf_calc(lambda, lambdaf, F_final, F_upwind, F_corr, n_bfaces) #questionable at boundaries
    i = @index(Global)
    mesh = lambda.mesh
    face = mesh.faces[i]

    (; ownerCells) = face

    if i <= n_bfaces
        lambdaf[i] = 1.0
        F_final[i] = F_upwind[i]
    else
        lambda_P = lambda[ownerCells[1]]
        lambda_N = lambda[ownerCells[2]]
        lambdaf[i] = min(lambda_P, lambda_N)
        F_final[i] = F_upwind[i] + lambdaf[i]*F_corr[i]
    end
end



function update_phase_thermodynamics!(EoS::AbstractEosModel, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    return nothing
end

function update_phase_thermodynamics!(EoS::ConstEos, phaseIndex::Val{N}, nueff, T, model, config) where {N} #EoS::Union{ConstEos, PerfectGas}
    phase = model.fluid.phases[N]
    phase.density(phase, model, config)
    phase.mu(phase, model)
end


## VELOCITY CORRECTION

function correct_velocity_rgh!(U, Hv, ∇p, rD, phi_g, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(U)
    kernel! = _correct_velocity_rgh!(_setup(backend, workgroup, ndrange)...)
    kernel!(U, Hv, ∇p, rD, phi_g)
    # # KernelAbstractions.synchronize(backend)
end

@kernel function _correct_velocity_rgh!(U, Hv, ∇p, rD, phi_g)
    i = @index(Global)

    @uniform begin
        Ux, Uy, Uz = U.x, U.y, U.z
        Hvx, Hvy, Hvz = Hv.x, Hv.y, Hv.z
        dpdx, dpdy, dpdz = ∇p.x, ∇p.y, ∇p.z
        phi_gx, phi_gy, phi_gz = phi_g.x, phi_g.y, phi_g.z
        rDvalues = rD.values
    end

    @inbounds begin
        rDvalues_i = rDvalues[i]

        Ux[i] = Hvx[i] + dpdx[i] * rDvalues_i
        Uy[i] = Hvy[i] + dpdy[i] * rDvalues_i
        Uz[i] = Hvz[i] + dpdz[i] * rDvalues_i
    end
end


# FaceScalarField into VectorField
function reconstruct_operation!(phi::VectorField, psif::FaceScalarField, config)
    mesh = phi.mesh
    (; cells, cell_nsign, cell_faces, faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    F = _get_float(mesh)

    # Launch main calculation kernel
    ndrange = length(cells)
    kernel! = _reconstruct_operation_internal!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, F, cell_faces, cell_nsign, faces, phi, psif)
end


@kernel function _reconstruct_operation_internal!(
    cells::AbstractArray{Cell{TF,SV,UR}}, F, cell_faces, cell_nsign, faces, phi, psif
) where {TF,SV,UR}

    i = @index(Global)

    @inbounds begin
        (; faces_range) = cells[i]

        m11 = zero(TF); m12 = zero(TF); m22 = zero(TF)
        b1  = zero(TF); b2  = zero(TF)

        for fi ∈ faces_range
            fID = cell_faces[fi]
            # nsign is optional IF psif is consistent with faces[fID].normal orientation (it seems like it is!)

            (; area, normal) = faces[fID]
            nx = normal[1]; ny = normal[2]

            m11 += area * nx * nx
            m12 += area * nx * ny
            m22 += area * ny * ny

            # b += n * ssf
            ssf = psif[fID]
            b1 += nx * ssf
            b2 += ny * ssf
        end

        det = m11*m22 - m12*m12
        invdet = one(TF) / (det + eps(TF))  # eps guard

        ux = ( m22*b1 - m12*b2) * invdet
        uy = (-m12*b1 + m11*b2) * invdet

        phi[i] = @SVector [ux, uy, zero(TF)]
    end
end

function pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    faces = ∇p_rghf_deconstructed.mesh.faces

    ndrange = length(∇p_rghf_deconstructed)
    kernel! = _pressure_grad!(_setup(backend, workgroup, ndrange)...)
    kernel!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, faces)
end
@kernel function _pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, faces)
    i = @index(Global)
    face = faces[i]
    (; area, normal, ownerCells, delta) = face

    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    p1 = p_rgh[cID1]
    p2 = p_rgh[cID2]
    face_grad = area*(p2 - p1)/delta

    ∇p_rghf_deconstructed[i] = (phi_gf[i] - (face_grad*rDf[i]))  / (rDf[i] + eps())
end


"""
    blend_properties!(property_field, alpha_field, property_0, property_1)

Blends a property between two phases using the phase fraction `alpha_field`.
Formula: `prop = (prop0 * alpha) + (prop1 * (1 - alpha))`
"""
function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0.values * alpha_field.values) + (property_1.values * (1.0 - alpha_field.values))
    nothing
end

function blend_properties_at_faces!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0 * alpha_field.values) + (property_1 * (1.0 - alpha_field.values))
    nothing
end

"""
    compute_gh!(gh, g, config)

Computes the dot product of gravity vector and cell centres: `g . x`.
Used for hydrostatic pressure reconstruction.
"""
function compute_gh!(gh, g, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    cells = gh.mesh.cells

    ndrange = length(gh)
    kernel! = _compute_gh!(_setup(backend, workgroup, ndrange)...)
    kernel!(gh, g, cells)
end
@kernel inbounds=true function _compute_gh!(gh, g, cells)
    i = @index(Global)
    (; centre) = cells[i]

    gh[i] = (g ⋅ (centre))
end

"""
    compute_ghf!(ghf, g, config)

Computes the dot product of gravity vector and face centres: `g . x_f`.
"""
function compute_ghf!(ghf, g, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    faces = ghf.mesh.faces

    ndrange = length(ghf)
    kernel! = _compute_ghf!(_setup(backend, workgroup, ndrange)...)
    kernel!(ghf, g, faces)
end
@kernel inbounds=true function _compute_ghf!(ghf, g, faces)
    i = @index(Global)
    (; centre) = faces[i]

    ghf[i] = (g ⋅ (centre))
end

"""
    compute_p_rgh!(p_rgh, gh, p, rho, config)

Computes dynamic pressure `p_rgh` from absolute pressure `p` and hydrostatic head.
Formula: `p_rgh = p - rho * (g . x)`
"""
function compute_p_rgh!(p_rgh, gh, p, rho, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(p_rgh)
    kernel! = _compute_p_rgh!(_setup(backend, workgroup, ndrange)...)
    kernel!(p_rgh, gh, p, rho)
end
@kernel inbounds=true function _compute_p_rgh!(p_rgh, gh, p, rho)
    i = @index(Global)

    p_rgh[i] = p[i] - (rho[i] * gh[i])
end

"""
    compute_p_rghf!(p_rghf, ghf, pf, rhof, config)

Computes dynamic pressure at faces.
"""
function compute_p_rghf!(p_rghf, ghf, pf, rhof, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(p_rghf)
    kernel! = _compute_p_rghf!(_setup(backend, workgroup, ndrange)...)
    kernel!(p_rghf, ghf, pf, rhof)
end
@kernel inbounds=true function _compute_p_rghf!(p_rghf, ghf, pf, rhof)
    i = @index(Global)

    p_rghf[i] = pf[i] - (rhof[i] * ghf[i])
end

"""
    phi_gf!(phi_gf, rho, ghf, rDf, model, config)

Computes the gravity contribution to the face flux for pressure-velocity coupling.
"""
function phi_gf!(phi_gf, rho, ghf, rDf, model, config)
    (; faces, cells, boundary_cellsID) = model.domain
    (; hardware) = config
    (; backend, workgroup) = hardware

    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_faces
    kernel! = _phi_gf!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_gf, rho, ghf, rDf, faces, cells, n_bfaces)
end

@kernel function _phi_gf!(phi_gf, rho, ghf, rDf, faces, cells, n_bfaces)
    fID = @index(Global)

    @inbounds begin 
        face = faces[fID]
        (; area, normal, ownerCells, delta) = face
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        rho1 = rho[cID1]
        rho2 = rho[cID2]

        face_grad = area*(rho2 - rho1)/delta

        phi_gf[fID] = -ghf[fID] * face_grad * rDf[fID]
    end
end