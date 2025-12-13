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

    (; solvers, schemes, runtime, hardware, boundaries, postprocess) = config

    @info "Extracting configuration and input fields..."

    (; U, p, Uf, pf) = model.momentum
    (; nu, nuf, rho, rhof, alpha, alphaf, p_rgh, p_rghf) = model.fluid
    mesh = model.domain
    
    phases = model.fluid.phases
    phase_eos = [phases[1].eosModel, phases[2].eosModel]

    @info "Pre-allocating fields..."
    
    ∇p = Grad{schemes.p.gradient}(p)
    mdotf = FaceScalarField(mesh)
    # rhorDf = FaceScalarField(mesh)
    # initialise!(rhorDf, 1.0)
    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)
    mueff = FaceScalarField(mesh)
    mueffgradUt = VectorField(mesh)
    divHv = ScalarField(mesh)
    ddtrho = ScalarField(mesh)
    psidpdt = ScalarField(mesh)
    divmdotf = ScalarField(mesh)
    psi = ScalarField(mesh)

    phi_g = VectorField(mesh)
    phi_gf = FaceScalarField(mesh)

    update_phase_thermodynamics!(phase_eos[1], Val(1), 0.0, ConstantScalar(300.0), model, config)
    update_phase_thermodynamics!(phase_eos[2], Val(2), 0.0, ConstantScalar(300.0), model, config)

    blend_properties!(rho, alpha, phases[1].rho, phases[2].rho) # mixture.correct(); in OpenFOAM
    blend_viscosity!(alpha, phases, nu) # mixture.correct() in OpenFOAM

    interpolate!(rhof, rho, config)
    interpolate!(nuf, nu, config)
    
    gh = model.fluid.physics_properties.gravity.gh   # gh field
    ghf = model.fluid.physics_properties.gravity.ghf # ghf field
    g = model.fluid.physics_properties.gravity.g     # gravity vector


    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)
    compute_p_rgh!(p_rgh, gh, p, rho, config)
    compute_p_rghf!(p_rghf, ghf, pf, rhof, config)
    
    ∇p_rgh = Grad{schemes.p_rgh.gradient}(p_rgh)
    grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)
    
    ∇alpha = Grad{schemes.alpha.gradient}(alpha)
    grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
    limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)
    
    ∇rho = Grad{schemes.p_rgh.gradient}(rho)
    grad!(∇rho, rhof, rho, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

    phi_g!(phi_g, gh, ∇rho, config)

    divCompressionFlux = ScalarField(mesh)
    compressionFlux = FaceVectorField(mesh)
    
    nhat = VectorField(mesh)
    nhatf = FaceVectorField(mesh)
    U_c = FaceVectorField(mesh)

    compute_compression_flux_step1!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha, config)
    interpolate!(nhatf, nhat, config)
    compute_compression_flux_step2!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha, config)
    div!(divCompressionFlux, compressionFlux, config)


    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(mdotf, U) 
        - Laplacian{schemes.U.laplacian}(mueff, U) 
        == 
        - Source(∇p_rgh.result)
        + Source(phi_g)
    ) → VectorEquation(U, boundaries.U)

    p_eqn = (
        # Time{schemes.p.time}(psi, p)
        - Laplacian{schemes.p.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)

    alpha_eqn = (
        Time{schemes.alpha.time}(rho, alpha)
        + Divergence{schemes.alpha.divergence}(mdotf, alpha)
        == 
        - Source(divCompressionFlux)
        # Source(ConstantScalar(0.0))
    ) → ScalarEquation(alpha, boundaries.alpha)

    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p_rgh.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."
     
    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p_rgh.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))

    # @info "Initialising energy model..."
    # energyModel = initialise(model.energy, model, mdotf, rho, p_eqn, config)

    # @info "Initialising turbulence model..."
    # turbulenceModel, config = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        # model, ∇p, ∇p_rgh, ∇rho, ∇alpha, U_eqn, p_eqn, alpha_eqn, gh, ghf, phi_g, phi_gf, config;
        model, ∇p, ∇p_rgh, ∇rho, ∇alpha, U_eqn, p_eqn, alpha_eqn, gh, ghf, phi_g, phi_gf, compressionFlux, divCompressionFlux, nhat, nhatf, U_c, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops)

    return residuals    
end # end function

function MULTIPHASE(
    model, ∇p, ∇p_rgh, ∇rho, ∇alpha, U_eqn, p_eqn, alpha_eqn, gh, ghf, phi_g, phi_gf, compressionFlux, divCompressionFlux, nhat, nhatf, U_c, config; 
    # model, ∇p, ∇p_rgh, ∇rho, ∇alpha, U_eqn, p_eqn, alpha_eqn, gh, ghf, phi_g, phi_gf, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    )
    
    (; U, p, Uf, pf) = model.momentum
    (; nu, nuf, rho, rhof, alpha, alphaf, p_rgh, p_rghf) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries, postprocess) = config
    (; iterations, write_interval, dt) = runtime
    (; backend) = hardware

    phases = model.fluid.phases
    
    postprocess = convert_time_to_iterations(postprocess,model,dt,iterations)
    mdotf = get_flux(U_eqn, 2)
    # nueff = get_flux(U_eqn, 3)
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)

    @info "Allocating working memory..."

    # Define aux fields 
    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    Uf = FaceVectorField(mesh)
    S = StrainRate(gradU, gradUT, U, Uf)

    n_cells = length(mesh.cells)
    pf = FaceScalarField(mesh)
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
    
    # Initial calculations
    time = zero(TF) # assuming time=0
    interpolate!(Uf, U, config)   
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    # flux!(mdotf, Uf, rhof, config) 
    flux!(mdotf, Uf, config) #OLD ONE


    # grad!(∇p, pf, p, boundaries.p, time, config)
    # limit_gradient!(schemes.p.limiter, ∇p, p, config)

    # grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
    # limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)

    phase_eos = [phases[1].eosModel, phases[2].eosModel]

    # update_nueff!(nueff, nu, model.turbulence, config)

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting MULTIPHASE loops..."

    progress = Progress(iterations; dt=1.0, showspeed=true)


    @time for iteration ∈ 1:iterations
        time = iteration *dt

        update_phase_thermodynamics!(phase_eos[1], Val(1), 0.0, ConstantScalar(300.0), model, config)
        update_phase_thermodynamics!(phase_eos[2], Val(2), 0.0, ConstantScalar(300.0), model, config)

        blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)
        blend_viscosity!(alpha, phases, nu)
        interpolate!(rhof, rho, config)
        interpolate!(nuf, nu, config)

        blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)
        blend_viscosity!(alpha, phases, nu) 

        interpolate!(rhof, rho, config)
        interpolate!(nuf, nu, config)

        grad!(∇rho, rhof, rho, time, config)

        # println("SOLVING U_EQN")
        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config, rhof, p_eqn, boundaries.p_rgh; time=time) # rhof, p_eqn


        ∇alpha = Grad{schemes.alpha.gradient}(alpha)
        grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)
        

        compute_compression_flux_step1!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha, config)
        interpolate!(nhatf, nhat, config)
        compute_compression_flux_step2!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha, config)
        div!(divCompressionFlux, compressionFlux, config)
        

        ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config, rhof, p_eqn, boundaries.U; time=time)
        interpolate!(alphaf, alpha, config)
        correct_boundaries!(alphaf, alpha, boundaries.alpha, time, config)
          
        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config)
        interpolate!(rDf, rD, config)
        remove_pressure_source!(U_eqn, ∇p_rgh, config)
        
        rp = 0.0
        for i ∈ 1:inner_loops
            H!(Hv, U, U_eqn, config)
            
            interpolate!(Uf, Hv, config) # Careful: reusing Uf for interpolation
            correct_boundaries!(Uf, Hv, boundaries.U, time, config)

            # flux!(mdotf, Uf, rhof, config)
            flux!(mdotf, Uf, config)

            phi_g!(phi_g, gh, ∇rho, config)
            phi_gf!(phi_gf, rho, ghf, rDf, model, config)
            @. mdotf.values += phi_gf.values

            div!(divHv, mdotf, config)
            
            @. prev = p_rgh.values
            # println("SOLVING P_rgh_EQN")
            rp = solve_equation!(p_eqn, p_rgh, boundaries.p_rgh, solvers.p_rgh, config, rhof, U_eqn, boundaries.U; ref=pref, time=time)
            if i == inner_loops
                explicit_relaxation!(p_rgh, prev, 1.0, config)
            else
                explicit_relaxation!(p_rgh, prev, solvers.p_rgh.relax, config)
            end

            grad!(∇p_rgh, pf, p_rgh, boundaries.p_rgh, time, config) 
            limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)

            # for i ∈ 1:ncorrectors
            #     discretise!(p_eqn, p_rgh, config)       
            #     apply_boundary_conditions!(p_eqn, boundaries.p, nothing, time, config)
            #     setReference!(p_eqn, pref, 1, config)
            #     nonorthogonal_face_correction(p_eqn, ∇p, rDf, config)
            #     update_preconditioner!(p_eqn.preconditioner, p.mesh, config)
            #     rp = solve_system!(p_eqn, solvers.p, p, nothing, config)

            #     if i == ncorrectors
            #         explicit_relaxation!(p, prev, 1.0, config)
            #     else
            #         explicit_relaxation!(p, prev, solvers.p.relax, config)
            #     end
            #     grad!(∇p, pf, p, boundaries.p, time, config) 
            #     limit_gradient!(schemes.p.limiter, ∇p, p, config)
            # end
            
            correct_mass_flux(mdotf, p_rgh, rDf, config)
            # correct_velocity!(U, Hv, ∇p_rgh, rD, config)
            #### phi_g is not updated!!!!!!!!!!!!!!!!!!

            correct_velocity_multiphase_test!(U, Hv, ∇p_rgh, rD, phi_g, config)

        end # corrector loop end
        
    @. p.values = p_rgh.values + (rho.values * gh.values)

    # if iteration % 500 == 0
    #     mean_alpha = sum(alpha.values) / length(alpha.values)
    #     @info "Mean alpha at iteration $iteration: $mean_alpha"
    # end

    # @info "Iteration"


    # turbulence!(turbulenceModel, model, S, prev, time, config) 
    # update_nueff!(nueff, nu, model.turbulence, config)

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
            # turbulenceModel.state.residuals...
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


function update_phase_thermodynamics!(EoS::AbstractEosModel, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    return nothing
end

function update_phase_thermodynamics!(EoS::Union{ConstEos, PerfectGas}, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    phase = model.fluid.phases[N]
    phase.eosModel(phase, model, config)
    phase.viscosityModel(phase, T)
end


function blend_viscosity!(alpha, phases, nu) # for viscosity specifically (does the conersion back-n-forth)
    mu_1 = phases[1].mu
    mu_2 = phases[2].mu

    rho_1 = phases[1].rho
    rho_2 = phases[2].rho

    @. mu_1.values = mu_1.values / rho_1.values # nu
    @. mu_2.values = mu_2.values / rho_2.values # nu

    blend_properties!(nu, alpha, mu_1, mu_2)

    @. mu_1.values = mu_1.values * rho_1.values # mu
    @. mu_2.values = mu_2.values * rho_2.values # mu
end

function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0.values * alpha_field.values) + (property_1.values * (1.0 - alpha_field.values))
    nothing
end


function compute_compression_flux_step1!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(nhat)
    kernel! = _compute_compression_flux_step1!(_setup(backend, workgroup, ndrange)...)
    kernel!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha)
end
@kernel inbounds=true function _compute_compression_flux_step1!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha)
    i = @index(Global)

    nhat[i] = ∇alpha[i] / (norm(∇alpha[i]) + eps())
end


function compute_compression_flux_step2!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(nhatf)
    kernel! = _compute_compression_flux_step2!(_setup(backend, workgroup, ndrange)...)
    kernel!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha)
end
@kernel inbounds=true function _compute_compression_flux_step2!(nhat, nhatf, compressionFlux, U_c, Uf, alphaf, ∇alpha)
    i = @index(Global)

    U_c[i] = 1.0 * nhatf[i] * norm(Uf[i])

    compressionFlux[i] = alphaf[i] * (1.0 - alphaf[i]) * U_c[i]
end



function compute_gh!(gh, g, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(gh)
    kernel! = _compute_gh!(_setup(backend, workgroup, ndrange)...)
    kernel!(gh, g)
end
@kernel inbounds=true function _compute_gh!(gh, g)
    i = @index(Global)

    gh[i] = dot(g, gh.mesh.cells[i].centre)
end


function compute_ghf!(ghf, g, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(ghf)
    kernel! = _compute_ghf!(_setup(backend, workgroup, ndrange)...)
    kernel!(ghf, g)
end
@kernel inbounds=true function _compute_ghf!(ghf, g)
    i = @index(Global)

    ghf[i] = dot(g, ghf.mesh.faces[i].centre)
end


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


function phi_g!(phi_g, gh, ∇rho, config)  # CELL-CENTRED FOR `U_EQN`
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(phi_g)
    kernel! = _phi_g!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_g, gh, ∇rho)
end
@kernel inbounds=true function _phi_g!(phi_g, gh, ∇rho)
    i = @index(Global)

    phi_g[i] = -gh[i] * ∇rho.result[i] # VECTOR FIELD (CELL-CENTRED)
end


function phi_gf!(phi_gf, rho, ghf, rDf, model, config)
    (; faces, cells, boundary_cellsID) = model.domain
    (; hardware) = config
    (; backend, workgroup) = hardware

    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_ifaces
    kernel! = _phi_gf!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_gf, rho, ghf, rDf, faces, cells, n_bfaces, model)
end

@kernel function _phi_gf!(phi_gf, rho, ghf, rDf, faces, cells, n_bfaces, model)
    i = @index(Global)
    fID = i + n_bfaces

    @inbounds begin 
        face = faces[fID]
        (; area, normal, ownerCells, delta) = face 
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        rho1 = rho[cID1]
        rho2 = rho[cID2]

        face_grad = area*(rho2 - rho1)/delta #area

        phi_gf[fID] = -ghf[fID] * face_grad * rDf[fID] # Important line (rDf multiplication)
    end
end







function correct_velocity_multiphase_test!(U, Hv, ∇p_rgh, rD, phi_g, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(U)
    kernel! = _correct_velocity_multiphase_test!(_setup(backend, workgroup, ndrange)...)
    kernel!(U, Hv, ∇p_rgh, rD, phi_g)
end

@kernel function _correct_velocity_multiphase_test!(U, Hv, ∇p_rgh, rD, phi_g)
    i = @index(Global)

    @uniform begin
        Ux, Uy, Uz = U.x, U.y, U.z
        Hvx, Hvy, Hvz = Hv.x, Hv.y, Hv.z
        dpdx, dpdy, dpdz = ∇p_rgh.result.x, ∇p_rgh.result.y, ∇p_rgh.result.z
        phi_g_x, phi_g_y, phi_g_z = phi_g.x, phi_g.y, phi_g.z
        rDvalues = rD.values
    end

    @inbounds begin
        rD_i = rDvalues[i]

        x_diff = (phi_g_x[i] - dpdx[i])
        y_diff = (phi_g_y[i] - dpdy[i])
        z_diff = (phi_g_z[i] - dpdz[i])

        # x_diff = (dpdx[i])
        # y_diff = (dpdy[i])
        # z_diff = (dpdz[i])

        Ux[i] = Hvx[i] + x_diff * rD_i
        Uy[i] = Hvy[i] + y_diff * rD_i
        Uz[i] = Hvz[i] + z_diff * rD_i
    end
end



# function correct_mass_flux_multiphase2(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, phi_gf, phi_force, config, model)
#     (; faces, cells, boundary_cellsID) = rhof.mesh
#     (; hardware) = config
#     (; backend, workgroup) = hardware

#     n_faces = length(faces)
#     n_bfaces = length(boundary_cellsID)
#     n_ifaces = n_faces - n_bfaces

#     ndrange = n_ifaces 
#     kernel! = _correct_mass_flux_multiphase(_setup(backend, workgroup, ndrange)...)
#     kernel!(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, phi_gf, phi_force, faces, cells, n_bfaces, model)
# end


# @kernel function _correct_mass_flux_multiphase2(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, phi_gf, phi_force,faces, cells, n_bfaces, model)
#     i = @index(Global)
#     fID = i + n_bfaces

#     @inbounds begin 
#         face = faces[fID]
#         (; area, normal, ownerCells, delta) = face 
#         cID1 = ownerCells[1]
#         cID2 = ownerCells[2]
#         p1 = p_rgh[cID1]
#         p2 = p_rgh[cID2]

#         face_grad = area * (p2-p1) / delta
#         mdotf[fID] -= face_grad*rDf[fID]
#         p_flux = face_grad * rDf[fID]
#         phi[fID] = phiHbyA[fID] - p_flux

#         rhoPhi[fID] = phi[fID] * rhof[fID]
#         phi_force[fID] = (phi_gf[fID] - p_flux)
        
#     end
# end


#     @inbounds begin 
#         face = faces[fID]
#         (; area, normal, ownerCells, delta) = face 
#         cID1 = ownerCells[1]
#         cID2 = ownerCells[2]
#         p1 = p[cID1]
#         p2 = p[cID2]
#         face_grad = area*(p2 - p1)/delta # best option so far!
#         mdotf[fID] -= face_grad*rDf[fID]
#     end