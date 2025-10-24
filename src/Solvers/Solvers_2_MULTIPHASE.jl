export multiphase!

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

    (; U, Uf, p, pf) = model.momentum
    (; alpha, alphaf, rho, rhof, nu, nuf, p_rgh, p_rghf) = model.fluid

    phases = model.fluid.phases

    backend = hardware.backend
    workgroup = hardware.workgroup
    mesh = model.domain

    @info "Pre-allocating fields..."

    TF = _get_float(mesh)
    time = zero(TF) # assuming time=0

    
    # ∇p_rgh = Grad{schemes.p_rgh.gradient}(p_rgh) # Define grad(p_rgh) 
    # grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
    # limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)
    
    # ∇rho = Grad{schemes.p_rgh.gradient}(rho)
    # grad!(∇rho, rhof, rho, time, config) # No BC given for this one
    # limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)


    ∇U = Grad{schemes.U.gradient}(U)
    grad!(∇U, Uf, U, boundaries.U, time, config)
    limit_gradient!(schemes.U.limiter, ∇U, U, config)
    
    ∇alpha = Grad{schemes.alpha.gradient}(alpha)
    grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
    limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

    mdotf = FaceScalarField(mesh)
    mdotf_VOF = FaceScalarField(mesh)

    rhoPhi = FaceScalarField(mesh)
    phi = FaceScalarField(mesh)

    
    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)
    nueff = FaceScalarField(mesh)
    divHv = ScalarField(mesh)

    interpolate!(alphaf, alpha, config)

    rho_l = phases[1].rho
    rho_v = phases[2].rho
    rhof_l = FaceScalarField(mesh)
    interpolate!(rhof_l, rho_l, config)


    @info "Computing Fluid Properties..."

    phase_eos = [phases[1].eosModel, phases[2].eosModel]

    update_phase_thermodynamics!(phase_eos[1], Val(1), nueff, ConstantScalar(300.0), model, config)
    update_phase_thermodynamics!(phase_eos[2], Val(2), nueff, ConstantScalar(300.0), model, config)

    blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)
    blend_viscosity!(alpha, phases, nu)
    interpolate!(rhof, rho, config)
    interpolate!(nuf, nu, config)



    phi_g = VectorField(mesh)
    phi_gf = FaceScalarField(mesh)
    
    gh = model.fluid.physics_properties.gravity.gh   # gh field
    ghf = model.fluid.physics_properties.gravity.ghf # ghf field
    g = model.fluid.physics_properties.gravity.g     # gravity vector

    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)
    compute_p_rgh!(p_rgh, gh, p, rho, config)
    compute_p_rghf!(p_rghf, ghf, pf, rhof, config)
    
    ∇p_rgh = Grad{schemes.p_rgh.gradient}(p_rgh) # Define grad(p_rgh) 
    grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)
    
    ∇rho = Grad{schemes.p_rgh.gradient}(rho)
    grad!(∇rho, rhof, rho, time, config) # No BC given for this one
    limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

    phi_g!(phi_g, gh, ∇rho, config)

    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(rhoPhi, U)  #????????
        - Laplacian{schemes.U.laplacian}(nueff, U) 
        == 
        - Source(∇p_rgh.result)
        - Source(phi_g)
    ) → VectorEquation(U, boundaries.U)

    alpha_eqn = (
        Time{schemes.alpha.time}(alpha)
        + Divergence{schemes.alpha.divergence}(phi, alpha) 
        ==
        Source(ConstantScalar(0.0))
    ) → ScalarEquation(alpha, boundaries.alpha)


    p_eqn = (
        - Laplacian{schemes.p_rgh.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)

    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p_rgh.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."

    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p_rgh.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))

    @info "Initialising turbulence model..."
    turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        model, turbulenceModel, ∇p_rgh, ∇rho, ∇U, ∇alpha, U_eqn, p_eqn, phi, rhoPhi, phi_g, phi_gf, alpha_eqn, rho_l, rho_v, rhof_l, config; 
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops)

    return residuals
end # end function

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

    phi_g[i] = gh[i] * ∇rho.result[i] # VECTOR FIELD (CELL-CENTRED)
end



function phi_gf!(phi_gf, rho, ghf, rDf, model, config)  # FACE-CENTRED FOR `P_EQN`
    mesh = model.domain
    (; faces, boundary_cellsID) = mesh

    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup
    
    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_ifaces
    kernel! = _phi_gf!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_gf, rho, ghf, rDf, faces, n_bfaces)
end
@kernel inbounds=true function _phi_gf!(phi_gf, rho, ghf, rDf, faces, n_bfaces)
    i = @index(Global)
    fID = i + n_bfaces

    face = faces[fID]
    (; area, normal, ownerCells, delta) = face
    
    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    rho1 = rho[cID1]
    rho2 = rho[cID2]

    face_grad = area*(rho2 - rho1)/delta

    phi_gf[fID] = -ghf[fID] * face_grad * rDf[fID]
end



function MULTIPHASE(
    model, turbulenceModel, ∇p_rgh, ∇rho, ∇U, ∇alpha, U_eqn, p_eqn, phi, rhoPhi, phi_g, phi_gf, alpha_eqn, rho_l, rho_v, rhof_l, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    )
    
    # Extract model variables and configuration
    (; U, Uf, p, pf) = model.momentum
    (; alpha, alphaf, rho, rhof, nu, nuf, p_rgh, p_rghf) = model.fluid
    # (; nu) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries) = config
    (; iterations, write_interval) = runtime
    (; backend) = hardware
    

    phases = model.fluid.phases
    props = model.fluid.physics_properties

    backend = hardware.backend
    workgroup = hardware.workgroup

    gh = model.fluid.physics_properties.gravity.gh   # gh field
    ghf = model.fluid.physics_properties.gravity.ghf # ghf field
    phiHbyA = FaceScalarField(mesh) # used for p_eqn divergence as mdotf+phi_g
    
    phase_eos = [phases[1].eosModel, phases[2].eosModel]

    mdotf = get_flux(U_eqn, 2)
    mdotf_VOF = get_flux(alpha_eqn, 2)
    nueff = get_flux(U_eqn, 3)
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)
    
    @info "Allocating working memory..."

    # Define aux fields 
    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    S = StrainRate(gradU, gradUT, U, Uf)

    n_cells = length(mesh.cells)
    Hv = VectorField(mesh)
    Hvf = FaceVectorField(mesh)
    rD = ScalarField(mesh)

    # Pre-allocate auxiliary variables
    TF = _get_float(mesh)
    # prev = _convert_array!(prev, backend) 
    prev = KernelAbstractions.zeros(backend, TF, n_cells) 

    # Pre-allocate vectors to hold residuals 
    R_ux = ones(TF, iterations)
    R_uy = ones(TF, iterations)
    R_uz = ones(TF, iterations)
    R_p = ones(TF, iterations)
    R_alpha = ones(TF, iterations)
    cellsCourant =adapt(backend, zeros(TF, length(mesh.cells)))
    
    # Initial calculations
    time = zero(TF) # assuming time=0
    interpolate!(Uf, U, config)   
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    # flux!(mdotf, Uf, config)

    update_nueff!(nueff, nuf, model.turbulence, config)

    @info "Starting MULTIPHASE loops..."

    progress = Progress(iterations; dt=1.0, showspeed=true)

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    for iteration ∈ 1:iterations
        time = iteration

        update_phase_thermodynamics!(phase_eos[1], Val(1), nueff, ConstantScalar(300.0), model, config)
        update_phase_thermodynamics!(phase_eos[2], Val(2), nueff, ConstantScalar(300.0), model, config)

        blend_properties!(rho, alpha, phases[1].rho, phases[2].rho) # mixture.correct(); in OpenFOAM
        blend_viscosity!(alpha, phases, nu) # mixture.correct() in OpenFOAM

        interpolate!(rhof, rho, config)
        interpolate!(nuf, nu, config)

        grad!(∇rho, rhof, rho, time, config)
        limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

        phi_g!(phi_g, gh, ∇rho, config)

        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config; time=time) # U_eqn.solve(); in OpenFOAM
          
        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config) #rAU = 1.0/UEqn.A(); in OpenFOAM
        interpolate!(rDf, rD, config)
        remove_pressure_source!(U_eqn, ∇p_rgh, config)
        
        rp = 0.0
        for i ∈ 1:inner_loops
            H!(Hv, U, U_eqn, config) #HbyA in OpenFOAM
            interpolate!(Hvf, Hv, config)
            
            # Interpolate faces
            interpolate!(Uf, Hv, config) # Careful: reusing Uf for interpolation
            correct_boundaries!(Uf, Hv, boundaries.U, time, config)

            flux!(phiHbyA, Hvf, config)

            phi_gf!(phi_gf, rho, ghf, rDf, model, config)            
            flux!(phi_gf, config) #multiplies by area..... because if we multiply by normal then we get a vector...
            
            # flux!(mdotf, Uf, config)
            # flux!(phi_gf, config) #multiplies by area..... because if we multiply by normal then we get a vector...
            
            @. phiHbyA.values += phi_gf.values

            # div!(divHv, mdotf, config) #requires + phi_gf # [COMMENT OUT]
            div!(divHv, phiHbyA, config)
            
            @. prev = p_rgh.values
            rp = solve_equation!(p_eqn, p_rgh, boundaries.p_rgh, solvers.p_rgh, config; ref=pref, time=time) # p_rghEqn.solve();
            explicit_relaxation!(p_rgh, prev, 1.0, config)

            grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config) 
            limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)

            # correct_mass_flux(mdotf, p_rgh, rDf, config) # phi = phiHbyA - p_rghEqn.flux(); # [COMMENT OUT]
            # correct_mass_flux(phiHbyA, p_rgh, rDf, config)
            correct_mass_flux_multiphase(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, config)

            # correct_velocity!(U, Hv, ∇p_rgh, rD, config) # U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf); # [COMMENT OUT]
            correct_velocity_multiphase!(U, Hv, ∇p_rgh, rD, phi_g, config)

        end # corrector loop end
    
        @. p.values = p_rgh.values + (rho.values * gh.values)

        interpolate!(rhof_l, rho_l, config)
        # turbulence!(turbulenceModel, model, S, prev, time, config) # GIVES ERROR!!!
        update_nueff!(nueff, nuf, model.turbulence, config)

        # @. mdotf_VOF.values = mdotf.values #/ rhof.values [We actually never used rho in flux....]
        
        ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config; time=time)
        interpolate!(alphaf, alpha, config)

        maxCourant = max_courant_number!(cellsCourant, model, config)
        continuity_err!(rhof, Uf, model)

        R_ux[iteration] = rx
        R_uy[iteration] = ry
        R_uz[iteration] = rz
        R_p[iteration] = rp
        R_alpha[iteration] = ralpha

        Uz_convergence = true
        if typeof(mesh) <: Mesh3
            Uz_convergence = rz <= solvers.U.convergence
        end

        ProgressMeter.next!(
            progress, showvalues = [
                (:time, iteration*runtime.dt),
                (:Courant, maxCourant),
                (:Ux, R_ux[iteration]),
                (:Uy, R_uy[iteration]),
                (:Uz, R_uz[iteration]),
                (:p, R_p[iteration]),
                (:alpha, R_alpha[iteration])
                # turbulenceModel.state.residuals...
                ]
            )

        if iteration%write_interval + signbit(write_interval) == 0      
            save_output(model, outputWriter, iteration, time, config)
        end

    end # end for loop
    
    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p)
end



function correct_mass_flux_multiphase(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, config)
    # sngrad = FaceScalarField(mesh)
    (; faces, cells, boundary_cellsID) = rhof.mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_ifaces # length(n_ifaces) was a BUG! should be n_ifaces only!!!!
    kernel! = _correct_mass_flux_multiphase(_setup(backend, workgroup, ndrange)...)
    kernel!(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, faces, cells, n_bfaces)
    # KernelAbstractions.synchronize(backend)
end

@kernel function _correct_mass_flux_multiphase(phi, phiHbyA, rhoPhi, rhof, p_rgh, rDf, faces, cells, n_bfaces)
    i = @index(Global)
    fID = i + n_bfaces

    @inbounds begin 
        face = faces[fID]
        (; area, normal, ownerCells, delta) = face 
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        p1 = p_rgh[cID1]
        p2 = p_rgh[cID2]
        face_grad = area*(p2 - p1)/delta


        phi[fID] = phiHbyA[fID] - face_grad*rDf[fID]
        rhoPhi[fID] = phi[fID] * rhof[fID]
    end
end



function correct_velocity_multiphase!(U, Hv, ∇p_rgh, rD, phi_g, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(U)
    kernel! = _correct_velocity_multiphase!(_setup(backend, workgroup, ndrange)...)
    kernel!(U, Hv, ∇p_rgh, rD, phi_g)
end

@kernel function _correct_velocity_multiphase!(U, Hv, ∇p_rgh, rD, phi_g)
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

        Ux[i] = Hvx[i] + (phi_g_x[i] - dpdx[i]) * rD_i
        Uy[i] = Hvy[i] + (phi_g_y[i] - dpdy[i]) * rD_i
        Uz[i] = Hvz[i] + (phi_g_z[i] - dpdz[i]) * rD_i
    end
end





function continuity_err!(rhof, Uf, model)
    mesh = model.domain
    #total error ?
    # define local error = dt * (div of flux)
    # compute total error from locals

    for (cell_id, cell) in enumerate(mesh.cells)
        faces_range = cell.faces_range #???? mesh.cells[cell].faces_range

        mdotf = 0.0

        for faceID in faces_range
            face = mesh.faces[faceID]
            
            rho_local = rhof[faceID]
            U_local = Uf[faceID]

            
            (; area, normal) = face

            mdotf += dot(U_local, rho_local * area * normal)
        end

        if mdotf > 1.0e-3
            @warn "[CONTINUITY ERROR] at cell_ID=$cell_id >> mdtof_sum=$mdotf"
        end
    end
end




function blend_viscosity!(alpha, phases, nu)
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

function update_phase_thermodynamics!(EoS::AbstractEosModel, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    return nothing
end

function update_phase_thermodynamics!(EoS::Union{ConstEos, PerfectGas}, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    phase = model.fluid.phases[N]
    phase.eosModel(phase, model, config)
    phase.viscosityModel(phase, T)
end