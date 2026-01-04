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
    isInit = true

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

    ∇alphaf = FaceVectorField(mesh)
    interpolate!(∇alphaf, ∇alpha.result, config)

    mdotf = FaceScalarField(mesh)
    rhoPhi = FaceScalarField(mesh)    
    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)
    nueff = FaceScalarField(mesh)
    divHv = ScalarField(mesh)
    ddtCorr = VectorField(mesh)

    divCompressionFlux = ScalarField(mesh)
    compressionFlux = FaceScalarField(mesh)
    
    phi_g = VectorField(mesh)
    CSF = VectorField(mesh)
    phi_gf = FaceScalarField(mesh)
    phi_gf_flux = FaceScalarField(mesh)
    nf = FaceVectorField(mesh)
    kp = ScalarField(mesh)
    
    Srho = VectorField(mesh)
    
    interpolate!(alphaf, alpha, config)

    compute_compression!(compressionFlux, mdotf, alphaf, ∇alphaf, config)
    div!(divCompressionFlux, compressionFlux, config)

    compute_nf!(nf, ∇alphaf, config)
    div!(kp, nf, config)
    @. kp.values = -kp.values
    sigma = 0.03
    compute_CSF!(CSF, kp, sigma, ∇alpha, config)

    # Need to be defined before energyModel
    p_eqn = (
        - Laplacian{schemes.p.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)


    @info "Computing Fluid Properties..."


    phase_eos = [phases[1].eosModel, phases[2].eosModel]
    T_field = model.energy.T

    update_phase_thermodynamics!(phase_eos[1], Val(1), nueff, T_field, model, config)
    update_phase_thermodynamics!(phase_eos[2], Val(2), nueff, T_field, model, config)

    blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)
    blend_properties!(nu, alpha, phases[1].nu, phases[2].nu)

    interpolate_upwind!(rhof, rho, mdotf, config)
    interpolate_upwind!(nuf, nu, mdotf, config)

    gh = model.fluid.physics_properties.gravity.gh
    ghf = model.fluid.physics_properties.gravity.ghf
    g = model.fluid.physics_properties.gravity.g


    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)
    # compute_p_rgh!(p_rgh, gh, p, rho, config)
    # compute_p_rghf!(p_rghf, ghf, pf, rhof, config)

    ∇rho = Grad{schemes.p_rgh.gradient}(rho)
    grad!(∇rho, rhof, rho, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(rhoPhi, U)
        # + Divergence{schemes.U.divergence}(mdotf, U)
        - Laplacian{schemes.U.laplacian}(nueff, U) 
        ==
        - Source(∇p_rgh.result)
        + Source(Srho)
        # + Source(phi_g)
        # + Source(CSF)
    ) → VectorEquation(U, boundaries.U)

    alpha_eqn = (
        Time{schemes.alpha.time}(alpha)
        + Divergence{schemes.alpha.divergence}(mdotf, alpha)
        # + Divergence{schemes.alpha.divergence}(rhoPhi, alpha)
        == 
        - Source(divCompressionFlux)
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
        model, turbulenceModel, ∇p, ∇p_rgh, ∇rho, ∇alpha, ∇alphaf, Srho, U_eqn, p_eqn, alpha_eqn, mdotf, rhoPhi, gh, ghf, phi_g, phi_gf, phi_gf_flux, compressionFlux, divCompressionFlux, CSF, kp, sigma, nf, config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops)

    return residuals
end # end function



function MULTIPHASE(
    model, turbulenceModel, ∇p, ∇p_rgh, ∇rho, ∇alpha, ∇alphaf, Srho, U_eqn, p_eqn, alpha_eqn, mdotf, rhoPhi, gh, ghf, phi_g, phi_gf, phi_gf_flux, compressionFlux, divCompressionFlux, CSF, kp, sigma, nf, config;

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
    # mdotf = get_flux(U_eqn, 2)
    nueff = get_flux(U_eqn, 3)
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)

    @info "Allocating working memory..."

    # Define aux fields 
    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    Uf = FaceVectorField(mesh)
    S = StrainRate(gradU, gradUT, U, Uf)
    # p_grad_field = ScalarField(mesh)
    ∇p_rghf_deconstructed = FaceScalarField(mesh)
    ∇p_rghf_reconstructed = VectorField(mesh)
    # reconstructed_U = VectorField(mesh)

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
    
    # Initial calculations
    time = zero(TF) # assuming time=0
    interpolate!(Uf, U, config)   
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    
    flux!(mdotf, Uf, config)
    @. rhoPhi.values = mdotf.values * rhof.values

    phase_eos = [phases[1].eosModel, phases[2].eosModel]
    T_field = model.energy.T

    update_nueff!(nueff, nuf, model.turbulence, config)

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting multiphase loops..."

    progress = Progress(iterations; dt=1.0, showspeed=true)

    rho_prev = ScalarField(mesh)
    impulse = ScalarField(mesh)

    @time for iteration ∈ 1:iterations
        time = iteration *dt

        grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)
        interpolate!(∇alphaf, ∇alpha.result, config)

        compute_compression!(compressionFlux, mdotf, alphaf, ∇alphaf, config)
        div!(divCompressionFlux, compressionFlux, config)

        correct_boundaries!(alphaf, alpha, boundaries.alpha, time, config)
        ralpha = solve_equation!(alpha_eqn, alpha, boundaries.alpha, solvers.alpha, config; time=time)
        interpolate!(alphaf, alpha, config)
        correct_boundaries!(alphaf, alpha, boundaries.alpha, time, config)

        update_phase_thermodynamics!(phase_eos[1], Val(1), nueff, T_field, model, config)
        update_phase_thermodynamics!(phase_eos[2], Val(2), nueff, T_field, model, config)

        @. rho_prev.values = rho.values
        blend_properties!(rho, alpha, phases[1].rho, phases[2].rho)
        blend_properties!(nu, alpha, phases[1].nu, phases[2].nu)

        interpolate_upwind!(rhof, rho, mdotf, config)
        interpolate_upwind!(nuf, nu, mdotf, config)

        grad!(∇rho, rhof, rho, time, config)
        limit_gradient!(schemes.p_rgh.limiter, ∇rho, rho, config)

        @. rhoPhi.values = mdotf.values * rhof.values

        compute_nf!(nf, ∇alphaf, config)
        div!(kp, nf, config)
        @. kp.values = -kp.values
        compute_CSF!(CSF, kp, sigma, ∇alpha, config)

        compute_Srho!(Srho, rho, rho_prev, U, config, dt)

        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config; time=time)

        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config)
        interpolate!(rDf, rD, config)

        remove_pressure_source!(U_eqn, ∇p_rgh, config)
        
        rp = 0.0
        # println("<<< NEW ITERATION >>>")
        for i ∈ 1:inner_loops
            # println("\n\npiso\n")
            H!(Hv, U, U_eqn, config)
            
            interpolate!(Uf, Hv, config)
            correct_boundaries!(Uf, Hv, boundaries.U, time, config)

            flux!(mdotf, Uf, config)
            # @. rhoPhi.values = mdotf.values * rhof.values

            phi_gf!(phi_gf, rho, ghf, rDf, model, config)
            reconstruct_operation!(phi_g, phi_gf, config)

            @. mdotf.values += phi_gf.values
            # @. rhoPhi.values = mdotf.values * rhof.values

            correct_boundaries!(p_rghf, p_rgh, boundaries.p_rgh, time, config)

            div!(divHv, mdotf, config)
            
            @. prev = p_rgh.values
            rp = solve_equation!(p_eqn, p_rgh, boundaries.p_rgh, solvers.p_rgh, config; ref=0.0, time=time)

            if i == inner_loops
                explicit_relaxation!(p_rgh, prev, 1.0, config)
            else
                explicit_relaxation!(p_rgh, prev, solvers.p_rgh.relax, config)
            end
            correct_boundaries!(p_rghf, p_rgh, boundaries.p_rgh, time, config)

            grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config) 
            limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)
            
            correct_mass_flux(mdotf, p_rgh, rDf, config)
            # @. rhoPhi.values = mdotf.values * rhof.values
            
            pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)
            reconstruct_operation!(∇p_rghf_reconstructed, ∇p_rghf_deconstructed, config)
            correct_velocity_rgh!(U, Hv, ∇p_rghf_reconstructed, rD, phi_g, config)
        end # corrector loop end
    interpolate!(Uf, U, config)
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    @. rhoPhi.values = mdotf.values * rhof.values


    # divv = ScalarField(mesh)
    # div!(divv, Uf, config)
    # imax = argmax(divv.values)
    # println("Max U divergence = $(divv.values[imax]) at index $imax; and MEAN FIELD: $(mean(divv.values))")
    
        
    @. p.values = p_rgh.values + (rho.values * gh.values)
    interpolate_upwind!(pf, p, mdotf, config)

    turbulence!(turbulenceModel, model, S, prev, time, config)
    update_nueff!(nueff, nuf, model.turbulence, config)

    # if iteration % 500 == 0
    #     # println("not volume weighted but mean of Ux: $(mean(U.x.values))")
    #     println("U_x max: $(maximum(U.x.values))")
    #     println("U_y max: $(maximum(U.y.values))")
    # end

    maxCourant = max_courant_number!(cellsCourant, model, config)

    R_ux[iteration] = rx
    R_uy[iteration] = ry
    R_uz[iteration] = rz
    R_p[iteration] = rp
    # R_alpha[iteration] = ralpha

    ProgressMeter.next!(
        progress, showvalues = [
            (:time, iteration*runtime.dt),
            (:Courant, maxCourant),
            (:Ux, R_ux[iteration]),
            (:Uy, R_uy[iteration]),
            (:Uz, R_uz[iteration]),
            (:p_rgh, R_p[iteration]),
            # (:alpha, R_alpha[iteration]),
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


function compute_compression!(compressionFlux, mdotf, alphaf, ∇alphaf, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    faces = mdotf.mesh.faces

    ndrange = length(mdotf)
    kernel! = _compute_compression!(_setup(backend, workgroup, ndrange)...)
    kernel!(compressionFlux, mdotf, alphaf, ∇alphaf, faces)
end
@kernel inbounds=true function _compute_compression!(compressionFlux, mdotf, alphaf, ∇alphaf, faces)
    i = @index(Global)

    (; area, normal) = faces[i]

    # C_alpha = 1.0
    C_alpha = 0.0

    # try abs(mdotf)/abs(Sf)

    compressionFlux[i] = C_alpha * abs(mdotf[i]) * (alphaf[i] * (1.0 - alphaf[i])) * dot((∇alphaf[i] / (norm(∇alphaf[i]) + eps())), (normal)) #non GPU #area
end



function compute_Srho!(Srho, rho, rho_prev, U, config, dt)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(rho)
    kernel! = _compute_Srho!(_setup(backend, workgroup, ndrange)...)
    kernel!(Srho, rho, rho_prev, U, config, dt)
end

@kernel function _compute_Srho!(Srho, rho, rho_prev, U, config, dt)
    i = @index(Global)

    @uniform begin
        Ux, Uy, Uz = U.x, U.y, U.z
        Srhox, Srhoy, Srhoz = Srho.x, Srho.y, Srho.z
    end

    @inbounds begin
        drho_dt = (rho[i] - rho_prev[i]) / dt # + eps()

        Srhox[i] = -Ux[i] * drho_dt
        Srhoy[i] = -Uy[i] * drho_dt
        Srhoz[i] = -Uz[i] * drho_dt
    end
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
        # dpdx, dpdy, dpdz = ∇p.result.x, ∇p.result.y, ∇p.result.z
        phi_gx, phi_gy, phi_gz = phi_g.x, phi_g.y, phi_g.z
        rDvalues = rD.values
    end

    @inbounds begin
        rDvalues_i = rDvalues[i]

        if i == 940
        #     # println("Ux: $(Ux[i])")
            # println("dpdy 940: $(dpdy[i])")
        #     println("phigy 940: $(phi_gy[i])")
        #     # println("Hvy: $(Hvy[i])")
        #     # diff = abs(-dpdy[i]*rDvalues_i - phi_gy[i]*rDvalues_i)
        #     # println("Diff: $diff %")
        end
        if i == 941
            # println("dpdy 941: $(dpdy[i])")
        #     println("phigy 941: $(phi_gy[i])")
            
        #     # diff = abs(-dpdy[i]*rDvalues_i - phi_gy[i]*rDvalues_i)
        #     # println("Diff: $diff %")
        #     # println("dpdy of 941: $(dpdy[i]*rDvalues_i)")
        #     # println("phi_gy: $(phi_gy[i]*rDvalues_i)")
        end


        Ux[i] = Hvx[i] + dpdx[i] * rDvalues_i #+ phi_gx[i] * rDvalues_i
        Uy[i] = Hvy[i] + dpdy[i] * rDvalues_i #- phi_gy[i] * rDvalues_i
        Uz[i] = Hvz[i] + dpdz[i] * rDvalues_i #- phi_gz[i] * rDvalues_i
    end
end


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
            # nsign is optional IF psif is consistent with faces[fID].normal orientation
            # keep it available if you later discover a sign convention mismatch
            # s = cell_nsign[fi]

            (; area, normal) = faces[fID]
            nx = normal[1]; ny = normal[2]

            # M += area * (n ⊗ n)
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



# function reconstruct_operation!(phi, psif, config) # phi = VectorField, psif = FaceScalarField
#     (; hardware) = config
#     (; backend, workgroup) = hardware

#     cells = phi.mesh.cells
#     faces = phi.mesh.faces

#     ndrange = length(phi_gf)
#     kernel! = _reconstruct_operation!(_setup(backend, workgroup, ndrange)...)
#     kernel!(phi, psif, faces)
# end
# @kernel function _reconstruct_operation!(phi, psif, faces)
#     i = @index(Global)
#     face = faces[i]
#     (; area, normal, ownerCells, delta) = face

#     cID1 = ownerCells[1]
#     cID2 = ownerCells[2]
#     p1 = p[cID1]
#     p2 = p[cID2]
#     face_grad = area*(p2 - p1)/delta # best option so far!

#     deconstructed_U[i] = phi_gf[i]
# end


function deconstructed_U_compute!(deconstructed_U, phi_gf, p, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    faces = p.mesh.faces

    ndrange = length(phi_gf)
    kernel! = _deconstructed_U_compute!(_setup(backend, workgroup, ndrange)...)
    kernel!(deconstructed_U, phi_gf, p, faces)
end
@kernel function _deconstructed_U_compute!(deconstructed_U, phi_gf, p, faces)
    i = @index(Global)
    face = faces[i]
    (; area, normal, ownerCells, delta) = face

    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    p1 = p[cID1]
    p2 = p[cID2]
    face_grad = area*(p2 - p1)/delta # best option so far!

    deconstructed_U[i] = phi_gf[i]
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

function compute_nf!(nf, ∇alphaf, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(nf.x)
    kernel! = _compute_nf!(_setup(backend, workgroup, ndrange)...)
    kernel!(nf, ∇alphaf)
end
@kernel inbounds=true function _compute_nf!(nf, ∇alphaf)
    i = @index(Global)

    nf[i] = ∇alphaf[i] / (norm(∇alphaf[i]) + eps())
end


function compute_CSF!(CSF, kp, sigma, ∇alpha, config)
    (; hardware) = config
    backend = hardware.backend
    workgroup = hardware.workgroup

    ndrange = length(kp)
    kernel! = _compute_CSF!(_setup(backend, workgroup, ndrange)...)
    kernel!(CSF, kp, sigma, ∇alpha)
end
@kernel inbounds=true function _compute_CSF!(CSF, kp, sigma, ∇alpha)
    i = @index(Global)

    CSF[i] = ∇alpha[i] * sigma * kp[i]
end


function update_phase_thermodynamics!(EoS::AbstractEosModel, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    return nothing
end

function update_phase_thermodynamics!(EoS::Union{ConstEos, PerfectGas}, phaseIndex::Val{N}, nueff, T, model, config) where {N}
    phase = model.fluid.phases[N]
    phase.eosModel(phase, model, config)
    phase.viscosityModel(phase, model)
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
    # gh[i] = (g ⋅ (centre-[0.0,0.5,0.0]))
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
    # ghf[i] = (g ⋅ (centre-[0.0,0.5,0.0]))
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

    # println("Number of faces: $n_faces")
    # ndrange = n_ifaces
    ndrange = n_faces #try all faces!
    kernel! = _phi_gf!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_gf, rho, ghf, rDf, faces, cells, n_bfaces)
end

@kernel function _phi_gf!(phi_gf, rho, ghf, rDf, faces, cells, n_bfaces)
    # i = @index(Global)
    # fID = i + n_bfaces
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

        # if face_grad > 1e-5
        #     println("face_grad: $face_grad; ghf: $(ghf[fID]); rDf=$(rDf[fID]); TOTAL=$(-ghf[fID] * face_grad * rDf[fID])")
        # end
    end
end



function reconstruct!(phi::VectorField, psif::FaceScalarField, config)
    mesh = phi.mesh
    (; cells, cell_nsign, cell_faces, faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    F = _get_float(mesh)

    # Launch main calculation kernel
    ndrange = length(cells)
    kernel! = _reconstruct_internal!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, F, cell_faces, cell_nsign, faces, phi, psif)

    # Retrieve number of boundary faces
    nbfaces = length(mesh.boundary_cellsID)

    # Launch boundary faces contribution kernel
    ndrange = nbfaces
    kernel! = _reconstruct_boundaries!(_setup(backend, workgroup, ndrange)...)
    kernel!(faces, cells, phi, psif)
end



@kernel function _reconstruct_internal!(cells::AbstractArray{Cell{TF,SV,UR}}, F, cell_faces, cell_nsign, faces, phi, psif) where {TF,SV,UR}
    i = @index(Global)
    
    @inbounds begin
        (; volume, faces_range) = cells[i]
        
        reduction = zero(SV)
        
        for fi ∈ faces_range
            fID = cell_faces[fi]
            nsign = cell_nsign[fi]

            (; area, normal) = faces[fID]

            reduction += psif[fID]*nsign*normal*area
        end
        phi[i] = reduction#/volume
    end
end

@kernel function _reconstruct_boundaries!(faces, cells, phi, psif)
    i = @index(Global)
    
    @inbounds begin
        cID = faces[i].ownerCells[1]
        volume = cells[cID].volume
        (; area, normal) = faces[i]

        phi[cID] = phi[cID] + normal * psif[i]*area #/ volume
    end
end