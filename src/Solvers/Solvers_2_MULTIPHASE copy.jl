export multiphase!

"""
    multiphase!(model, config;
        output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2)

Multiphase (VOF) solver baseline — PISO structure with gravity injected as a
direct face flux and reconstructed to the cell-centre gradient used in velocity
correction. Phase transport and surface tension will be layered on top of this.

# Input arguments

- `model` reference to a `Physics` model defined by the user.
- `config` Configuration structure defined by the user with solvers, schemes, runtime and hardware structures configuration details.
- `output` select the format used for simulation results from `VTK()` or `OpenFOAM` (default = `VTK()`)
- `pref` Reference pressure value for cases that do not have a pressure defining BC (default = `nothing`)
- `ncorrectors` number of non-orthogonality correction loops (default = `0`)
- `inner_loops` number of PISO inner loops (default = `2`)
"""
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
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    )

    (; solvers, schemes, runtime, hardware, boundaries) = config

    @info "Extracting configuration and input fields..."

    (; U, p) = model.momentum
    (; alpha, alphaf, rho, rhof, nu, nuf, p_rgh, p_rghf) = model.fluid

    phases = model.fluid.phases
    volume_fraction = model.fluid.volume_fraction
    main = volume_fraction
    secondary = 3 - volume_fraction

    mesh = model.domain

    @info "Pre-allocating fields..."

    TF = _get_float(mesh)
    time = zero(TF)

    ∇p = Grad{schemes.p_rgh.gradient}(p)

    ∇p_rgh = Grad{schemes.p_rgh.gradient}(p_rgh)
    grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
    limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)

    mdotf = FaceScalarField(mesh)
    rhoPhi = FaceScalarField(mesh)
    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)
    nueff = FaceScalarField(mesh)
    divHv = ScalarField(mesh)

    phi_g = VectorField(mesh)       # reconstructed gravity gradient (cell-centred)
    phi_gf = FaceScalarField(mesh)  # gravity face flux injection

    @info "Computing fluid properties..."

    blend_properties!(rho, alpha, phases[main].rho[1], phases[secondary].rho[1])
    blend_properties!(rhof, alphaf, phases[main].rho[1], phases[secondary].rho[1])
    blend_properties!(nuf, alphaf, phases[main].mu[1] / phases[main].rho[1], phases[secondary].mu[1] / phases[secondary].rho[1])

    gh = model.fluid.physics_properties.gravity.gh
    ghf = model.fluid.physics_properties.gravity.ghf
    g = model.fluid.physics_properties.gravity.g

    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)

    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(rhoPhi, U)
        - Laplacian{schemes.U.laplacian}(nueff, U)
        ==
        - Source(∇p_rgh.result)
    ) → VectorEquation(U, boundaries.U)

    p_eqn = (
        - Laplacian{schemes.p.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)

    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p_rgh.preconditioner, p_eqn)

    @info "Pre-allocating solvers..."

    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p_rgh.solver, _b(p_eqn))

    @info "Initialising turbulence model..."
    turbulenceModel, config = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals = solver_variant(
        model, turbulenceModel, ∇p, ∇p_rgh, U_eqn, p_eqn,
        mdotf, rhoPhi, gh, ghf, phi_g, phi_gf, config;
        output=output, pref=pref,
        ncorrectors=ncorrectors, inner_loops=inner_loops)

    return residuals
end

function MULTIPHASE(
    model, turbulenceModel, ∇p, ∇p_rgh, U_eqn, p_eqn,
    mdotf, rhoPhi, gh, ghf, phi_g, phi_gf, config;
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=2
    )

    (; U, p) = model.momentum
    (; nu, nuf, rho, rhof, alpha, alphaf, p_rgh, p_rghf) = model.fluid
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries, postprocess) = config
    (; iterations, write_interval) = runtime
    (; backend) = hardware

    dt_cpu = zeros(_get_float(mesh), 1)
    copyto!(dt_cpu, config.runtime.dt)

    postprocess = convert_time_to_iterations(postprocess, model, dt_cpu[1], iterations)

    nueff = FaceScalarField(mesh)
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)

    @info "Allocating working memory..."

    # Define aux fields
    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    Uf = FaceVectorField(mesh)
    S = StrainRate(gradU, gradUT, U, Uf)

    # Working fields for pressure-gradient decomposition and reconstruction
    ∇p_rghf_deconstructed = FaceScalarField(mesh)
    ∇p_rghf_reconstructed = VectorField(mesh)

    # Alpha transport working fields (MULES-bounded, interFOAM-style compression)
    alpha_prev      = ScalarField(mesh)
    div_alpha       = ScalarField(mesh)
    div_mdotf       = ScalarField(mesh)       # ∇·mdotf for div-corrected α transport
    alpha_fluxf     = FaceScalarField(mesh)   # F_final in interFOAM — limited alpha face flux
    alphaf_upwind   = FaceScalarField(mesh)
    alphaf_HO       = FaceScalarField(mesh)
    phirf           = FaceScalarField(mesh)   # compression face flux (interFOAM phir*Sf)
    phiLf           = FaceScalarField(mesh)   # bounded low-order flux
    phiHf           = FaceScalarField(mesh)   # high-order flux (advection + compression)
    phiAf           = FaceScalarField(mesh)   # anti-diffusive correction (HO − LO), limited in place
    Pplus           = ScalarField(mesh)
    Pminus          = ScalarField(mesh)
    Qplus           = ScalarField(mesh)
    Qminus          = ScalarField(mesh)
    Rplus           = ScalarField(mesh)
    Rminus          = ScalarField(mesh)
    C_alpha         = 1.0                    # interFOAM cAlpha (sharpening strength)

    # Surface tension (CSF) working fields — σ κ ∇α injected into phi_gf
    sigma = 2.0
    # sigma = 0.0025 # RTI for the new RTI or 0.0025

    # Surface tension (CSF) working fields — σ κ ∇α injected into phi_gf
    # .sigma = 2.0
    ∇alpha          = Grad{schemes.alpha.gradient}(alpha)
    ∇alphaf         = FaceVectorField(mesh)
    nhatf_prep      = FaceVectorField(mesh)
    kappa           = ScalarField(mesh)
    kappaf          = FaceScalarField(mesh)
    ∇kappa          = Grad{schemes.alpha.gradient}(kappa)
    alpha_smooth    = ScalarField(mesh)       # Laplacian-smoothed α for curvature
    alpha_smooth_tmp = ScalarField(mesh)      # ping-pong buffer for smoothing passes
    n_smooth        = 0                       # number of Laplacian smoothing passes

    phases = model.fluid.phases
    volume_fraction = model.fluid.volume_fraction
    main = volume_fraction
    secondary = 3 - volume_fraction

    n_cells = length(mesh.cells)
    Hv = VectorField(mesh)
    rD = ScalarField(mesh)
    rho_prev = ScalarField(mesh)

    TF = _get_float(mesh)
    TI = _get_int(mesh)
    prev = KernelAbstractions.zeros(backend, TF, n_cells)

    R_ux = ones(TF, iterations)
    R_uy = ones(TF, iterations)
    R_uz = ones(TF, iterations)
    R_p  = ones(TF, iterations)
    cellsCourant = KernelAbstractions.zeros(backend, TF, n_cells)
    cellsAlphaCourant = KernelAbstractions.zeros(backend, TF, n_cells)

    time = zero(TF)
    interpolate!(Uf, U, config)
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    flux!(mdotf, Uf, config)

    # Initialise rhoPhi before first U solve — phase state assumed static in this baseline
    @. rhoPhi.values = mdotf.values * rhof.values

    update_nueff!(nueff, nuf, model.turbulence, config)

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting multiphase PISO loops..."

    progress = Progress(iterations; dt=1.0, showspeed=true)

    @time for iteration ∈ 1:iterations
        copyto!(dt_cpu, config.runtime.dt)
        time += dt_cpu[1]

        @. rho_prev.values = rho.values

        # α transport with interFOAM-style compression + MULES bounded update.
        # HO advection: van Leer interpolation of α onto faces using current ∇α.
        # Compression:  phirf = min(Cα|φ|/|Sf|, max|φ|/|Sf|) · (n̂f · Sf), with n̂f from ∇αf.
        # Anti-diffusive flux Δf = φ_HO − φ_LO is Zalesak-limited per face so no cell
        # leaves [0,1]; final α = α_LO − dt/V · div(limited Δf).
        @. alpha_prev.values = alpha.values

        grad!(∇alpha, alphaf, alpha, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)
        interpolate!(∇alphaf, ∇alpha.result, config)

        interpolate_upwind!(alphaf_upwind, alpha, mdotf, config)
        interpolate_vanleer!(alphaf_HO, alpha, ∇alpha, mdotf, config)

        compression_flux!(phirf, ∇alphaf, mdotf, C_alpha, config)

        @. phiLf.values = mdotf.values * alphaf_upwind.values
        @. phiHf.values = mdotf.values * alphaf_HO.values +
                          phirf.values * alphaf_upwind.values * (1.0 - alphaf_upwind.values)
        @. phiAf.values = phiHf.values - phiLf.values

        # Zalesak cannot limit boundary faces (only one owner), so drop the
        # anti-diffusive correction there and let the upwind LO flux pass through.
        zero_boundary_faces!(phiAf, config)

        mules_limit!(phiAf, alpha_prev, phiLf,
                     Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
                     dt_cpu[1], mesh, config)

        @. alpha_fluxf.values = phiLf.values + phiAf.values
        div!(div_alpha, alpha_fluxf, config)
        div!(div_mdotf, mdotf, config)
        @. alpha.values = alpha_prev.values -
            dt_cpu[1] * (div_alpha.values - alpha_prev.values * div_mdotf.values)

        # Refresh cell-centred and face properties from the new α.
        # interpolate_upwind!(alphaf, alpha, mdotf, config)
        interpolate_vanleer!(alphaf, alpha, ∇alpha, mdotf, config) # the shape is so much better with this

        blend_properties!(rho,  alpha,  phases[main].rho[1], phases[secondary].rho[1])
        blend_properties!(rhof, alphaf, phases[main].rho[1], phases[secondary].rho[1])
        blend_properties!(nuf,  alphaf, phases[main].mu[1] / phases[main].rho[1],
                                        phases[secondary].mu[1] / phases[secondary].rho[1])

        # interFOAM rhoPhi: F_final·(ρ1 − ρ2) + mdotf·ρ2
        # alpha_fluxf is F_final (the MULES-limited alpha face flux).
        # This correctly weights the density flux by the phase composition at each face.
        rho1 = phases[main].rho[1]
        rho2 = phases[secondary].rho[1]
        @. rhoPhi.values = alpha_fluxf.values * (rho1 - rho2) + mdotf.values * rho2

        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config, rho_prev; time=time)

        # Pressure correction
        inverse_diagonal!(rD, U_eqn, config)
        interpolate!(rDf, rD, config)

        remove_pressure_source!(U_eqn, ∇p_rgh, config)

        # Compute κ once per outer iteration from a Laplacian-smoothed α.
        # Smoothing reduces curvature noise that drives parasitic currents.
        laplacian_smooth!(alpha_smooth, alpha_smooth_tmp, alpha, n_smooth, config)
        grad!(∇alpha, alphaf, alpha_smooth, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha_smooth, config)
        interpolate!(∇alphaf, ∇alpha.result, config)
        nhat_prep!(nhatf_prep, alpha_smooth, ∇alphaf, config)
        div!(kappa, nhatf_prep, config)
        interpolate!(kappaf, kappa, config)

        rp = 0.0
        for i ∈ 1:inner_loops
            H!(Hv, U, U_eqn, config)

            interpolate!(Uf, Hv, config)
            correct_boundaries!(Uf, Hv, boundaries.U, time, config)

            flux!(mdotf, Uf, config)

            # Gravity + frozen CSF injected into face flux, then reconstructed.
            phi_gf!(phi_gf, rho, ghf, rDf, model, config)
            surface_tension_flux!(rDf, sigma, kappaf, alpha, ∇alphaf, phi_gf, config)

            reconstruct_operation!(phi_g, phi_gf, config)
            @. mdotf.values += phi_gf.values

            div!(divHv, mdotf, config)

            @. prev = p_rgh.values
            rp = solve_equation!(p_eqn, p_rgh, boundaries.p_rgh, solvers.p_rgh, config, rho_prev; ref=0.0, time=time)

            grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
            limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)

            correct_mass_flux(mdotf, p_eqn, config)

            interpolate_upwind!(alphaf, alpha, mdotf, config)
            blend_properties!(rhof, alphaf, phases[main].rho[1], phases[secondary].rho[1])
            @. rhoPhi.values = alpha_fluxf.values * (rho1 - rho2) + mdotf.values * rho2

            pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)
            reconstruct_operation!(∇p_rghf_reconstructed, ∇p_rghf_deconstructed, config)
            correct_velocity_rgh!(U, Hv, ∇p_rghf_reconstructed, rD, phi_g, config)
        end # corrector loop end

        @. p.values = p_rgh.values + (rho.values * gh.values)

        turbulence!(turbulenceModel, model, S, prev, time, config)
        update_nueff!(nueff, nuf, model.turbulence, config)

        courant = max_courant_number!(cellsCourant, model, config)
        alphaCourant = max_alpha_courant_number!(cellsAlphaCourant, alpha, mdotf, model, config, dt_cpu[1])
        update_dt!(config.runtime, courant, alphaCourant)

        R_ux[iteration] = rx
        R_uy[iteration] = ry
        R_uz[iteration] = rz
        R_p[iteration]  = rp

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
                turbulenceModel.state.residuals...
                ]
            )

        runtime_postprocessing!(postprocess, iteration, iterations)

        if iteration % write_interval + signbit(write_interval) == 0
            save_output(model, outputWriter, iteration, time, config)
            save_postprocessing(postprocess, iteration, time, mesh, outputWriter, config.boundaries)
        end
    end # end for loop

    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p)
end


## BASELINE HELPERS (gravity flux injection, reconstruction, blending)

function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0 * alpha_field.values) + (property_1 * (1.0 - alpha_field.values))
    nothing
end

"""
    compute_gh!(gh, g, config)

Computes `g . x` at cell centres. Used for hydrostatic pressure reconstruction.
"""
function compute_gh!(gh, g, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    cells = gh.mesh.cells

    ndrange = length(gh)
    kernel! = _compute_gh!(_setup(backend, workgroup, ndrange)...)
    kernel!(gh, g, cells)
end
@kernel inbounds=true function _compute_gh!(gh, g, cells)
    i = @index(Global)
    (; centre) = cells[i]
    gh[i] = (g ⋅ centre)
end

"""
    compute_ghf!(ghf, g, config)

Computes `g . x` at face centres.
"""
function compute_ghf!(ghf, g, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    faces = ghf.mesh.faces

    ndrange = length(ghf)
    kernel! = _compute_ghf!(_setup(backend, workgroup, ndrange)...)
    kernel!(ghf, g, faces)
end
@kernel inbounds=true function _compute_ghf!(ghf, g, faces)
    i = @index(Global)
    (; centre) = faces[i]
    ghf[i] = (g ⋅ centre)
end

"""
    phi_gf!(phi_gf, rho, ghf, rDf, model, config)

Gravity contribution to the face flux for pressure-velocity coupling.
Formula: `phi_gf = -ghf * (dρ/dn)_face * rDf`
"""
function phi_gf!(phi_gf, rho, ghf, rDf, model, config)
    (; faces, cells, boundary_cellsID) = model.domain
    (; hardware) = config
    (; backend, workgroup) = hardware

    n_bfaces = length(boundary_cellsID)

    ndrange = length(faces)
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

        face_grad = area * (rho2 - rho1) / delta

        phi_gf[fID] = -ghf[fID] * face_grad * rDf[fID]
    end
end

"""
    pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)

Per-face pressure gradient with the gravity contribution removed, ready for
reconstruction to a cell-centred vector field.
"""
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
    face_grad = area * (p2 - p1) / delta

    ∇p_rghf_deconstructed[i] = (phi_gf[i] - (face_grad * rDf[i])) / (rDf[i] + eps())
end

## Velocity correction (uses reconstructed ∇p_rgh that embeds gravity)

function correct_velocity_rgh!(U, Hv, ∇p, rD, phi_g, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(U)
    kernel! = _correct_velocity_rgh!(_setup(backend, workgroup, ndrange)...)
    kernel!(U, Hv, ∇p, rD, phi_g)
end
@kernel function _correct_velocity_rgh!(U, Hv, ∇p, rD, phi_g)
    i = @index(Global)

    @uniform begin
        Ux, Uy, Uz = U.x, U.y, U.z
        Hvx, Hvy, Hvz = Hv.x, Hv.y, Hv.z
        dpdx, dpdy, dpdz = ∇p.x, ∇p.y, ∇p.z
        rDvalues = rD.values
    end

    @inbounds begin
        rDvalues_i = rDvalues[i]
        Ux[i] = Hvx[i] + dpdx[i] * rDvalues_i
        Uy[i] = Hvy[i] + dpdy[i] * rDvalues_i
        Uz[i] = Hvz[i] + dpdz[i] * rDvalues_i
    end
end

## Compression flux (interFOAM-style phir_f · Sf)
##
## On each face:
##   phir_f · Sf = min(Cα · |φf|/|Sf|, Φmax) · (n̂f · Sf)
## with n̂f = (∇α)f / |(∇α)f|, and Φmax = max_f |φf|/|Sf| (global reference speed).
## The returned FaceScalarField already has the dot with Sf baked in, so it slots
## straight into the alpha flux as `phirf · α(1−α)`.

function compression_flux!(phirf, ∇alphaf, mdotf, C_alpha, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh  = phirf.mesh
    faces = mesh.faces
    TF = _get_float(mesh)

    # Global reference speed Φmax = max |φf|/|Sf|. Fill a scratch face buffer
    # with |φf|/|Sf| via a kernel, then reduce. This stays on-backend.
    phi_over_S_buf = similar(mdotf.values)
    ndrange = length(faces)
    kernel! = _fill_phi_over_S!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_over_S_buf, mdotf, faces)
    phimax = maximum(phi_over_S_buf)

    ndrange = length(faces)
    kernel! = _compression_flux!(_setup(backend, workgroup, ndrange)...)
    kernel!(phirf, ∇alphaf, mdotf, faces, TF(C_alpha), TF(phimax))
end

@kernel inbounds=true function _fill_phi_over_S!(buf, mdotf, faces)
    i = @index(Global)
    TF = eltype(buf)
    area = faces[i].area
    buf[i] = abs(mdotf[i]) / (area + eps(TF))
end

@kernel inbounds=true function _compression_flux!(phirf, ∇alphaf, mdotf, faces, C_alpha, phimax)
    i = @index(Global)
    face = faces[i]
    (; area, normal) = face

    Sf = area * normal
    g  = ∇alphaf[i]
    mag = norm(g)

    nhat = mag > eps(eltype(phirf.values)) ? g / mag : zero(g)

    phi_over_S = abs(mdotf[i]) / (area + eps(eltype(phirf.values)))
    compr_speed = min(C_alpha * phi_over_S, phimax)

    phirf[i] = compr_speed * (nhat ⋅ Sf)
end


## Boundary-face zeroing helper (used before MULES limiter)

function zero_boundary_faces!(phif::FaceScalarField, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = phif.mesh
    nbfaces = length(mesh.boundary_cellsID)
    if nbfaces > 0
        ndrange = nbfaces
        kernel! = _zero_boundary_faces!(_setup(backend, workgroup, ndrange)...)
        kernel!(phif)
    end
end
@kernel inbounds=true function _zero_boundary_faces!(phif)
    i = @index(Global)
    phif[i] = zero(eltype(phif.values))
end


## MULES limiter (Zalesak's flux-corrected transport)
##
## Given:
##   - α^n (alpha_prev)
##   - bounded low-order flux φLf (upwind)
##   - anti-diffusive flux Δf = φHf − φLf (passed in via phiAf, limited in place)
## Compute per-cell slack Q± = (bound − α*)·V/dt where α* is the low-order update,
## accumulate inflow/outflow of Δf into P±, form R± = min(1, Q±/P±), and for each
## internal face take λf = min(R+_down, R−_up) with the direction set by sign(Δf)
## relative to the owner. Boundary faces are left unlimited (λ=1): they don't
## threaten boundedness for scalar α under incompressible transport.

function mules_limit!(phiAf, alpha_prev, phiLf,
                      Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
                      dt, mesh, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    (; cells, cell_nsign, cell_faces, faces) = mesh

    n_cells  = length(cells)
    n_faces  = length(faces)
    nbfaces  = length(mesh.boundary_cellsID)

    # Reset accumulators
    fill!(Pplus.values,  0)
    fill!(Pminus.values, 0)

    # 1) Accumulate low-order α* implicitly by computing Q± from α_prev and the
    #    bounded flux, and P± from signed anti-diffusive contributions.
    ndrange = n_cells
    kernel! = _mules_cell_accum!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, cell_faces, cell_nsign, faces,
            alpha_prev, phiLf, phiAf,
            Pplus, Pminus, Qplus, Qminus, dt)

    # 2) Form R± = min(1, Q±/P±) per cell.
    ndrange = n_cells
    kernel! = _mules_ratios!(_setup(backend, workgroup, ndrange)...)
    kernel!(Pplus, Pminus, Qplus, Qminus, Rplus, Rminus)

    # 3) Apply per-face limiter: scale Δf in place (internal faces only).
    ndrange = n_faces
    kernel! = _mules_apply!(_setup(backend, workgroup, ndrange)...)
    kernel!(phiAf, faces, Rplus, Rminus, nbfaces)
end

@kernel inbounds=true function _mules_cell_accum!(
    cells::AbstractArray{Cell{TF,SV,UR}}, cell_faces, cell_nsign, faces,
    alpha_prev, phiLf, phiAf, Pplus, Pminus, Qplus, Qminus, dt
) where {TF,SV,UR}
    i = @index(Global)
    (; volume, faces_range) = cells[i]

    # α* from the low-order (bounded) update: α* = α_prev − dt/V · Σ nsign·φLf
    sum_L = zero(TF)
    sum_A_pos = zero(TF)  # Δf contribution that would INCREASE α in this cell
    sum_A_neg = zero(TF)  # Δf contribution that would DECREASE α in this cell

    for fi in faces_range
        fID   = cell_faces[fi]
        nsign = cell_nsign[fi]

        fL = phiLf[fID] * nsign
        fA = phiAf[fID] * nsign

        sum_L += fL
        if fA < zero(TF)
            sum_A_pos += -fA    # outgoing Δf from cell i → α increases
        else
            sum_A_neg += fA     # incoming Δf into cell i → α decreases
        end
    end

    alpha_star = alpha_prev[i] - dt / volume * sum_L

    # P±: total |anti-diffusive| flux tending to raise / lower α in this cell
    Pplus[i]  = sum_A_pos
    Pminus[i] = sum_A_neg

    # Q±: room to bounds [0,1]. Clamp to zero — if α* has already overshot the
    # bounds (e.g. due to non-orthogonal correction or boundary flux), R± would
    # otherwise go negative and the limiter would amplify Δf instead of damping.
    Qplus[i]  = max(zero(TF), (one(TF) - alpha_star)) * volume / dt
    Qminus[i] = max(zero(TF), alpha_star) * volume / dt
end

@kernel inbounds=true function _mules_ratios!(Pplus, Pminus, Qplus, Qminus, Rplus, Rminus)
    i = @index(Global)
    TF = eltype(Pplus.values)

    Pp = Pplus[i]; Pm = Pminus[i]
    Qp = Qplus[i]; Qm = Qminus[i]

    # R± ∈ [0, 1]: scale the anti-diffusive contribution so a cell cannot exceed
    # its remaining headroom. R± = 0 when α* is already at/past the bound.
    Rplus[i]  = Pp > eps(TF) ? clamp(Qp / Pp, zero(TF), one(TF)) : one(TF)
    Rminus[i] = Pm > eps(TF) ? clamp(Qm / Pm, zero(TF), one(TF)) : one(TF)
end

@kernel inbounds=true function _mules_apply!(phiAf, faces, Rplus, Rminus, nbfaces)
    i = @index(Global)
    TF = eltype(phiAf.values)

    if i > nbfaces
        face = faces[i]
        (; ownerCells) = face
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]

        fA = phiAf[i]
        # Flux sign convention: phiAf[i] > 0 means α transported from owner(1) → owner(2).
        # Owner gives α → must respect its R−. Neighbour receives α → must respect its R+.
        lambda = if fA > zero(TF)
            min(Rplus[cID2], Rminus[cID1])
        elseif fA < zero(TF)
            min(Rplus[cID1], Rminus[cID2])
        else
            one(TF)
        end

        phiAf[i] = lambda * fA
    end
end


## Laplacian smoothing for curvature computation

# Perform `n_passes` face-area-weighted Laplacian smoothing of `phi` into
# `phi_smooth`, using `tmp` as a ping-pong buffer.
function laplacian_smooth!(phi_smooth, tmp, phi, n_passes, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = phi.mesh
    (; cells, cell_faces, cell_nsign, faces) = mesh
    ndrange = length(cells)

    # Initialise smooth field from raw field
    @. phi_smooth.values = phi.values

    kernel! = _laplacian_smooth_pass!(_setup(backend, workgroup, ndrange)...)
    for _ in 1:n_passes
        kernel!(tmp, phi_smooth, cells, cell_faces, faces)
        @. phi_smooth.values = tmp.values
    end
end

@kernel inbounds=true function _laplacian_smooth_pass!(
        out, phi, cells::AbstractArray{Cell{TF,SV,UR}}, cell_faces, faces) where {TF,SV,UR}
    i = @index(Global)
    @inbounds begin
        (; faces_range) = cells[i]
        sum_w   = zero(TF)
        sum_wa  = zero(TF)
        for fi ∈ faces_range
            fID = cell_faces[fi]
            area = faces[fID].area
            # Use the neighbour cell value for internal faces; for boundary faces
            # ownerCells[2] == ownerCells[1], so the centre value is used (no-flux).
            oc = faces[fID].ownerCells
            j  = ifelse(oc[1] == i, oc[2], oc[1])
            sum_wa += area * phi.values[j]
            sum_w  += area
        end
        out.values[i] = sum_wa / sum_w
    end
end


## CSF surface tension helpers

function nhat_prep!(nhatf_prep, alpha, ∇alphaf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = alpha.mesh
    faces = mesh.faces
    nbfaces = length(mesh.boundary_cellsID)
    nfaces  = length(faces)

    # Zero boundary faces so div(n̂) does not receive a biased boundary injection.
    if nbfaces > 0
        kernel! = _nhat_zero_bfaces!(_setup(backend, workgroup, nbfaces)...)
        kernel!(nhatf_prep)
    end

    # Normalise internal faces only.
    ninternal = nfaces - nbfaces
    if ninternal > 0
        kernel! = _nhat_normalise_ifaces!(_setup(backend, workgroup, ninternal)...)
        kernel!(nhatf_prep, faces, ∇alphaf, nbfaces)
    end
end

@kernel inbounds=true function _nhat_zero_bfaces!(nhatf_prep)
    i = @index(Global)
    nhatf_prep[i] = SVector(0.0, 0.0, 0.0)
end

@kernel inbounds=true function _nhat_normalise_ifaces!(nhatf_prep, faces, ∇alphaf_, nbfaces)
    i = @index(Global)
    fID = i + nbfaces
    face = faces[fID]
    (; area) = face

    g = ∇alphaf_[fID]
    mag = norm(g)

    deltaN = 1e-8 / area

    nhatf_prep[fID] = mag < deltaN ? SVector(0.0, 0.0, 0.0) : g / mag
end

function surface_tension_flux!(rDf, sigma, kappaf, alpha, ∇alphaf, phi_gf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    faces = phi_gf.mesh.faces

    ndrange = length(phi_gf)
    kernel! = _surface_tension_flux!(_setup(backend, workgroup, ndrange)...)
    kernel!(rDf, sigma, kappaf, alpha, ∇alphaf, phi_gf, faces)
end
@kernel inbounds=true function _surface_tension_flux!(rDf, sigma, kappaf, alpha, ∇alphaf, phi_gf, faces)
    i = @index(Global)
    face = faces[i]
    (; area, normal, ownerCells, delta) = face
    Sf = area * normal

    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    alpha1 = alpha[cID1]
    alpha2 = alpha[cID2]

    ∇alphaf_vec = normal * ((alpha2 - alpha1) / delta)

    phi_gf[i] -= sigma * kappaf[i] * (∇alphaf_vec ⋅ Sf) * rDf[i]
end


## Reconstruction: FaceScalarField -> cell-centred VectorField (least-squares)

function reconstruct_operation!(phi::VectorField, psif::FaceScalarField, config)
    mesh = phi.mesh
    (; cells, cell_nsign, cell_faces, faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    F = _get_float(mesh)

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
            (; area, normal) = faces[fID]
            nx = normal[1]; ny = normal[2]

            m11 += area * nx * nx
            m12 += area * nx * ny
            m22 += area * ny * ny

            ssf = psif[fID]
            b1 += nx * ssf
            b2 += ny * ssf
        end

        det = m11*m22 - m12*m12
        invdet = abs(det) > eps(TF) ? one(TF)/det : zero(TF)

        ux = ( m22*b1 - m12*b2) * invdet
        uy = (-m12*b1 + m11*b2) * invdet

        phi[i] = @SVector [ux, uy, zero(TF)]
    end
end
