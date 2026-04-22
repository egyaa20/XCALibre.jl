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
    mueff = FaceScalarField(mesh)
    divHv = ScalarField(mesh)

    phi_g = VectorField(mesh)
    phi_gf = FaceScalarField(mesh)

    slip_momentum_term = FaceTensorField(mesh)
    div_slip_momentum  = VectorField(mesh)

    # Explicit α transport via MULES (Zalesak FCT), adapted from the VOF solver.
    # HO flux = mdotf·αf_HO (van Leer) + drift (−αf·(1−αf)·Urdotf).
    # LO flux = mdotf·αf_upwind (bounded, monotone).
    # Anti-diffusive Δf = HO − LO is face-limited via MULES so α stays in
    # [αMin_local, αMax_local] ⊂ [0,1] without ad-hoc clamping.
    alpha_prev      = ScalarField(mesh)
    div_alpha       = ScalarField(mesh)
    div_mdotf       = ScalarField(mesh)
    alpha_fluxf     = FaceScalarField(mesh)
    alphaf_upwind   = FaceScalarField(mesh)
    alphaf_HO       = FaceScalarField(mesh)
    phiLf           = FaceScalarField(mesh)
    phiHf           = FaceScalarField(mesh)
    phiAf           = FaceScalarField(mesh)
    Pplus           = ScalarField(mesh)
    Pminus          = ScalarField(mesh)
    Qplus           = ScalarField(mesh)
    Qminus          = ScalarField(mesh)
    Rplus           = ScalarField(mesh)
    Rminus          = ScalarField(mesh)
    alphaMaxLocal   = ScalarField(mesh)
    alphaMinLocal   = ScalarField(mesh)

    @info "Computing fluid properties..."

    blend_properties!(rho,  alpha,  phases[main].rho[1], phases[secondary].rho[1])
    blend_properties!(rhof, alphaf, phases[main].rho[1], phases[secondary].rho[1])
    blend_properties!(nu,   alpha,  phases[main].mu[1] / phases[main].rho[1],
                                    phases[secondary].mu[1] / phases[secondary].rho[1])
    blend_properties!(nuf,  alphaf, phases[main].mu[1] / phases[main].rho[1],
                                    phases[secondary].mu[1] / phases[secondary].rho[1])
    @. mueff.values = rhof.values * nueff.values

    gh = model.fluid.physics_properties.gravity.gh
    ghf = model.fluid.physics_properties.gravity.ghf
    g = model.fluid.physics_properties.gravity.g

    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)

    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(rho, U)
        + Divergence{schemes.U.divergence}(rhoPhi, U)
        - Laplacian{schemes.U.laplacian}(mueff, U)
        ==
        - Source(∇p_rgh.result)
        - Source(div_slip_momentum)
    ) → VectorEquation(U, boundaries.U)

    p_eqn = (
        - Laplacian{schemes.p.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)

    # alpha_eqn is not used — MULES performs an explicit update below. Kept as
    # a minimal placeholder to avoid touching the solver_variant call signature.
    alpha_eqn = (
        Time{schemes.alpha.time}(alpha)
        + Divergence{schemes.alpha.divergence}(mdotf, alpha)
        ==
        Source(ConstantScalar(0.0))
    ) → ScalarEquation(alpha, boundaries.alpha)

    @info "Initialising preconditioners..."

    @reset U_eqn.preconditioner     = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner     = set_preconditioner(solvers.p_rgh.preconditioner, p_eqn)
    @reset alpha_eqn.preconditioner = set_preconditioner(solvers.alpha.preconditioner, alpha_eqn)

    @info "Pre-allocating solvers..."

    @reset U_eqn.solver     = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver     = _workspace(solvers.p_rgh.solver, _b(p_eqn))
    @reset alpha_eqn.solver = _workspace(solvers.alpha.solver, _b(alpha_eqn))

    @info "Initialising turbulence model..."
    turbulenceModel, config = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals = solver_variant(
        model, turbulenceModel, ∇p, ∇p_rgh, U_eqn, p_eqn, alpha_eqn,
        mdotf, rhoPhi, gh, ghf, phi_g, phi_gf,
        slip_momentum_term, div_slip_momentum,
        alpha_prev, div_alpha, div_mdotf, alpha_fluxf,
        alphaf_upwind, alphaf_HO, phiLf, phiHf, phiAf,
        Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
        alphaMaxLocal, alphaMinLocal,
        config;
        output=output, pref=pref,
        ncorrectors=ncorrectors, inner_loops=inner_loops)

    return residuals
end


function MULTIPHASE(
    model, turbulenceModel, ∇p, ∇p_rgh, U_eqn, p_eqn, alpha_eqn,
    mdotf, rhoPhi, gh, ghf, phi_g, phi_gf,
    slip_momentum_term, div_slip_momentum,
    alpha_prev, div_alpha, div_mdotf, alpha_fluxf,
    alphaf_upwind, alphaf_HO, phiLf, phiHf, phiAf,
    Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
    alphaMaxLocal, alphaMinLocal,
    config;
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
    mueff = get_flux(U_eqn, 3)
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)

    @info "Allocating working memory..."

    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    Uf = FaceVectorField(mesh)
    S = StrainRate(gradU, gradUT, U, Uf)

    ∇p_rghf_deconstructed = FaceScalarField(mesh)
    ∇p_rghf_reconstructed = VectorField(mesh)

    phases = model.fluid.phases
    volume_fraction = model.fluid.volume_fraction
    main = volume_fraction
    secondary = 3 - volume_fraction

    rho1_val = phases[main].rho[1]
    rho2_val = phases[secondary].rho[1]
    mu1_val  = phases[main].mu[1]
    # Multiphase sub-model (VOF or Mixture) carries its own knobs. Pulled from
    # `model.fluid.model` — set via `model = Mixture(diameter=...)` or
    # `model = VOF(sigma=..., cAlpha=...)` on `Fluid{Multiphase}(...)`.
    mp_model = model.fluid.model
    diameter = mp_model isa Mixture ? mp_model.diameter : 1.0e-6

    # ------------------------------------------------------------------
    # Example: CPISO-style dispatch on the multiphase sub-model. Use the
    # same pattern wherever the VOF and Mixture paths need to diverge
    # inside this solver (momentum sources, α transport, flux correction,
    # …). Left here commented for reference — not wired up yet.
    #
    # if typeof(mp_model) <: VOF
    #     # VOF-specific assembly, e.g. surface-tension source
    #     sigma = mp_model.sigma
    #     # U_eqn = (...  - Source(sigma * kappa * ∇alpha)) → VectorEquation(U, boundaries.U)
    # elseif typeof(mp_model) <: Mixture
    #     # Mixture/drift-flux assembly with slip momentum source
    #     # U_eqn = (...  - Source(div_slip_momentum))       → VectorEquation(U, boundaries.U)
    # end
    # ------------------------------------------------------------------
    # Buoyancy form for the drift-velocity closure in compute_Ur!:
    #   true  → standard Manninen  : (ρ_d - ρ_m) / ρ_d
    #   false → legacy form        : (ρ_d - ρ_m) / ρ_m
    # Standard Manninen is the published drift-flux closure and matches the
    # hand-calculated u_t ≈ 21 mm/s for 200 µm glass beads (see
    # test/mixture_testing/README.md). Keep the legacy form only for
    # backward reproduction of earlier results.
    standard_manninen = true
    Sc_t     = 0.7
    g_vec    = model.fluid.physics_properties.gravity.g
    tau_d    = (rho2_val * diameter^2) / (18.0 * mu1_val + eps())

    rho1f = FaceScalarField(mesh)
    rho2f = FaceScalarField(mesh)
    initialise!(rho1f, rho1_val)
    initialise!(rho2f, rho2_val)

    ∇alpha     = Grad{schemes.alpha.gradient}(alpha)
    ∇alphaf    = FaceVectorField(mesh)

    U_prev = VectorField(mesh)
    DUmDt  = VectorField(mesh)
    Ur     = VectorField(mesh)
    Urf    = FaceVectorField(mesh)
    Urdotf = FaceScalarField(mesh)
    ∇U     = Grad{schemes.U.gradient}(U)

    n_cells = length(mesh.cells)
    Hv = VectorField(mesh)
    rD = ScalarField(mesh)
    rho_prev = ScalarField(mesh)

    TF = _get_float(mesh)
    TI = _get_int(mesh)
    prev = KernelAbstractions.zeros(backend, TF, n_cells)

    R_ux    = ones(TF, iterations)
    R_uy    = ones(TF, iterations)
    R_uz    = ones(TF, iterations)
    R_p     = ones(TF, iterations)
    R_alpha = ones(TF, iterations)
    cellsCourant      = KernelAbstractions.zeros(backend, TF, n_cells)
    cellsAlphaCourant = KernelAbstractions.zeros(backend, TF, n_cells)

    time = zero(TF)
    interpolate!(Uf, U, config)
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    flux!(mdotf, Uf, config)

    @. rhoPhi.values = mdotf.values * rhof.values

    update_nueff!(nueff, nuf, model.turbulence, config)
    @. mueff.values = rhof.values * nueff.values

    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    @info "Starting mixture multiphase loops..."

    diag_enable = get(ENV, "XCALIBRE_MMP_DIAG", "0") == "1"
    diag_every  = parse(Int, get(ENV, "XCALIBRE_MMP_DIAG_EVERY", "50"))

    progress = Progress(iterations; dt=1.0, showspeed=true)

    @time for iteration ∈ 1:iterations
        copyto!(dt_cpu, config.runtime.dt)
        time += dt_cpu[1]

        @. alpha_prev.values = alpha.values
        @. U_prev.x.values   = U.x.values
        @. U_prev.y.values   = U.y.values
        @. U_prev.z.values   = U.z.values
        @. rho_prev.values   = rho.values

        grad!(∇U, Uf, U, time, config)
        grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
        limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

        compute_DUmDt!(DUmDt, U, U_prev, ∇U, dt_cpu[1], config)
        compute_Ur!(Ur, alpha, rho, g_vec, DUmDt,
                    rho1_val, rho2_val, mu1_val, diameter, tau_d,
                    standard_manninen, config)
        # initialise!(Ur, [0.0, 0.0, 0.0])

        turbulent_dispersion!(Ur, alpha, ∇alpha, model.turbulence, Sc_t, config)

        interpolate_vanleer!(Urf, Ur, mdotf, config)
        zero_wall_drift_velocity!(Urf, config)
        face_dot_Sf!(Urdotf, Urf, config)

        # --- MULES-bounded explicit α transport --------------------------
        # HO flux includes the drift: −αf·(1−αf)·Urdotf on top of mdotf·αf_HO.
        # LO flux is plain upwind on mdotf only — monotone, bounded.
        # Anti-diffusive Δf = HO − LO is Zalesak-limited so α stays in the
        # per-cell neighbour bounds ∩ [0,1] without a global clamp.
        interpolate_upwind!(alphaf_upwind, alpha, mdotf, config)
        correct_boundaries!(alphaf_upwind, alpha, boundaries.alpha, time, config)
        interpolate_vanleer!(alphaf_HO, alpha, ∇alpha, mdotf, config)
        correct_boundaries!(alphaf_HO, alpha, boundaries.alpha, time, config)

        @. phiLf.values = mdotf.values * alphaf_upwind.values
        # Drift contribution to HO flux: tracked-α equation has −∇·(α(1−α)Ur_d).
        # Use alphaf_upwind for the α(1−α) factor (interFOAM convention for
        # the compression-style anti-diffusive term).
        @. phiHf.values = mdotf.values * alphaf_HO.values -
                          Urdotf.values * alphaf_upwind.values *
                          (1.0 - alphaf_upwind.values)
        @. phiAf.values = phiHf.values - phiLf.values

        # Zalesak can't limit boundary faces (single owner); let the LO flux
        # pass through unmodified there.
        zero_boundary_faces!(phiAf, config)

        mules_limit!(phiAf, alpha_prev, phiLf,
                     Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
                     alphaMaxLocal, alphaMinLocal,
                     dt_cpu[1], mesh, config)

        @. alpha_fluxf.values = phiLf.values + phiAf.values
        div!(div_alpha, alpha_fluxf, config)
        div!(div_mdotf, mdotf, config)
        # div-corrected explicit update: α^{n+1} = α^n − dt/V · (div(F) − α^n·div(mdotf))
        @. alpha.values = alpha_prev.values -
            dt_cpu[1] * (div_alpha.values - alpha_prev.values * div_mdotf.values)

        alpha_min_pre = minimum(alpha.values)
        alpha_max_pre = maximum(alpha.values)
        # No hard clamp — MULES keeps α bounded. Re-interpolate αf from the
        # new α with van Leer for downstream density blending.
        interpolate_vanleer!(alphaf, alpha, ∇alpha, mdotf, config)
        correct_boundaries!(alphaf, alpha, boundaries.alpha, time, config)

        if diag_enable
            mules_bc_diagnostics!(
                iteration, diag_every, alphaf, alpha, Urdotf,
                boundaries.alpha, mesh)
        end

        ralpha = zero(TF)

        blend_properties!(rho,  alpha,  rho1_val, rho2_val)
        blend_properties!(rhof, alphaf, rho1_val, rho2_val)
        blend_properties!(nu,   alpha,  mu1_val / rho1_val,
                                        phases[secondary].mu[1] / rho2_val)
        blend_properties!(nuf,  alphaf, mu1_val / rho1_val,
                                        phases[secondary].mu[1] / rho2_val)
        @. mueff.values = rhof.values * nueff.values

        compute_tensor_term!(alphaf, rhof, rho1f, rho2f, Urf, slip_momentum_term, config)
        div_tensor!(div_slip_momentum, slip_momentum_term, config)

        @. rhoPhi.values = mdotf.values * rhof.values

        rx, ry, rz = solve_equation!(
            U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config, rho_prev; time=time)

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
            reconstruct_operation!(phi_g, phi_gf, config)

            update_fixedfluxpressure_gradient!(boundaries.p_rgh, phi_gf, rDf, mesh)

            @. mdotf.values += phi_gf.values

            div!(divHv, mdotf, config)

            @. prev = p_rgh.values
            rp = solve_equation!(p_eqn, p_rgh, boundaries.p_rgh, solvers.p_rgh,
                                 config, rho_prev; ref=pref, time=time)

            grad!(∇p_rgh, p_rghf, p_rgh, boundaries.p_rgh, time, config)
            limit_gradient!(schemes.p_rgh.limiter, ∇p_rgh, p_rgh, config)

            correct_mass_flux1(mdotf, p_eqn, config)

            pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)
            reconstruct_operation!(∇p_rghf_reconstructed, ∇p_rghf_deconstructed, config)
            correct_velocity_rgh!(U, Hv, ∇p_rghf_reconstructed, rD, phi_g, config)
        end

        @. p.values = p_rgh.values + (rho.values * gh.values)

        turbulence!(turbulenceModel, model, S, prev, time, config)
        update_nueff!(nueff, nuf, model.turbulence, config)
        @. mueff.values = rhof.values * nueff.values

        courant      = max_courant_number!(cellsCourant, model, config)
        alphaCourant = max_alpha_courant_number!(cellsAlphaCourant, alpha, mdotf, model, config, dt_cpu[1])
        update_dt!(config.runtime, courant, alphaCourant)

        if diag_enable
            multiphase_diagnostics!(
                iteration, diag_every, time, dt_cpu[1],
                alpha, alpha_min_pre, alpha_max_pre,
                U, rho, mdotf, phi_gf, Ur,
                boundaries.U, boundaries.alpha, mesh)
        end

        R_ux[iteration]    = rx
        R_uy[iteration]    = ry
        R_uz[iteration]    = rz
        R_p[iteration]     = rp
        R_alpha[iteration] = ralpha

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
                (:alpha, R_alpha[iteration]),
                turbulenceModel.state.residuals...
                ]
            )

        runtime_postprocessing!(postprocess, iteration, iterations)

        if iteration % write_interval + signbit(write_interval) == 0
            save_output(model, outputWriter, iteration, time, config)
            save_postprocessing(postprocess, iteration, time, mesh, outputWriter, config.boundaries)
        end
    end

    return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p, alpha=R_alpha)
end


function compute_DUmDt!(DUmDt, U, U_prev, gradU, dt, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(DUmDt)
    kernel! = _compute_DUmDt!(_setup(backend, workgroup, ndrange)...)
    kernel!(DUmDt, U, U_prev, gradU.result, dt)
end

@kernel inbounds=true function _compute_DUmDt!(DUmDt, U, U_prev, gradU_result, dt)
    i = @index(Global)
    TF = eltype(U.x)

    dUdt = (U[i] - U_prev[i]) / dt

    ux = U.x[i]; uy = U.y[i]; uz = U.z[i]
    dudx = gradU_result.xx[i]; dudy = gradU_result.xy[i]; dudz = gradU_result.xz[i]
    dvdx = gradU_result.yx[i]; dvdy = gradU_result.yy[i]; dvdz = gradU_result.yz[i]
    dwdx = gradU_result.zx[i]; dwdy = gradU_result.zy[i]; dwdz = gradU_result.zz[i]

    conv_x = ux*dudx + uy*dudy + uz*dudz
    conv_y = ux*dvdx + uy*dvdy + uz*dvdz
    conv_z = ux*dwdx + uy*dwdy + uz*dwdz

    DUmDt[i] = dUdt + @SVector [conv_x, conv_y, conv_z]
end

function compute_Ur!(Ur, alpha, rho, g, DUmDt, rho1, rho2, mu1, d, tau_d,
                     standard_manninen::Bool, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(Ur)
    kernel! = _compute_Ur!(_setup(backend, workgroup, ndrange)...)
    kernel!(Ur, alpha, rho, g, DUmDt, rho1, rho2, mu1, d, tau_d, standard_manninen)
end

@kernel inbounds=true function _compute_Ur!(Ur, alpha, rho, g, DUmDt, rho1, rho2, mu1, d, tau_d, standard_manninen)
    i = @index(Global)
    TF = eltype(rho.values)

    rho_m = rho[i]
    rho_c = rho1
    rho_d = rho2
    mu_c  = mu1

    Ur_mag = norm(Ur[i])
    Re_p   = rho_c * Ur_mag * d / (mu_c + eps(TF))

    f_drag = ifelse(
        Re_p < TF(1000),
        one(TF) + TF(0.15) * Re_p^TF(0.687),
        TF(0.0183) * Re_p
    )
    f_drag = max(f_drag, one(TF))

    a_eff    = g - DUmDt[i]
    # Standard Manninen drift closure uses ρ_d in the denominator (published
    # form). Legacy form used ρ_m (mixture density) and inflated u_t by ~1.7×
    # at dilute limit for heavy dispersed phase.
    buoyancy = ifelse(
        standard_manninen,
        (rho_d - rho_m) / (rho_d + eps(TF)),
        (rho_d - rho_m) / (rho_m + eps(TF))
    )

    U_dm = (tau_d / (f_drag + eps(TF))) * buoyancy * a_eff

    alpha_c = max(alpha[i], TF(1e-3))
    Ur[i] = U_dm / alpha_c
end

turbulent_dispersion!(Ur, alpha, ∇alpha, turbulence::Laminar, Sc_t, config) = nothing

function turbulent_dispersion!(Ur, alpha, ∇alpha, turbulence, Sc_t, config)
    hasproperty(turbulence, :nut) || return nothing
    (; hardware) = config
    (; backend, workgroup) = hardware
    nut = turbulence.nut

    ndrange = length(Ur)
    kernel! = _turbulent_dispersion!(_setup(backend, workgroup, ndrange)...)
    kernel!(Ur, alpha, ∇alpha.result, nut, Sc_t)
end

@kernel inbounds=true function _turbulent_dispersion!(Ur, alpha, gradA, nut, Sc_t)
    i = @index(Global)
    TF = eltype(alpha.values)

    α_c = alpha[i]
    α_c_safe = max(α_c, TF(1e-3))
    α_d_safe = max(one(TF) - α_c, TF(1e-3))

    D_t   = nut[i] / TF(Sc_t)
    denom = α_c_safe * α_d_safe + eps(TF)
    coef  = D_t / denom

    gx = gradA.x[i]; gy = gradA.y[i]; gz = gradA.z[i]

    Ur[i] = Ur[i] + @SVector [coef*gx, coef*gy, coef*gz]
end

function compute_tensor_term!(alphaf, rhof, rho1f, rho2f, Urf, T, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(alphaf)
    kernel! = _compute_tensor_term!(_setup(backend, workgroup, ndrange)...)
    kernel!(alphaf, rhof, rho1f, rho2f, Urf, T)
end

@kernel inbounds=true function _compute_tensor_term!(alphaf, rhof, rho1f, rho2f, Urf, T)
    i = @index(Global)
    TF = eltype(alphaf)

    af = alphaf[i]
    ρm = rhof[i]
    ρ1 = rho1f[i]
    ρ2 = rho2f[i]

    ux = Urf.x[i]; uy = Urf.y[i]; uz = Urf.z[i]

    coeff = (af * (one(TF) - af) * (ρ1 + ρ2)) / (ρm + eps(TF))

    T.xx[i] = coeff * ux * ux
    T.xy[i] = coeff * ux * uy
    T.xz[i] = coeff * ux * uz

    T.yx[i] = coeff * uy * ux
    T.yy[i] = coeff * uy * uy
    T.yz[i] = coeff * uy * uz

    T.zx[i] = coeff * uz * ux
    T.zy[i] = coeff * uz * uy
    T.zz[i] = coeff * uz * uz
end

function build_drift_alpha_flux!(flux, Urf, alphaf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(alphaf)
    kernel! = _build_drift_alpha_flux!(_setup(backend, workgroup, ndrange)...)
    kernel!(flux, Urf, alphaf)
end

@kernel inbounds=true function _build_drift_alpha_flux!(flux, Urf, alphaf)
    i = @index(Global)
    TF = eltype(alphaf)
    af = alphaf[i]
    w  = -af * (one(TF) - af)
    flux.x[i] = Urf.x[i] * w
    flux.y[i] = Urf.y[i] * w
    flux.z[i] = Urf.z[i] * w
end

function face_dot_Sf!(phidotf, phif, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh  = phidotf.mesh
    faces = mesh.faces

    ndrange = length(faces)
    kernel! = _face_dot_Sf!(_setup(backend, workgroup, ndrange)...)
    kernel!(phidotf, phif, faces)
end

@kernel inbounds=true function _face_dot_Sf!(phidotf, phif, faces)
    i = @index(Global)
    (; area, normal) = faces[i]
    phidotf[i] = area * (phif.x[i]*normal[1] + phif.y[i]*normal[2] + phif.z[i]*normal[3])
end

function compute_U_eff_alpha_dotSf!(U_eff_dotSf, alpha, Urdotf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh  = U_eff_dotSf.mesh
    faces = mesh.faces
    avals = alpha.values

    ndrange = length(faces)
    kernel! = _compute_U_eff_alpha_dotSf!(_setup(backend, workgroup, ndrange)...)
    kernel!(U_eff_dotSf, avals, Urdotf, faces)
end

@kernel inbounds=true function _compute_U_eff_alpha_dotSf!(U_eff_dotSf, avals, Urdotf, faces)
    i = @index(Global)
    TF = eltype(U_eff_dotSf)
    face = faces[i]
    (; ownerCells) = face
    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    Urdn = Urdotf[i]
    # Upwind (1-α) by sign of Urdotf. At boundary faces cID1==cID2 and
    # Urdn=0 (set by zero_wall_drift_velocity!), so U_eff_dotSf=0 trivially.
    # Sign matches original drift_alpha_flux = -αf(1-αf)·Urf: the LHS term
    # is +div(U_eff·α) with U_eff = -(1-α_up)·Urdotf.
    α_up = Urdn >= zero(TF) ? avals[cID1] : avals[cID2]
    U_eff_dotSf[i] = -(one(TF) - α_up) * Urdn
end

function zero_wall_drift_velocity!(Urf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    nbfaces = length(Urf.mesh.boundary_cellsID)
    if nbfaces > 0
        ndrange = nbfaces
        kernel! = _zero_wall_drift_velocity!(_setup(backend, workgroup, ndrange)...)
        kernel!(Urf)
    end
end

@kernel inbounds=true function _zero_wall_drift_velocity!(Urf)
    i = @index(Global)
    TF = eltype(Urf.x)
    Urf.x[i] = zero(TF)
    Urf.y[i] = zero(TF)
    Urf.z[i] = zero(TF)
end


function multiphase_diagnostics!(
    iteration, every, time, dt,
    alpha, alpha_min_pre, alpha_max_pre,
    U, rho, mdotf, phi_gf, Ur,
    U_BCs, alpha_BCs, mesh)

    overshoot_lo = alpha_min_pre < -1e-6
    overshoot_hi = alpha_max_pre >  1.0 + 1e-6
    should_print = (iteration % every == 0) || (iteration == 1)
    should_print || return

    mesh_boundaries = mesh.boundaries isa AbstractArray ?
                      Array(mesh.boundaries) : mesh.boundaries
    mvals   = mdotf.values
    pgfvals = phi_gf.values
    avals   = alpha.values
    boundary_cellsID = mesh.boundary_cellsID
    F = eltype(mvals)

    amin = minimum(avals); amax = maximum(avals)
    rmin = minimum(rho.values); rmax = maximum(rho.values)
    Uxm  = maximum(abs, U.x.values)
    Uym  = maximum(abs, U.y.values)
    Uzm  = maximum(abs, U.z.values)
    Urxm = maximum(abs, Ur.x.values)
    Urym = maximum(abs, Ur.y.values)
    Urzm = maximum(abs, Ur.z.values)

    header = "=== MMP DIAG iter=$iteration t=$(round(time,sigdigits=4)) dt=$(round(dt,sigdigits=3)) ==="
    lines = String[header]

    push!(lines, "α (pre-clamp): min=$(round(alpha_min_pre,sigdigits=6)) max=$(round(alpha_max_pre,sigdigits=6))")
    push!(lines, "α (post-clamp): min=$(round(amin,sigdigits=6)) max=$(round(amax,sigdigits=6))")
    push!(lines, "ρ: [$(round(rmin,sigdigits=5)), $(round(rmax,sigdigits=5))]  |U|max=(x:$(round(Uxm,sigdigits=3)), y:$(round(Uym,sigdigits=3)), z:$(round(Uzm,sigdigits=3)))  |Ur|max=(x:$(round(Urxm,sigdigits=3)), y:$(round(Urym,sigdigits=3)), z:$(round(Urzm,sigdigits=3)))")

    if overshoot_lo || overshoot_hi
        push!(lines, "  [WARNING] alpha overshoot outside [0,1] pre-clamp — numerical unboundedness in alpha transport")
    end

    total_in  = zero(F); total_out = zero(F)
    for BC in U_BCs
        name = mesh_boundaries[BC.ID].name
        net  = zero(F); nrev = 0; nfwd = 0
        pgsum = zero(F); pgmax = zero(F)
        amin_b = F(Inf); amax_b = F(-Inf)
        @inbounds for fID in BC.IDs_range
            v = mvals[fID]
            net += v
            if v < 0; nrev += 1; else; nfwd += 1; end
            p = pgfvals[fID]; pgsum += p
            ap = abs(p); ap > pgmax && (pgmax = ap)
            cID = boundary_cellsID[fID]
            av = avals[cID]
            av < amin_b && (amin_b = av)
            av > amax_b && (amax_b = av)
        end
        if net < 0; total_in += net; else; total_out += net; end
        push!(lines, "  :$name  mdot_net=$(round(net,sigdigits=4)) (rev:$nrev/fwd:$nfwd)  φ_g[sum=$(round(pgsum,sigdigits=3)) |max|=$(round(pgmax,sigdigits=3))]  α_adj=[$(round(amin_b,sigdigits=4)), $(round(amax_b,sigdigits=4))]")
    end
    push!(lines, "  TOTALS  in=$(round(total_in,sigdigits=4))  out=$(round(total_out,sigdigits=4))  imbalance=$(round(total_in+total_out,sigdigits=3))")

    @info join(lines, "\n")
end

function mules_bc_diagnostics!(
    iteration, every, alphaf, alpha, Urdotf, alpha_BCs, mesh)

    should_print = (iteration % every == 0) || (iteration == 1)
    should_print || return

    mesh_boundaries = mesh.boundaries isa AbstractArray ?
                      Array(mesh.boundaries) : mesh.boundaries
    boundary_cellsID = mesh.boundary_cellsID
    avals  = Array(alpha.values)
    afvals = Array(alphaf.values)
    urvals = Array(Urdotf.values)
    F = eltype(avals)

    lines = String["=== MMP BC-DIAG iter=$iteration ==="]
    push!(lines, "  α-face vs α-owner-cell deviation per patch")
    push!(lines, "  (≈0 ⇒ zero-gradient BC applied correctly; ≈αc ⇒ face defaulting to 0)")
    for BC in alpha_BCs
        name = mesh_boundaries[BC.ID].name
        bctype = nameof(typeof(BC))
        dmax = zero(F); dsum = zero(F); n = 0
        afmin = F(Inf); afmax = F(-Inf)
        acmin = F(Inf); acmax = F(-Inf)
        @inbounds for fID in BC.IDs_range
            cID = boundary_cellsID[fID]
            af = afvals[fID]; ac = avals[cID]
            d  = abs(af - ac)
            d > dmax && (dmax = d); dsum += d; n += 1
            af < afmin && (afmin = af); af > afmax && (afmax = af)
            ac < acmin && (acmin = ac); ac > acmax && (acmax = ac)
        end
        davg = n > 0 ? dsum / n : zero(F)
        push!(lines,
            "  :$name [$bctype]  |αf-αc| max=$(round(dmax,sigdigits=4)) avg=$(round(davg,sigdigits=4))  "
            * "αf∈[$(round(afmin,sigdigits=4)),$(round(afmax,sigdigits=4))]  "
            * "αc∈[$(round(acmin,sigdigits=4)),$(round(acmax,sigdigits=4))]")
    end

    push!(lines, "  |Urdotf| per patch (should = 0 at walls after zero_wall_drift_velocity!)")
    for BC in alpha_BCs
        name = mesh_boundaries[BC.ID].name
        urmax = zero(F); ursum = zero(F); n = 0
        @inbounds for fID in BC.IDs_range
            u = abs(urvals[fID])
            u > urmax && (urmax = u); ursum += u; n += 1
        end
        uravg = n > 0 ? ursum / n : zero(F)
        push!(lines,
            "  :$name  |Urdotf| max=$(round(urmax,sigdigits=4)) avg=$(round(uravg,sigdigits=4))")
    end

    @info join(lines, "\n")
end

function update_fixedfluxpressure_gradient!(p_BCs, phi_gf, rDf, mesh)
    phi_vals = phi_gf.values
    rDf_vals = rDf.values
    faces_cpu = mesh.faces
    for BC in p_BCs
        if BC isa FixedFluxPressure
            grad = BC.value
            start = BC.IDs_range.start
            @inbounds for fID in BC.IDs_range
                area = faces_cpu[fID].area
                denom = rDf_vals[fID]*area
                grad[fID - start + 1] = phi_vals[fID]/denom
            end
        end
    end
end

function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0 * alpha_field.values) + (property_1 * (1.0 - alpha_field.values))
    nothing
end

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


## -------------------------------------------------------------------
## MULES (Zalesak FCT) — ported from Solvers_2_MULTIPHASE_vof.jl.
## Adapted here with drift flux in phiHf instead of compression flux.
## -------------------------------------------------------------------

function zero_boundary_faces!(phif::FaceScalarField, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = phif.mesh
    nbfaces = length(mesh.boundary_cellsID)
    if nbfaces > 0
        ndrange = nbfaces
        kernel! = _mmp_zero_boundary_faces!(_setup(backend, workgroup, ndrange)...)
        kernel!(phif)
    end
end
@kernel inbounds=true function _mmp_zero_boundary_faces!(phif)
    i = @index(Global)
    phif[i] = zero(eltype(phif.values))
end

function mules_limit!(phiAf, alpha_prev, phiLf,
                      Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
                      alphaMaxLocal, alphaMinLocal,
                      dt, mesh, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    (; cells, cell_nsign, cell_faces, faces) = mesh

    n_cells  = length(cells)
    n_faces  = length(faces)
    nbfaces  = length(mesh.boundary_cellsID)

    fill!(Pplus.values,  0)
    fill!(Pminus.values, 0)

    ndrange = n_cells
    kernel! = _mmp_mules_local_bounds!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, cell_faces, faces, alpha_prev, alphaMaxLocal, alphaMinLocal)

    ndrange = n_cells
    kernel! = _mmp_mules_cell_accum!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, cell_faces, cell_nsign, faces,
            alpha_prev, phiLf, phiAf,
            Pplus, Pminus, Qplus, Qminus,
            alphaMaxLocal, alphaMinLocal, dt)

    ndrange = n_cells
    kernel! = _mmp_mules_ratios!(_setup(backend, workgroup, ndrange)...)
    kernel!(Pplus, Pminus, Qplus, Qminus, Rplus, Rminus)

    ndrange = n_faces
    kernel! = _mmp_mules_apply!(_setup(backend, workgroup, ndrange)...)
    kernel!(phiAf, faces, Rplus, Rminus, nbfaces)
end

@kernel inbounds=true function _mmp_mules_local_bounds!(
    cells::AbstractArray{Cell{TF,SV,UR}}, cell_faces, faces,
    alpha_prev, alphaMaxLocal, alphaMinLocal
) where {TF,SV,UR}
    # Hard physical bounds [0, 1] for the dispersed/continuous mixture.
    # The VOF-style local-extrema bound (no new local extrema from α_prev
    # neighbours) is wrong here: a uniformly suspended initial state
    # (α = α_c0 everywhere) would pin αMaxLocal = αMinLocal = α_c0 and
    # the limiter would kill every anti-diffusive flux, including the
    # physical drift that drives settling. For a mixture of two miscible-
    # in-volume-fraction phases, [0,1] is the only defensible bound.
    i = @index(Global)
    alphaMaxLocal[i] = one(TF)
    alphaMinLocal[i] = zero(TF)
end

@kernel inbounds=true function _mmp_mules_cell_accum!(
    cells::AbstractArray{Cell{TF,SV,UR}}, cell_faces, cell_nsign, faces,
    alpha_prev, phiLf, phiAf, Pplus, Pminus, Qplus, Qminus,
    alphaMaxLocal, alphaMinLocal, dt
) where {TF,SV,UR}
    i = @index(Global)
    (; volume, faces_range) = cells[i]

    sum_L = zero(TF)
    sum_A_pos = zero(TF)
    sum_A_neg = zero(TF)

    for fi in faces_range
        fID   = cell_faces[fi]
        nsign = cell_nsign[fi]
        fL = phiLf[fID] * nsign
        fA = phiAf[fID] * nsign
        sum_L += fL
        if fA < zero(TF)
            sum_A_pos += -fA
        else
            sum_A_neg += fA
        end
    end

    alpha_star = alpha_prev[i] - dt / volume * sum_L

    Pplus[i]  = sum_A_pos
    Pminus[i] = sum_A_neg

    Qplus[i]  = max(zero(TF), (alphaMaxLocal[i] - alpha_star)) * volume / dt
    Qminus[i] = max(zero(TF), (alpha_star - alphaMinLocal[i])) * volume / dt
end

@kernel inbounds=true function _mmp_mules_ratios!(Pplus, Pminus, Qplus, Qminus, Rplus, Rminus)
    i = @index(Global)
    TF = eltype(Pplus.values)
    Pp = Pplus[i]; Pm = Pminus[i]
    Qp = Qplus[i]; Qm = Qminus[i]
    Rplus[i]  = Pp > eps(TF) ? clamp(Qp / Pp, zero(TF), one(TF)) : one(TF)
    Rminus[i] = Pm > eps(TF) ? clamp(Qm / Pm, zero(TF), one(TF)) : one(TF)
end

@kernel inbounds=true function _mmp_mules_apply!(phiAf, faces, Rplus, Rminus, nbfaces)
    i = @index(Global)
    TF = eltype(phiAf.values)
    if i > nbfaces
        face = faces[i]
        (; ownerCells) = face
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        fA = phiAf[i]
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
