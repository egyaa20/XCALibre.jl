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

multiphase_extras(::VOF, mesh) = ()
function multiphase_extras(::Mixture, mesh)
    slip_momentum_term = FaceTensorField(mesh)
    div_slip_momentum  = VectorField(mesh)
    return (slip_momentum_term, div_slip_momentum)
end

function setup_multiphase_solvers(
    solver_variant, model, config;
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    )

    (; solvers, schemes, runtime, hardware, boundaries) = config

    @info "Extracting configuration and input fields..."

    (; U, p) = model.momentum
    (; alpha, alphaf, rho, rhof, nu, nuf, p_rgh, p_rghf) = model.fluid

    mp_model = model.fluid.model

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

    extra_models = multiphase_extras(mp_model, mesh)

    mules = (
        alpha_prev    = ScalarField(mesh),
        div_alpha     = ScalarField(mesh),
        div_mdotf     = ScalarField(mesh),
        alpha_fluxf   = FaceScalarField(mesh),
        alphaf_upwind = FaceScalarField(mesh),
        alphaf_HO     = FaceScalarField(mesh),
        phiLf         = FaceScalarField(mesh),
        phiHf         = FaceScalarField(mesh),
        phiAf         = FaceScalarField(mesh),
        Pplus         = ScalarField(mesh),
        Pminus        = ScalarField(mesh),
        Qplus         = ScalarField(mesh),
        Qminus        = ScalarField(mesh),
        Rplus         = ScalarField(mesh),
        Rminus        = ScalarField(mesh),
        alphaMaxLocal = ScalarField(mesh),
        alphaMinLocal = ScalarField(mesh),
    )

    @info "Computing fluid properties..."

    blend_properties!(rho,  alpha,  phases[main].rho[1], phases[secondary].rho[1])
    blend_properties!(rhof, alphaf, phases[main].rho[1], phases[secondary].rho[1])
    blend_mixture_nu!(nu,  alpha,  rho,  phases[main].mu[1], phases[secondary].mu[1])
    blend_mixture_nu!(nuf, alphaf, rhof, phases[main].mu[1], phases[secondary].mu[1])
    @. mueff.values = rhof.values * nueff.values

    gh = model.fluid.physics_properties.gravity.gh
    ghf = model.fluid.physics_properties.gravity.ghf
    g = model.fluid.physics_properties.gravity.g

    compute_gh!(gh, g, config)
    compute_ghf!(ghf, g, config)

    @info "Defining models..."

    if typeof(mp_model) <: Mixture

        div_slip_momentum = extra_models[2]

        U_eqn = (
            Time{schemes.U.time}(rho, U)
            + Divergence{schemes.U.divergence}(rhoPhi, U)
            - Laplacian{schemes.U.laplacian}(mueff, U)
            ==
            - Source(∇p_rgh.result)
            - Source(div_slip_momentum)
        ) → VectorEquation(U, boundaries.U)

    elseif typeof(mp_model) <: VOF

        U_eqn = (
            Time{schemes.U.time}(rho, U)
            + Divergence{schemes.U.divergence}(rhoPhi, U)
            - Laplacian{schemes.U.laplacian}(mueff, U)
            ==
            - Source(∇p_rgh.result)
        ) → VectorEquation(U, boundaries.U)

    end

    p_eqn = (
        - Laplacian{schemes.p.laplacian}(rDf, p_rgh)
        ==
        - Source(divHv)
    ) → ScalarEquation(p_rgh, boundaries.p_rgh)

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

    @info "Initialising energy model..."
    energyModel = init_multiphase_energy(model.energy, model, mdotf, rho, p_eqn, config)

    residuals = solver_variant(
        model, turbulenceModel, energyModel, ∇p, ∇p_rgh, U_eqn, p_eqn, alpha_eqn,
        mdotf, rhoPhi, gh, ghf, phi_g, phi_gf,
        extra_models, mules,
        config;
        output=output, pref=pref,
        ncorrectors=ncorrectors, inner_loops=inner_loops)

    return residuals
end

# Energy-model dispatch helpers — keep the multiphase loop agnostic to whether
# the user picked Isothermal (no energy equation) or MultiphaseTemperature.
init_multiphase_energy(::Isothermal, model, mdotf, rho, p_eqn, config) = nothing
init_multiphase_energy(e::MultiphaseTemperature, model, mdotf, rho, p_eqn, config) =
    initialise(e, model, mdotf, rho, p_eqn, config)

step_multiphase_energy!(::Nothing, model, prev, mdotf, rho, time, config) = nothing
step_multiphase_energy!(em::MultiphaseTemperatureModel, model, prev, mdotf, rho, time, config) =
    energy!(em, model, prev, mdotf, rho, time, config)

# Phase-2 live property refresh. Only fires when the user has both a
# `HelmholtzTable` driving properties AND a `MultiphaseTemperature`
# energy model providing a per-cell T field. Without one of those, the
# function is a no-op.
function live_update_phase_properties!(model)
    fp = get(model.fluid.physics_properties, :fluid_properties, nothing)
    fp isa HelmholtzTable || return nothing
    energy = model.energy
    energy isa MultiphaseTemperature || return nothing
    update_phase_properties_from_table!(model.fluid.phases, fp, energy.T)
    return nothing
end


function MULTIPHASE(
    model, turbulenceModel, energyModel, ∇p, ∇p_rgh, U_eqn, p_eqn, alpha_eqn,
    mdotf, rhoPhi, gh, ghf, phi_g, phi_gf,
    extra_models::Tuple, mules::NamedTuple,
    config;
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=3
    )

    (; alpha_prev, div_alpha, div_mdotf, alpha_fluxf,
       alphaf_upwind, alphaf_HO, phiLf, phiHf, phiAf,
       Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
       alphaMaxLocal, alphaMinLocal) = mules

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
    # Buffer for the well-balanced predictor pressure-gradient reconstruction
    pressure_force_face   = FaceScalarField(mesh)

    phases = model.fluid.phases
    volume_fraction = model.fluid.volume_fraction
    main = volume_fraction
    secondary = 3 - volume_fraction

    rho1_val = phases[main].rho[1]
    rho2_val = phases[secondary].rho[1]
    mu1_val  = phases[main].mu[1]
    mu2_val  = phases[secondary].mu[1]

    # Phase-2B field mode flag — set when `HelmholtzTable` drove
    # `build_phase_table_mode` so per-phase rho/mu live on ScalarField
    # cells. In this mode the property kernels read the per-cell fields
    # rather than the captured scalars above.
    field_mode = phases[main].rho isa ScalarField

    # Per-phase face fields. Always allocated so the same arithmetic
    # path can be used by the rhoPhi formula and the Mixture tensor
    # term. In constant mode they're filled once with the snapshot
    # scalars; in field mode they're refreshed each iteration via
    # interpolate! from the live cell fields.
    rho1f = FaceScalarField(mesh)
    rho2f = FaceScalarField(mesh)
    mu1f  = FaceScalarField(mesh)
    mu2f  = FaceScalarField(mesh)
    if field_mode
        interpolate!(rho1f, phases[main].rho,      config)
        interpolate!(rho2f, phases[secondary].rho, config)
        interpolate!(mu1f,  phases[main].mu,       config)
        interpolate!(mu2f,  phases[secondary].mu,  config)
    else
        initialise!(rho1f, rho1_val)
        initialise!(rho2f, rho2_val)
        initialise!(mu1f,  mu1_val)
        initialise!(mu2f,  mu2_val)
    end

    mp_model = model.fluid.model

    ∇alpha  = Grad{schemes.alpha.gradient}(alpha)
    ∇alphaf = FaceVectorField(mesh)

    if typeof(mp_model) <: Mixture
        slip_momentum_term, div_slip_momentum = extra_models

        Sc_t     = 0.7
        g_vec    = model.fluid.physics_properties.gravity.g
        diameter = mp_model.diameter
        # Constant-mode τ_d (single scalar). Field mode below allocates a
        # per-cell ScalarField that is refreshed every outer iteration.
        tau_d    = (rho2_val * diameter^2) / (18.0 * mu1_val + eps())
        tau_d_field = ScalarField(mesh)
        initialise!(tau_d_field, tau_d)

        # rho1f / rho2f already allocated above (used by both modes).

        U_prev = VectorField(mesh)
        DUmDt  = VectorField(mesh)
        Ur     = VectorField(mesh)
        Urf    = FaceVectorField(mesh)
        Urdotf = FaceScalarField(mesh)
        ∇U     = Grad{schemes.U.gradient}(U)
    end

    n_alpha_subcycles = 1   # default; overridden inside VOF block from mp_model.cycles

    if typeof(mp_model) <: VOF
        sigma   = mp_model.sigma
        C_alpha = mp_model.cAlpha

        phirf            = FaceScalarField(mesh)
        nhatf_prep       = FaceVectorField(mesh)
        kappa            = ScalarField(mesh)
        kappaf           = FaceScalarField(mesh)
        grad_alpha_mag   = ScalarField(mesh)
        alpha_smooth     = ScalarField(mesh)
        alpha_smooth_tmp = ScalarField(mesh)
        n_smooth         = 0
        use_csf_smoothed = false

        # Compression-normal α smoothing (mp_model.smooth = 0 disables it).
        # ∇alpha_smooth is reused every step; ∇alpha is preserved for the
        # van-Leer HO α flux which still uses the unsmoothed gradient.
        n_smooth_compression = mp_model.smooth
        ∇alpha_smooth        = Grad{schemes.alpha.gradient}(alpha_smooth)

        # MULES sub-cycling (mp_model.cycles = 1 disables it). N > 1 splits
        # the outer dt into N MULES α-updates with dt_sub = dt/N. The U/p
        # coupling stays at full dt; the adaptive-dt controller is informed
        # via alphaCourant / N so the outer dt grows up to N × maxAlphaCo.
        n_alpha_subcycles = mp_model.cycles
    end

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

    println("Piso loops: $(inner_loops)")
    @info "Starting multiphase loops..."
    

    progress = Progress(iterations; dt=1.0, showspeed=true)

    mean_U_csv_path = "mean_U.csv"
    open(mean_U_csv_path, "w") do io
        println(io, "iteration,time,mean_Umag,max_Umag")
    end

    @time for iteration ∈ 1:iterations
        # Optional time-based termination — exits as soon as the simulated
        # time reaches `runtime.t_end`. Skipped when `t_end === nothing`.
        if config.runtime.t_end !== nothing && time >= config.runtime.t_end
            @info "Reached t_end = $(config.runtime.t_end) s after $(iteration-1) iterations — stopping."
            break
        end

        copyto!(dt_cpu, config.runtime.dt)
        time += dt_cpu[1]

        @. rho_prev.values = rho.values

        if typeof(mp_model) <: Mixture
            @. U_prev.x.values = U.x.values
            @. U_prev.y.values = U.y.values
            @. U_prev.z.values = U.z.values

            grad!(∇U, Uf, U, time, config)
        end

        # Mixture drift terms are evaluated once per outer step using the
        # outer α field. They stay constant during α sub-cycling, mirroring
        # mdotf which is also held fixed inside the sub-cycle loop.
        if typeof(mp_model) <: Mixture
            grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
            limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

            compute_DUmDt!(DUmDt, U, U_prev, ∇U, dt_cpu[1], config)
            # rho1/rho2/mu1/tau_d below are either ConstantScalar (legacy)
            # or per-cell ScalarField (Phase-2B field mode). The kernel
            # reads them as `x[i]`, which is well-defined for both.
            tau_d_arg = field_mode ? tau_d_field : ConstantScalar(tau_d)
            compute_Ur!(Ur, alpha, rho, g_vec, DUmDt,
                        phases[main].rho, phases[secondary].rho, phases[main].mu,
                        diameter, tau_d_arg, config)
            turbulent_dispersion!(Ur, alpha, ∇alpha, model.turbulence, Sc_t, config)
            interpolate_vanleer!(Urf, Ur, mdotf, config)
            zero_wall_drift_velocity!(Urf, config)
            face_dot_Sf!(Urdotf, Urf, config)
        end

        dt_sub = dt_cpu[1] / n_alpha_subcycles

        for _ in 1:n_alpha_subcycles
            @. alpha_prev.values = alpha.values

            grad!(∇alpha, alphaf, alpha, boundaries.alpha, time, config)
            limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha, config)

            if typeof(mp_model) <: VOF
                if n_smooth_compression > 0
                    laplacian_smooth!(alpha_smooth, alpha_smooth_tmp, alpha,
                                      n_smooth_compression, config)
                    grad!(∇alpha_smooth, alphaf, alpha_smooth, boundaries.alpha, time, config)
                    limit_gradient!(schemes.alpha.limiter, ∇alpha_smooth, alpha_smooth, config)
                    interpolate!(∇alphaf, ∇alpha_smooth.result, config)
                else
                    interpolate!(∇alphaf, ∇alpha.result, config)
                end
                compression_flux!(phirf, ∇alphaf, mdotf, C_alpha, config)
            end

            interpolate_upwind!(alphaf_upwind, alpha, mdotf, config)
            correct_boundaries!(alphaf_upwind, alpha, boundaries.alpha, time, config)
            interpolate_vanleer!(alphaf_HO, alpha, ∇alpha, mdotf, config)
            correct_boundaries!(alphaf_HO, alpha, boundaries.alpha, time, config)

            @. phiLf.values = mdotf.values * alphaf_upwind.values

            if typeof(mp_model) <: Mixture
                @. phiHf.values = mdotf.values * alphaf_HO.values -
                                  Urdotf.values * alphaf_upwind.values *
                                  (1.0 - alphaf_upwind.values)
            else
                @. phiHf.values = mdotf.values * alphaf_HO.values +
                                  phirf.values * alphaf_HO.values *
                                  (1.0 - alphaf_HO.values)
            end
            @. phiAf.values = phiHf.values - phiLf.values

            zero_boundary_faces!(phiAf, config)

            mules_limit!(mp_model,
                         phiAf, alpha_prev, phiLf,
                         Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
                         alphaMaxLocal, alphaMinLocal,
                         dt_sub, mesh, config)

            @. alpha_fluxf.values = phiLf.values + phiAf.values
            div!(div_alpha, alpha_fluxf, config)
            div!(div_mdotf, mdotf, config)
            @. alpha.values = alpha_prev.values -
                dt_sub * (div_alpha.values - alpha_prev.values * div_mdotf.values)

            interpolate_vanleer!(alphaf, alpha, ∇alpha, mdotf, config)
            correct_boundaries!(alphaf, alpha, boundaries.alpha, time, config)
        end

        ralpha = zero(TF)

        if field_mode
            # Phase-2B: refresh per-cell phase rho/mu from the current
            # T field (via the HelmholtzTable) and interpolate to faces
            # so the property blending below sees up-to-date values.
            live_update_phase_properties!(model)
            interpolate!(rho1f, phases[main].rho,      config)
            interpolate!(rho2f, phases[secondary].rho, config)
            interpolate!(mu1f,  phases[main].mu,       config)
            interpolate!(mu2f,  phases[secondary].mu,  config)

            blend_properties!(rho,  alpha,  phases[main].rho, phases[secondary].rho)
            blend_properties!(rhof, alphaf, rho1f,            rho2f)
            blend_mixture_nu!(nu,  alpha,  rho,  phases[main].mu, phases[secondary].mu)
            blend_mixture_nu!(nuf, alphaf, rhof, mu1f,            mu2f)

            if typeof(mp_model) <: Mixture
                # τ_d is per-cell now (depends on rho_secondary and mu_main).
                @. tau_d_field.values = (phases[secondary].rho.values * diameter^2) /
                                        (18.0 * phases[main].mu.values + eps())
            end
        else
            blend_properties!(rho,  alpha,  rho1_val, rho2_val)
            blend_properties!(rhof, alphaf, rho1_val, rho2_val)
            blend_mixture_nu!(nu,  alpha,  rho,  mu1_val, mu2_val)
            blend_mixture_nu!(nuf, alphaf, rhof, mu1_val, mu2_val)
        end
        @. mueff.values = rhof.values * nueff.values

        if typeof(mp_model) <: Mixture
            compute_tensor_term!(alphaf, rhof, rho1f, rho2f, Urf, slip_momentum_term, config)
            div_tensor!(div_slip_momentum, slip_momentum_term, config)
        end

        if typeof(mp_model) <: Mixture
            @. rhoPhi.values = mdotf.values * rhof.values
        elseif field_mode
            # VOF rhoPhi with per-face per-phase densities
            @. rhoPhi.values = alpha_fluxf.values * (rho1f.values - rho2f.values) +
                               mdotf.values * rho2f.values
        else
            @. rhoPhi.values = alpha_fluxf.values * (rho1_val - rho2_val) +
                               mdotf.values * rho2_val
        end

        # κ refresh — must happen BEFORE the U predictor so the well-balanced
        # pressure source uses the current step's κ. Mirrors interFoam's
        # `mixture.correct()` between alphaEqn and UEqn in PIMPLE.
        if typeof(mp_model) <: VOF
            laplacian_smooth!(alpha_smooth, alpha_smooth_tmp, alpha, n_smooth, config)
            grad!(∇alpha, alphaf, alpha_smooth, boundaries.alpha, time, config)
            limit_gradient!(schemes.alpha.limiter, ∇alpha, alpha_smooth, config)
            interpolate!(∇alphaf, ∇alpha.result, config)
            nhat_prep!(nhatf_prep, alpha_smooth, ∇alphaf, config)
            div!(kappa, nhatf_prep, config)
            cell_grad_magnitude!(grad_alpha_mag, ∇alpha, config)
            interpolate_weighted!(kappaf, kappa, grad_alpha_mag, config)
        end

        # Override the predictor pressure-source with a well-balanced
        # reconstruction that includes gravity body force (and surface
        # tension for VOF) at the face level. Without this, the predictor
        # produces spurious flow at variable-ρ interfaces because its source
        # has only `-∇p_rgh` (cell-Gauss) with no body-force counter-term.
        if typeof(mp_model) <: VOF
            well_balanced_pressure_grad!(
                ∇p_rgh.result, pressure_force_face,
                p_rgh, rho, ghf, mesh, config;
                sigma=sigma, kappaf=kappaf, alpha=alpha)
        else
            well_balanced_pressure_grad!(
                ∇p_rgh.result, pressure_force_face,
                p_rgh, rho, ghf, mesh, config)
        end

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

            if typeof(mp_model) <: VOF
                if use_csf_smoothed
                    csf_flux_smoothed!(phi_gf, rDf, sigma, kappaf, alpha_smooth, config)
                else
                    surface_tension_flux!(rDf, sigma, kappaf, alpha, phi_gf, config)
                end
            end

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

            # if typeof(mp_model) <: VOF
            #     interpolate_vanleer!(alphaf, alpha, ∇alpha, mdotf, config)
            #     correct_boundaries!(alphaf, alpha, boundaries.alpha, time, config)
            #     blend_properties!(rhof, alphaf, rho1_val, rho2_val)
            #     @. rhoPhi.values = mdotf.values *
            #                        (alphaf.values * (rho1_val - rho2_val) + rho2_val)
            # end

            pressure_grad!(p_rgh, ∇p_rghf_deconstructed, phi_gf, rDf, config)
            reconstruct_operation!(∇p_rghf_reconstructed, ∇p_rghf_deconstructed, config)
            correct_velocity_rgh!(U, Hv, ∇p_rghf_reconstructed, rD, phi_g, config)
        end

        @. p.values = p_rgh.values + (rho.values * gh.values)

        turbulence!(turbulenceModel, model, S, prev, time, config)
        update_nueff!(nueff, nuf, model.turbulence, config)
        @. mueff.values = rhof.values * nueff.values

        # Note: live property refresh now happens at the top of the
        # outer iteration (inside the `field_mode` branch above) so the
        # U-equation, the rhoPhi flux and the energy equation all see
        # the same T-driven properties this step.

        # Solve the mixture-temperature equation (no-op when the energy model
        # is Isothermal). Uses the post-pressure-correction mdotf and the
        # blended ρ at cells, with cp/k blended inside `energy!`.
        step_multiphase_energy!(energyModel, model, prev, mdotf, rho, time, config)

        courant           = max_courant_number!(cellsCourant, model, config)
        alphaCourantOuter = max_alpha_courant_number!(cellsAlphaCourant, alpha, mdotf, model, config, dt_cpu[1])
        # The MULES stability criterion applies to the per-sub-step α-CFL, not
        # the outer one. With N sub-cycles the per-step value is divided by N;
        # this is what update_dt! limits against maxAlphaCo and what we display.
        alphaCourant = alphaCourantOuter / n_alpha_subcycles
        update_dt!(config.runtime, courant, alphaCourant)

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

        if iteration % 100 == 0
            Ux_cpu = Array(U.x.values)
            Uy_cpu = Array(U.y.values)
            Uz_cpu = Array(U.z.values)
            Umag_cpu = sqrt.(Ux_cpu.^2 .+ Uy_cpu.^2 .+ Uz_cpu.^2)
            mean_Umag = sum(Umag_cpu) / length(Umag_cpu)
            max_Umag  = maximum(Umag_cpu)
            open(mean_U_csv_path, "a") do io
                println(io, "$iteration,$time,$mean_Umag,$max_Umag")
            end
        end

        runtime_postprocessing!(postprocess, iteration, iterations)

        if iteration % write_interval + signbit(write_interval) == 0
            save_output(model, outputWriter, iteration, time, config)
            save_postprocessing(postprocess, iteration, time, mesh, outputWriter, config.boundaries)
        end
    end

    # End-of-run velocity magnitude summary.
    Ux_final = Array(U.x.values)
    Uy_final = Array(U.y.values)
    Uz_final = Array(U.z.values)
    Umag_final = sqrt.(Ux_final.^2 .+ Uy_final.^2 .+ Uz_final.^2)
    mean_Umag_final = sum(Umag_final) / length(Umag_final)
    max_Umag_final  = maximum(Umag_final)
    @info "Final velocity summary" mean_Umag=mean_Umag_final max_Umag=max_Umag_final

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

function compute_Ur!(Ur, alpha, rho, g, DUmDt, rho1, rho2, mu1, d, tau_d, config)
    (; hardware) = config
    (; backend, workgroup) = hardware

    ndrange = length(Ur)
    kernel! = _compute_Ur!(_setup(backend, workgroup, ndrange)...)
    kernel!(Ur, alpha, rho, g, DUmDt, rho1, rho2, mu1, d, tau_d)
end

@kernel inbounds=true function _compute_Ur!(Ur, alpha, rho, g, DUmDt, rho1, rho2, mu1, d, tau_d)
    i = @index(Global)
    TF = eltype(rho.values)

    # Per-cell phase properties — these arguments are passed as
    # `ConstantScalar` (uniform constants) or `ScalarField` (per-cell,
    # populated each iteration by the HelmholtzTable live update).
    # Both implement `[i]`, returning the local value.
    rho_m = rho[i]
    rho_c = rho1[i]
    rho_d = rho2[i]
    mu_c  = mu1[i]
    tau   = tau_d[i]

    Ur_mag = norm(Ur[i])
    Re_p   = rho_c * Ur_mag * d / (mu_c + eps(TF))

    f_drag = ifelse(
        Re_p < TF(1000),
        one(TF) + TF(0.15) * Re_p^TF(0.687),
        TF(0.0183) * Re_p
    )
    f_drag = max(f_drag, one(TF))

    a_eff    = g - DUmDt[i]
    buoyancy = (rho_d - rho_m) / (rho_d + eps(TF))

    U_dm = (tau / (f_drag + eps(TF))) * buoyancy * a_eff

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


function update_fixedfluxpressure_gradient!(p_BCs, phi_gf, rDf, mesh)
    phi_vals = phi_gf.values
    rDf_vals = rDf.values
    faces_cpu = mesh.faces
    for BC in p_BCs
        if typeof(BC) <: FixedFluxPressure
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

"""
    blend_properties!(property_field, alpha_field, property_0, property_1)

Linearly blend a per-phase scalar property using the volume-fraction field:

    property_field = property_0 * alpha + property_1 * (1 - alpha)

`property_0` is the primary (tracked) phase value (`alpha = 1`), `property_1`
is the secondary phase value (`alpha = 0`).
"""
function blend_properties!(property_field, alpha_field, property_0, property_1)
    @. property_field.values = (property_0 * alpha_field.values) + (property_1 * (1.0 - alpha_field.values))
    nothing
end

# Field-aware variant — used when the per-phase property is itself a
# spatially-varying field (Phase-2B live mode driven by a HelmholtzTable).
# Both `property_0` and `property_1` carry per-cell values that get
# blended cell-by-cell with `alpha_field`.
function blend_properties!(property_field, alpha_field,
                           property_0::AbstractScalarField,
                           property_1::AbstractScalarField)
    @. property_field.values = (property_0.values * alpha_field.values) +
                               (property_1.values * (1.0 - alpha_field.values))
    nothing
end

"""
    blend_mixture_nu!(nu_field, alpha_field, rho_field, mu_0, mu_1)

Mixture kinematic viscosity built from a *linear dynamic-viscosity blend*:

    nu = (mu_0 * alpha + mu_1 * (1 - alpha)) / rho_blend

Using this for `nu`/`nuf` makes the downstream `mueff = rhof * nueff` equal to
the interFOAM-style linear μ-blend exactly, instead of the cross-product
`rhof * (ν₀·α + ν₁·(1−α))` which overshoots μ at the interface whenever
`ν₀ ≠ ν₁`.
"""
function blend_mixture_nu!(nu_field, alpha_field, rho_field, mu_0, mu_1)
    @. nu_field.values = (mu_0 * alpha_field.values + mu_1 * (1.0 - alpha_field.values)) / rho_field.values
    nothing
end

# Field-aware variant for Phase-2B: per-cell μ_0, μ_1 fields.
function blend_mixture_nu!(nu_field, alpha_field, rho_field,
                           mu_0::AbstractScalarField,
                           mu_1::AbstractScalarField)
    @. nu_field.values = (mu_0.values * alpha_field.values +
                          mu_1.values * (1.0 - alpha_field.values)) / rho_field.values
    nothing
end

"""
    compute_gh!(gh, g, config)

Computes `g · x` at cell centres. Used for hydrostatic pressure reconstruction.
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

Computes `g · x` at face centres.
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

    if typeof(mesh) <: Mesh2
        kernel! = _reconstruct_operation_2D!(_setup(backend, workgroup, ndrange)...)
    else
        kernel! = _reconstruct_operation_3D!(_setup(backend, workgroup, ndrange)...)
    end
    kernel!(cells, F, cell_faces, cell_nsign, faces, phi, psif)
end

@kernel function _reconstruct_operation_2D!(
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

@kernel function _reconstruct_operation_3D!(
    cells::AbstractArray{Cell{TF,SV,UR}}, F, cell_faces, cell_nsign, faces, phi, psif
) where {TF,SV,UR}
    i = @index(Global)
    @inbounds begin
        (; faces_range) = cells[i]

        m11 = zero(TF); m12 = zero(TF); m13 = zero(TF)
                        m22 = zero(TF); m23 = zero(TF)
                                        m33 = zero(TF)
        b1 = zero(TF);  b2 = zero(TF);  b3 = zero(TF)

        for fi ∈ faces_range
            fID = cell_faces[fi]
            (; area, normal) = faces[fID]
            nx = normal[1]; ny = normal[2]; nz = normal[3]

            m11 += area * nx * nx
            m12 += area * nx * ny
            m13 += area * nx * nz
            m22 += area * ny * ny
            m23 += area * ny * nz
            m33 += area * nz * nz

            ssf = psif[fID]
            b1 += nx * ssf
            b2 += ny * ssf
            b3 += nz * ssf
        end

        A11 = m22*m33 - m23*m23
        A12 = m13*m23 - m12*m33
        A13 = m12*m23 - m13*m22
        det = m11*A11 + m12*A12 + m13*A13
        invdet = abs(det) > eps(TF) ? one(TF)/det : zero(TF)

        ux = (A11*b1 + A12*b2 + A13*b3) * invdet
        uy = (A12*b1 + (m11*m33 - m13*m13)*b2 + (m13*m12 - m11*m23)*b3) * invdet
        uz = (A13*b1 + (m13*m12 - m11*m23)*b2 + (m11*m22 - m12*m12)*b3) * invdet

        phi[i] = @SVector [ux, uy, uz]
    end
end


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

"""
    mules_limit!(mp_model::AbstractMultiphaseModel,
                 phiAf, alpha_prev, phiLf,
                 Pplus, Pminus, Qplus, Qminus, Rplus, Rminus,
                 alphaMaxLocal, alphaMinLocal,
                 dt, mesh, config)

Unified MULES (Zalesak FCT) flux-limiter for explicit α transport, shared by
both `VOF` and `Mixture` multiphase sub-models. Limits the anti-diffusive
face flux `phiAf = phiHf - phiLf` so the explicit update
`α^{n+1} = α^n - dt/V · div(phiLf + phiAf)` stays bounded.

Pipeline:

  1. `mules_set_bounds!(mp_model, ...)` — sets per-cell `αMax`, `αMin`.
      Dispatched:
        - `VOF`     → local neighbour-stencil extrema of `α_prev` ∩ [0,1];
                      prevents new local extrema, preserves sharp interfaces.
        - `Mixture` → hard `[0,1]`; suspensions must tolerate uniform α
                      without freezing the anti-diffusive drift flux.
  2. Cell accumulation — build low-order α*, P±, Q± (shared).
  3. Ratios — R± = clamp(Q±/P±, 0, 1) (shared).
  4. Per-face apply — λ_f scales phiAf in place (shared).

Boundary faces are left unlimited (λ = 1): they have a single owner and
don't threaten boundedness for scalar α under incompressible transport.
"""
function mules_limit!(mp_model::AbstractMultiphaseModel,
                      phiAf, alpha_prev, phiLf,
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

    mules_set_bounds!(mp_model, alphaMaxLocal, alphaMinLocal, alpha_prev, mesh, config)

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

"""
    mules_set_bounds!(::Mixture, aMax, aMin, alpha_prev, mesh, config)

Hard physical bounds `[0, 1]` for a drift-flux mixture. The VOF-style
local-extrema bound is *wrong* here: a uniformly suspended initial state
(α = α₀ everywhere) would pin αMax = αMin = α₀ and the limiter would kill
every anti-diffusive flux, including the physical drift that drives
settling. For a mixture of two miscible-in-volume-fraction phases, `[0,1]`
is the only defensible bound.
"""
function mules_set_bounds!(::Mixture, alphaMaxLocal, alphaMinLocal, alpha_prev, mesh, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(mesh.cells)
    kernel! = _mmp_mules_hard_bounds!(_setup(backend, workgroup, ndrange)...)
    kernel!(alphaMaxLocal, alphaMinLocal)
end

"""
    mules_set_bounds!(::VOF, aMax, aMin, alpha_prev, mesh, config)

Local neighbour-stencil bounds for sharp-interface VOF:
`αMax_i = min(1, max(α_prev over i + its face neighbours))`,
`αMin_i = max(0, min(...))`. Prevents new local extrema so bulk cells with
`α_prev ≈ 0` cannot be pushed up by anti-diffusive noise — no phantom
droplets in the wake, no spurious peaks.
"""
function mules_set_bounds!(::VOF, alphaMaxLocal, alphaMinLocal, alpha_prev, mesh, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    (; cells, cell_faces, faces) = mesh
    ndrange = length(cells)
    kernel! = _mmp_mules_stencil_bounds!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, cell_faces, faces, alpha_prev, alphaMaxLocal, alphaMinLocal)
end

@kernel inbounds=true function _mmp_mules_hard_bounds!(alphaMaxLocal, alphaMinLocal)
    i = @index(Global)
    TF = eltype(alphaMaxLocal.values)
    alphaMaxLocal[i] = one(TF)
    alphaMinLocal[i] = zero(TF)
end

@kernel inbounds=true function _mmp_mules_stencil_bounds!(
    cells::AbstractArray{Cell{TF,SV,UR}}, cell_faces, faces,
    alpha_prev, alphaMaxLocal, alphaMinLocal
) where {TF,SV,UR}
    i = @index(Global)
    (; faces_range) = cells[i]
    aMax = alpha_prev[i]
    aMin = alpha_prev[i]
    for fi in faces_range
        fID = cell_faces[fi]
        oc = faces[fID].ownerCells
        j = ifelse(oc[1] == i, oc[2], oc[1])
        aj = alpha_prev[j]
        aMax = max(aMax, aj)
        aMin = min(aMin, aj)
    end
    alphaMaxLocal[i] = min(one(TF),  aMax)
    alphaMinLocal[i] = max(zero(TF), aMin)
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


"""
    compression_flux!(phirf, ∇alphaf, mdotf, C_alpha, config)

interFOAM-style anti-diffusive compression face flux:

    phir_f · Sf = min(Cα · |φf|/|Sf|, Φmax) · (n̂f · Sf)

with `n̂f = (∇α)f / |(∇α)f|` and `Φmax = max_f |φf|/|Sf|`. The dot with `Sf`
is baked in, so the result drops straight into the α flux as
`phirf · α(1−α)`.
"""
function compression_flux!(phirf, ∇alphaf, mdotf, C_alpha, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh  = phirf.mesh
    faces = mesh.faces
    TF = _get_float(mesh)

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
    (; area, normal, delta) = face
    TF = eltype(phirf.values)

    Sf = area * normal
    g  = ∇alphaf[i]
    mag = norm(g)

    deltaN = TF(1e-8) / delta
    nhat = mag > deltaN ? g / mag : zero(g)

    phi_over_S = abs(mdotf[i]) / (area + eps(TF))
    compr_speed = min(C_alpha * phi_over_S, phimax)

    phirf[i] = compr_speed * (nhat ⋅ Sf)
end

"""
    laplacian_smooth!(phi_smooth, tmp, phi, n_passes, config)

Face-area-weighted Laplacian smoothing of `phi` into `phi_smooth`, using
`tmp` as a ping-pong buffer. Used to dampen high-frequency α noise before
computing curvature κ for CSF.
"""
function laplacian_smooth!(phi_smooth, tmp, phi, n_passes, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = phi.mesh
    (; cells, cell_faces, cell_nsign, faces) = mesh
    ndrange = length(cells)

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
            oc = faces[fID].ownerCells
            j  = ifelse(oc[1] == i, oc[2], oc[1])
            sum_wa += area * phi.values[j]
            sum_w  += area
        end
        out.values[i] = sum_wa / sum_w
    end
end

"""
    cell_grad_magnitude!(mag_field, grad, config)

Cell-centred `|∇α|` from a gradient field, used as the |∇α|-weighted face
interpolation weight for κf so wall-adjacent cells with small |∇α| don't
contaminate the surface-tension force.
"""
function cell_grad_magnitude!(mag_field, grad, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = mag_field.mesh
    ndrange = length(mesh.cells)
    kernel! = _cell_grad_magnitude!(_setup(backend, workgroup, ndrange)...)
    kernel!(mag_field, grad.result)
end

@kernel inbounds=true function _cell_grad_magnitude!(mag_field, grad_result)
    i = @index(Global)
    mag_field.values[i] = norm(grad_result[i])
end

"""
    interpolate_weighted!(phif, phi, weight_field, config)

Weighted face interpolation: `phif[f] = (φ_c1·w_c1 + φ_c2·w_c2) /
(w_c1 + w_c2 + ε)`. With `w = |∇α|` it preserves κf at interfaces and zeros
it in bulk where κ is noisy.
"""
function interpolate_weighted!(phif::FaceScalarField, phi::ScalarField,
                                 weight_field::ScalarField, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = phif.mesh
    (; faces) = mesh
    ndrange = length(faces)
    kernel! = _interpolate_weighted!(_setup(backend, workgroup, ndrange)...)
    kernel!(phif, phi, weight_field, faces)
end

@kernel inbounds=true function _interpolate_weighted!(phif, phi, w, faces)
    i = @index(Global)
    face = faces[i]
    (; ownerCells) = face
    c1 = ownerCells[1]
    c2 = ownerCells[2]
    TF = eltype(phif.values)
    w1 = w.values[c1]
    w2 = w.values[c2]
    denom = w1 + w2 + eps(TF)
    phif[i] = (phi.values[c1] * w1 + phi.values[c2] * w2) / denom
end

"""
    nhat_prep!(nhatf_prep, alpha, ∇alphaf, config)

Builds a face-normal unit vector field `n̂f` from the face-interpolated α
gradient. Boundary faces are zeroed (avoids biased div(n̂) injection); on
internal faces, `|∇α|` below `1e-8/delta` is treated as numerical noise and
the corresponding `n̂` is set to zero.
"""
function nhat_prep!(nhatf_prep, alpha, ∇alphaf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    mesh = alpha.mesh
    faces = mesh.faces
    nbfaces = length(mesh.boundary_cellsID)
    nfaces  = length(faces)

    if nbfaces > 0
        kernel! = _nhat_zero_bfaces!(_setup(backend, workgroup, nbfaces)...)
        kernel!(nhatf_prep)
    end

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
    (; delta) = face

    g = ∇alphaf_[fID]
    mag = norm(g)
    deltaN = 1e-8 / delta

    nhatf_prep[fID] = mag < deltaN ? SVector(0.0, 0.0, 0.0) : g / mag
end

"""
    surface_tension_flux!(rDf, sigma, kappaf, alpha, phi_gf, config)

CSF surface-tension contribution to the face flux:
`phi_gf -= σ · κf · (∇αf · Sf) · rDf`.
"""
function surface_tension_flux!(rDf, sigma, kappaf, alpha, phi_gf, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    faces = phi_gf.mesh.faces

    ndrange = length(phi_gf)
    kernel! = _surface_tension_flux!(_setup(backend, workgroup, ndrange)...)
    kernel!(rDf, sigma, kappaf, alpha, phi_gf, faces)
end

@kernel inbounds=true function _surface_tension_flux!(rDf, sigma, kappaf, alpha, phi_gf, faces)
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

"""
    csf_flux_smoothed!(phi_gf, rDf, sigma, kappaf, alpha_smooth, config)

Experimental CSF that uses `alpha_smooth` for the face α gradient too.
"""
function csf_flux_smoothed!(phi_gf, rDf, sigma, kappaf, alpha_smooth, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    faces = phi_gf.mesh.faces
    ndrange = length(phi_gf)
    kernel! = _csf_flux_smoothed!(_setup(backend, workgroup, ndrange)...)
    kernel!(phi_gf, rDf, sigma, kappaf, alpha_smooth, faces)
end

@kernel inbounds=true function _csf_flux_smoothed!(phi_gf, rDf, sigma, kappaf, alpha_smooth, faces)
    i = @index(Global)
    face = faces[i]
    (; area, normal, ownerCells, delta) = face
    Sf = area * normal
    c1 = ownerCells[1]
    c2 = ownerCells[2]
    as1 = alpha_smooth[c1]
    as2 = alpha_smooth[c2]
    ∇αf_vec = normal * ((as2 - as1) / delta)
    phi_gf[i] -= sigma * kappaf[i] * (∇αf_vec ⋅ Sf) * rDf[i]
end


"""
    well_balanced_pressure_grad!(grad_field, face_buf, p_rgh, rho, ghf,
                                  mesh, config; sigma=0, kappaf=nothing,
                                  alpha=nothing)

Overrides `grad_field.values` with a *combined* face-snGrad reconstruction
of the predictor body-force terms:

    face_buf[f] = area_f · ( snGrad(p_rgh) + ghf · snGrad(ρ) - σ·κf · snGrad(α) )

Equivalent to interFoam's
`fvc::reconstruct((surfaceTensionForce - ghf·snGrad(rho) - snGrad(p_rgh))·magSf)`,
flipped in sign to match the `-Source(grad_field)` convention in U_eqn.

For a discrete hydrostatic state, `snGrad(p_rgh) = -ghf·snGrad(ρ)` (and σ
either 0 or balanced by curvature) at every face, so the face value vanishes
and the reconstructed cell vector is zero — the predictor produces no
spurious U. This is the "well-balanced" property that interFoam's predictor
has and XCALibre's plain cell-Gauss `grad!(∇p_rgh, ...)` lacks.

When `kappaf === nothing` or `sigma == 0`, the σ term is omitted (used for
Mixture and σ=0 VOF cases). For surface-tension-active VOF cases pass
`sigma`, `kappaf`, and the cell `alpha` field. `kappaf` may be one step
stale (computed at the end of the previous outer iteration); the lag is
fine because curvature changes slowly relative to dt.
"""
function well_balanced_pressure_grad!(
    grad_field, face_buf, p_rgh, rho, ghf, mesh, config;
    sigma=zero(eltype(p_rgh.values)),
    kappaf=nothing, alpha=nothing,
)
    (; hardware) = config
    (; backend, workgroup) = hardware
    faces = mesh.faces

    ndrange = length(faces)
    if kappaf === nothing || alpha === nothing || iszero(sigma)
        kernel! = _well_balanced_pressure_face!(_setup(backend, workgroup, ndrange)...)
        kernel!(face_buf, p_rgh, rho, ghf, faces)
    else
        kernel! = _well_balanced_pressure_face_csf!(_setup(backend, workgroup, ndrange)...)
        kernel!(face_buf, p_rgh, rho, alpha, ghf, kappaf, sigma, faces)
    end

    reconstruct_operation!(grad_field, face_buf, config)
end

@kernel inbounds=true function _well_balanced_pressure_face!(
    face_buf, p_rgh, rho, ghf, faces
)
    i = @index(Global)
    face = faces[i]
    (; area, ownerCells, delta) = face
    c1 = ownerCells[1]
    c2 = ownerCells[2]

    snGrad_p   = (p_rgh[c2] - p_rgh[c1]) / delta
    snGrad_rho = (rho[c2]   - rho[c1])   / delta

    face_buf[i] = area * (snGrad_p + ghf[i] * snGrad_rho)
end

@kernel inbounds=true function _well_balanced_pressure_face_csf!(
    face_buf, p_rgh, rho, alpha, ghf, kappaf, sigma, faces
)
    i = @index(Global)
    face = faces[i]
    (; area, ownerCells, delta) = face
    c1 = ownerCells[1]
    c2 = ownerCells[2]

    snGrad_p   = (p_rgh[c2]   - p_rgh[c1])   / delta
    snGrad_rho = (rho[c2]     - rho[c1])     / delta
    snGrad_a   = (alpha[c2]   - alpha[c1])   / delta

    face_buf[i] = area * (snGrad_p + ghf[i] * snGrad_rho - sigma * kappaf[i] * snGrad_a)
end

