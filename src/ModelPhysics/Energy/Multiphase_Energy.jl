export MultiphaseEnergy

struct MultiphaseEnergy{S,SF,V,F} <: AbstractEnergyModel where {
    S<:ScalarField,
    SF<:FaceScalarField,
    V<:AbstractVectorField,
    F<:AbstractFloat
}
    T::S
    cp_prev::S
    p_prev::S
    beta_mix::S
    cp_mix::S
    cpf_mix::SF
    betaDpDt::S
    DcpDt::S
    time_term::S
    div_term::SF
    vf_p::V
    vf_q::V
    k_eff::SF
    Pr_t::F
end
Adapt.@adapt_structure MultiphaseEnergy


struct MultiphaseEnergyModel{E1,State}
    energy_eqn::E1
    state::State
end
Adapt.@adapt_structure MultiphaseEnergyModel



# Model API constructor
Energy{MultiphaseEnergy}(; Pr_t) = begin
    coeffs = (; Pr_t)
    ARG = typeof(coeffs)
    Energy{MultiphaseEnergy,ARG}(coeffs)
end

# Functor as constructor
(energy::Energy{EnergyModel, ARG})(mesh, fluid) where {EnergyModel<:MultiphaseEnergy,ARG} = begin

    coeffs = energy.args
    (; Pr_t) = coeffs

    T = ScalarField(mesh)
    cp_prev = ScalarField(mesh)
    p_prev = ScalarField(mesh)
    beta_mix = ScalarField(mesh)
    cp_mix = ScalarField(mesh)
    cpf_mix = FaceScalarField(mesh)
    
    betaDpDt = ScalarField(mesh)
    DcpDt = ScalarField(mesh)

    
    time_term = ScalarField(mesh)
    k_eff = FaceScalarField(mesh)
    div_term = FaceScalarField(mesh)

    vf_p = FaceVectorField(mesh)
    vf_q = FaceVectorField(mesh)
    
    MultiphaseEnergy(T, cp_prev, p_prev, beta_mix, cp_mix, cpf_mix, betaDpDt, DcpDt, time_term, div_term, vf_p, vf_q, k_eff, Pr_t)
end


function initialise(
    energy::MultiphaseEnergy, model::Physics{T1,F,SO,M,Tu,E,D,BI}, nutf, ∇p, config, mesh, dt
    ) where {T1,F,SO,M,Tu,E,D,BI}

    (; solvers, schemes, runtime, hardware, boundaries) = config

    backend = config.hardware.backend
    workgroup = config.hardware.workgroup

    compute_energy_properties!(model, mesh, config)


    @. cp_prev.values = cp_mix.values
    @. p_prev.values = p.values


    energy_rhs = - Src(betaDpDt, 1) # If we make this a pre created field and update it then it should be just fine?
    energy_rhs = construct_RHS(model.energy, energy_rhs, model, props, alpha, rho, phases, config, mesh)


    energy_eqn = (
        Time{schemes.T.time}(time_term, T)
        + Divergence{schemes.T.divergence}(div_term, T) #this is not mdotf times cp!
        - Laplacian{schemes.T.laplacian}(k_eff, T)  #cell centre of faces???
        - Si(DcpDt, T)
        == energy_rhs
    ) → ScalarEquation(T, boundaries.T)


    @reset energy_eqn.preconditioner = set_preconditioner(solvers.T.preconditioner, energy_eqn)
    
    # preallocating solvers
    @reset energy_eqn.solver = _workspace(solvers.T.solver, _b(energy_eqn))


    init_residual = (:T, 1.0) # ?
    init_converged = false # ?
    state = ModelState(init_residual, init_converged) # ?

    return MultiphaseEnergyModel(energy_eqn, state)
end

# Cannot simplify div term towards mixture values - maths is not working...
# Computing from face values would imply another HelmholtzE call.....

function energy!(
    energy::MultiphaseEnergyModel, model::Physics{T1,F,SO,M,Tu,E,D,BI}, energy_eqn, nutf, ∇p, config, mesh, dt, time
    ) where {T1,F,SO,M,Tu,E,D,BI}

    mesh = model.domain

    (; state) = energy

    compute_energy_properties!(model, mesh, config)

    update_sources!(model.energy, model, model.fluid.physics_properties, energy_eqn, alpha, rho, model.fluid.phases, config, mesh)

    # Set up and solve energy equation
    
    @. cp_prev.values = cp_mix.values
    @. p_prev.values = p.values
    
    T_res = solve_equation!(energy_eqn, T, boundaries.T, solvers.T, config; time=time)

    residuals = (:T, T_res)
    converged = T_res <= solvers.T.convergence
    state.residuals = residuals
    state.converged = converged


    return nothing
end




@kernel inbounds=true function _compute_betaDpDt!(betaDpDt, beta_mix, p, p_prev, dt, U, ∇p)
    i = @index(Global)

    betaDpDt[i] = beta_mix[i] * ( ( (p[i]-p_prev[i])/dt ) + (U[i] ⋅ ∇p[i]) )
end


@kernel inbounds=true function _compute_DcpDt!(DcpDt, rho, cp_mix, cp_prev, dt, U, ∇cp)
    i = @index(Global)

    DcpDt[i] = rho[i] * ( ( (cp_mix[i]-cp_prev[i])/dt ) + (U[i] ⋅ ∇cp[i]) )
end






function compute_energy_properties!(model, mesh, config)

    (; solvers, schemes, runtime, hardware, boundaries) = config

    backend = config.hardware.backend
    workgroup = config.hardware.workgroup
    props = model.fluid.physics_properties
    phases = model.fluid.phases

    Pr_t = model.energy.Pr_t
    T = model.energy.T
    beta_mix = model.energy.beta_mix
    cp_mix = model.energy.cp_mix
    cpf_mix = model.energy.cpf_mix

    betaDpDt = model.energy.betaDpDt
    DcpDt = model.energy.DcpDt

    alpha = model.fluid.alpha
    alphaf = model.fluid.alphaf
    rho = model.fluid.rho
    rhof = model.fluid.rhof

    U = model.momentum.U
    Uf = model.momentum.Uf

    p = model.momentum.p
    
    rho_p = model.fluid.phases[1].rho
    rho_q = model.fluid.phases[2].rho
    
    rhof_p = FaceScalarField(mesh)
    rhof_q = FaceScalarField(mesh)

    interpolate!(rhof_p, rho_p, config)
    interpolate!(rhof_q, rho_q, config)

    cp_p = model.fluid.phases[1].cp
    cp_q = model.fluid.phases[2].cp
    
    cpf_p = FaceScalarField(mesh)
    cpf_q = FaceScalarField(mesh)

    interpolate!(cpf_p, cp_p, config)
    interpolate!(cpf_q, cp_q, config)

    beta_p = model.fluid.phases[1].cp
    beta_q = model.fluid.phases[2].cp

    blend_properties!(beta_mix, alpha, beta_p, beta_q)
    blend_properties!(cp_mix, alpha, cp_p, cp_q)
    
    interpolate!(cpf_mix, cp_mix, config)
    
    v_p = model.fluid.physics_properties.driftVelocity.v_p
    v_q = model.fluid.physics_properties.driftVelocity.v_q

    vf_p = model.energy.vf_p
    vf_q = model.energy.vf_q

    interpolate!(vf_p, v_p, config)
    interpolate!(vf_q, v_q, config)

    cp_prev = model.energy.cp_prev
    p_prev = model.energy.p_prev


    time_term = model.energy.time_term

    div_term = model.energy.div_term
    div_term_adition = FaceScalarField(mesh)

    @. time_term.values = ( (alpha.values * (rho_p.values*cp_p.values)) + ((1.0-alpha.values) * (rho_q.values*cp_q.values)) ) # * T

    @. div_term.values = (alphaf.values * (rhof_p.values*cpf_p.values)) # * T
    @. div_term_adition.values = ((1.0-alphaf.values) * (rhof_q.values*cpf_q.values)) # * T

    flux!(div_term, vf_p, config)
    flux!(div_term_adition, vf_q, config)

    @. div_term.values = div_term.values + div_term_adition.values


    #Calculate cpf from the face values (T, p)

    ∇cp = Grad{schemes.p.gradient}(cp)
    grad!(∇cp, cpf_mix, cp_mix, time, config)
    limit_gradient!(schemes.p.limiter, ∇cp, cp, config)

    
    ndrange = length(betaDpDt)
    kernel! = _compute_betaDpDt!(_setup(backend, workgroup, ndrange)...)
    kernel!(betaDpDt, beta_mix, p, p_prev, dt, U, ∇p)
    
    ndrange = length(DcpDt)
    kernel! = _compute_DcpDt!(_setup(backend, workgroup, ndrange)...)
    kernel!(DcpDt, rho, cp_mix, cp_prev, dt, U, ∇cp)

    k_mix = ScalarField(mesh)
    kf_mix = FaceScalarField(mesh)

    k_p = model.fluid.phases[1].k
    k_q = model.fluid.phases[2].k
    blend_properties!(k_mix, alpha, k_p, k_q)
    interpolate!(kf_mix, k_mix, config)

    k_eff = model.energy.k_eff
    @. k_eff.values = (rhof.values * nutf.values * cpf_mix.values) / Pr_t

    @. k_eff.values = kf_mix.values + k_eff.values # when laminar, nut is a scalar field of 0 thus k_t = 0
end