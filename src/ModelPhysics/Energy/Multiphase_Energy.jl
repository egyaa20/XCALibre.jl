export MultiphaseEnergy

struct MultiphaseEnergy{S1,S2,S3,S4,S5,F1} <: AbstractEnergyModel
    T::S1
    cp_prev::S2
    p_prev::S3
    beta_mix::S4
    cp_mix::S5
    cpf_mix::F1
end
Adapt.@adapt_structure MultiphaseEnergy


struct MultiphaseEnergyModel{E1,State}
    energy_eqn::E1 
    state::State
end
Adapt.@adapt_structure MultiphaseEnergyModel



# Model API constructor
Energy{MultiphaseEnergy}() = begin
    coeffs = (nothing)
    ARG = typeof(coeffs)
    Energy{MultiphaseEnergy,ARG}(coeffs)
end

# Functor as constructor
(energy::Energy{EnergyModel, ARG})(mesh, fluid) where {EnergyModel<:MultiphaseEnergy,ARG} = begin
    T = ScalarField(mesh)
    cp_prev = ScalarField(mesh)
    p_prev = ScalarField(mesh)
    beta_mix = ScalarField(mesh)
    cp_mix = ScalarField(mesh)
    cpf_mix = FaceScalarField(mesh)
    
    MultiphaseEnergy(h)
end


function initialise(
    energy::MultiphaseEnergy, model::Physics{T1,F,SO,M,Tu,E,D,BI}, nut, latentHeat, m_qp, m_pq, ∇p, config, mesh, dt #how to get Pr_t ?
    ) where {T1,F,SO,M,Tu,E,D,BI}

    Pr_t = 0.9 # Prob better to be defined by user via API - turbulent Pr number

    backend = config.hardware.backend
    workgroup = config.hardware.workgroup

    T = model.energy.T
    beta_mix = model.energy.beta_mix
    cp_mix = model.energy.cp_mix
    cpf_mix = model.energy.cpf_mix

    alpha = model.fluid.alpha
    alphaf = model.fluid.alphaf
    rho = model.fluid.rho

    U = model.momentum.U
    Uf = model.momentum.Uf
    
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

    cp_prev = model.energy.cp_prev
    p_prev = model.energy.p_prev


    time_term = ScalarField(mesh)

    div_term = FaceScalarField(mesh)
    div_term_adition = FaceScalarField(mesh)

    betaDpDt = ScalarField(mesh)
    DcpDt = ScalarField(mesh)

    @. time_term.values = ( (alpha.values * (rho_p.values*cp_p.values)) + ((1.0-alpha.values) * (rho_q.values*cp_q.values)) ) # * T

    @. div_term.values = (alphaf.values * (rhof_p.values*cpf_p.values)) # * T
    @. div_term_adition.values = ((1.0-alphaf.values) * (rhof_q.values*cpf_q.values)) # * T

    flux!(div_term, v_p, config)
    flux!(div_term_adition, v_q, config)

    @. div_term.values = div_term.values + div_term_adition.values



    ∇cp = Grad{schemes.p.gradient}(cp)
    grad!(∇cp, cpf_mix, cp_mix, boundaries.p, time, config) # WHAT TO PASS IN BCs ????
    limit_gradient!(schemes.p.limiter, ∇cp, cp, config)


    
    ndrange = length(betaDpDt)
    kernel! = _compute_betaDpDt!(_setup(backend, workgroup, ndrange)...)
    kernel!(betaDpDt, beta_mix, p, p_prev, dt, U, ∇p)
    
    ndrange = length(DcpDt)
    kernel! = _compute_DcpDt!(_setup(backend, workgroup, ndrange)...)
    kernel!(DcpDt, rho, cp_mix, cp_prev, dt, U, ∇cp)

    k_mix = ScalarField(mesh)
    k_p = model.fluid.phases[1].k
    k_q = model.fluid.phases[2].k
    blend_properties!(k_mix, alpha, k_p, k_q)

    k_eff = ScalarField(mesh)
    @. k_eff.values = (rho.values * nut.values * cp_mix.values) / Pr_t

    @. k_eff.values = k_mix.values + k_eff.values # comment out k_eff to remove turbulence contribution







    energy_eqn = (
        Time{schemes.T.time}(time_term, T)
        + Divergence{schemes.T.divergence}(div_term, T)
        - Laplacian{schemes.T.laplacian}(k_eff, T) # how does this term definition actually work? does it take divergence and grad on its own?
        - Si(DcpDt, T) # SEMI IMPLICIT
        == 
        - Source(betaDpDt)
        + Source(m_qp*latentHeat)
        - Source(m_pq*latentHeat)
    ) → eqn


    @reset energy_eqn.preconditioner = set_preconditioner(solvers.h.preconditioner, energy_eqn)
    
    # preallocating solvers
    @reset energy_eqn.solver = _workspace(solvers.h.solver, _b(energy_eqn))

    init_residual = (:h, 1.0)
    init_converged = false
    state = ModelState(init_residual, init_converged)

    return MultiphaseEnergyModel(energy_eqn, state) #remove state
end



function energy!(
    energy::MultiphaseEnergyModel, model::Physics{T1,F,SO,M,Tu,E,D,BI}, 
        alpha, p, rho_l, rho_v, u_l, u_v, U_m, U_m_prev, mdotf, mueff, time, config
    ) where {T1,F,SO,M,Tu,E,D,BI}

    # Fields required:
        # alpha
        # p
        # rho_k
        # u_k
        # U_m

        # prev U_m

    # Steps:
    # 1) ( U - prev_U ) / dt
    # dUdt = VectorField(mesh)
    # 2) grad (U)
    # 3) Compute a
    # 4) Compute Re
    # 5) Compute f_drag
    
    # Keep turbulence = 0 for now
    # 6) Compute rho_m
    # 7) Compute d_p
    # 8) Compute v_pq
    # 9) Compute v_dr
    # 10) Compute v_k

    # 11) Compute E_k
    # 12) Compute S_h

    # 13) Compute k_eff
    # 14) Compute tau_eff
    # 15) Compute laplacian

    # 16) Compute Lee Model field
    # 17) Compute RPI field

    # SOLVE ENERGY EQN.

    # TBD: MOMENTUM;    VOF;    Sketch algo;    PPT

    mesh = model.domain

    return nothing
end




@kernel inbounds=true function _compute_betaDpDt!(betaDpDt, beta_mix, p, p_prev, dt, U, ∇p)
    i = @index(Global)

    betaDpDt[i] = beta_mix[i] * ( ( (p[i]-p_prev[i])/dt ) + (U[i] * ∇p[i]) ) # is it the dot or not actually?
end


@kernel inbounds=true function _compute_DcpDt!(DcpDt, rho, cp_mix, cp_prev, dt, U, ∇cp)
    i = @index(Global)

    DcpDt[i] = rho[i] * ( ( (cp_mix[i]-cp_prev[i])/dt ) + (U[i] * ∇cp[i]) ) # is it the dot or not actually?
end



# export SensibleEnthalpy
# export Ttoh, htoT!, Ttoh!, thermo_Psi!

# # Model type definition
# """
#     SensibleEnthalpy <: AbstractEnergyModel

# Type that represents energy model, coefficients and respective fields.

# ### Fields
# - `h`: Sensible enthalpy ScalarField.
# - `T`: Terature ScalarField.
# - `hf`: Sensible enthalpy FaceScalarField.
# - `Tf`: Temperature FaceScalarField.
# - `K`: Specific kinetic energy ScalarField.
# - `dpdt`: Pressure time derivative ScalarField.
# - `coeffs`: A tuple of model coefficients.

# """
# # struct SensibleEnthalpy{S1,S2,F1,F2,S3,S4,F,C} <: AbstractEnergyModel
# struct SensibleEnthalpy{S1,S2,F1,F2,S3,S4,C} <: AbstractEnergyModel
#     h::S1
#     T::S2
#     hf::F1
#     Tf::F2
#     K::S3
#     dpdt::S4
#     # update_BC::F
#     coeffs::C
# end
# Adapt.@adapt_structure SensibleEnthalpy

# struct SensibleEnthalpyModel{E1,State}
#     energy_eqn::E1 
#     state::State
# end
# Adapt.@adapt_structure SensibleEnthalpyModel

# # Model API constructor
# Energy{SensibleEnthalpy}(; Tref) = begin
#     coeffs = (Tref=Tref, other=nothing)
#     ARG = typeof(coeffs)
#     Energy{SensibleEnthalpy,ARG}(coeffs)
# end

# # Functor as constructor
# (energy::Energy{EnergyModel, ARG})(mesh, fluid) where {EnergyModel<:SensibleEnthalpy,ARG} = begin
#     h = ScalarField(mesh)
#     T = ScalarField(mesh)
#     hf = FaceScalarField(mesh)
#     Tf = FaceScalarField(mesh)
#     K = ScalarField(mesh)
#     dpdt = ScalarField(mesh)
#     # update_BC =  return_thingy(EnergyModel, fluid, energy.args.Tref)
#     coeffs = energy.args
#     # SensibleEnthalpy(h, T, hf, Tf, K, dpdt, update_BC, coeffs)
#     SensibleEnthalpy(h, T, hf, Tf, K, dpdt, coeffs)
# end

# """
#     initialise(energy::SensibleEnthalpy, model::Physics{T1,F,SO,M,Tu,E,D,BI}, mdotf, rho, peqn, config
#     ) where {T1,F,SO,M,Tu,E,D,BI})

# Initialisation of energy transport equations.

# # Input
# - `energy`: Energy model.
# - `model`: Physics model defined by user.
# - `mdtof`: Face mass flow.
# - `rho`: Density ScalarField.
# - `peqn`: Pressure equation.
# - `config`: Configuration structure defined by user with solvers, schemes, runtime and 
#               hardware structures set.

# # Output
# - `SensibleEnthalpyModel`: Energy model struct containing energy equation.

# """
# function initialise(
#     energy::SensibleEnthalpy, model::Physics{T1,F,SO,M,Tu,E,D,BI}, mdotf, rho, peqn, config
#     ) where {T1,F,SO,M,Tu,E,D,BI}

#     (; h, T, dpdt) = energy
#     (; solvers, schemes, runtime, boundaries) = config
#     mesh = mdotf.mesh
#     eqn = peqn.equation
    
#     # rho = ScalarField(mesh)
#     keff_by_cp = FaceScalarField(mesh)
#     divK = ScalarField(mesh)
#     dKdt = ScalarField(mesh)

#     Ttoh!(model, T, h)

#     energy_eqn = (
#         Time{schemes.h.time}(rho, h)
#         + Divergence{schemes.h.divergence}(mdotf, h) 
#         - Laplacian{schemes.h.laplacian}(keff_by_cp, h) 
#         == 
#         -Source(divK)
#         -Source(dKdt)
#         +Source(dpdt)
#     ) → eqn
    
#     # Set up preconditioners
#     # @reset energy_eqn.preconditioner = set_preconditioner(
#     #             solvers.h.preconditioner, energy_eqn, boundaries.h, config)

#     @reset energy_eqn.preconditioner = set_preconditioner(solvers.h.preconditioner, energy_eqn)
    
#     # preallocating solvers
#     @reset energy_eqn.solver = _workspace(solvers.h.solver, _b(energy_eqn))

#     init_residual = (:h, 1.0)
#     init_converged = false
#     state = ModelState(init_residual, init_converged)

#     return SensibleEnthalpyModel(energy_eqn, state)
# end


# """
#     energy::SensibleEnthalpyModel, model::Physics{T1,F,SO,M,Tu,E,D,BI}, prev, mdotf, rho, mueff, time, config
#     ) where {T1,F,SO,M,Tu,E,D,BI,E1}

# Run energy transport equations.

# # Input
# - `energy`: Energy model.
# - `model`: Physics model defined by user.
# - `prev`: Previous energy cell values.
# - `mdtof`: Face mass flow.
# - `rho`: Density ScalarField.
# - `mueff`: Effective viscosity FaceScalarField.
# - `time`: Simulation runtime.
# - `config`: Configuration structure defined by user with solvers, schemes, runtime and hardware structures set.

# """
# function energy!(
#     energy::SensibleEnthalpyModel, model::Physics{T1,F,SO,M,Tu,E,D,BI}, prev, mdotf, rho, mueff, time, config
#     ) where {T1,F,SO,M,Tu,E,D,BI}

#     mesh = model.domain

#     (;U) = model.momentum
#     (;h, hf, T, K, dpdt) = model.energy
#     (;energy_eqn, state) = energy
#     (; solvers, runtime, hardware, boundaries) = config
#     (; iterations, write_interval) = runtime
#     (; backend) = hardware

#     # rho = get_flux(energy_eqn, 1)
#     keff_by_cp = get_flux(energy_eqn, 3)
#     divK = get_source(energy_eqn, 1)
#     dKdt = get_source(energy_eqn, 2)

#     Uf = FaceVectorField(mesh)
#     Kf = FaceScalarField(mesh)
#     # Kbounded = ScalarField(mesh)
#     Pr = model.fluid.Pr

#     dt = runtime.dt

#     # Pre-allocate auxiliary variables
#     TF = _get_float(mesh)
#     n_cells = length(mesh.cells)
#     # prev = zeros(TF, n_cells)
#     # prev = _convert_array!(prev, backend)
#     prev = KernelAbstractions.zeros(backend, TF, n_cells) 

#     volumes = getproperty.(mesh.cells, :volume)

#     @. keff_by_cp.values = mueff.values/Pr.values

#     @. prev = K.values
#     interpolate!(Uf, U, config)
#     correct_boundaries!(Uf, U, boundaries.U, time, config)
#     for i ∈ eachindex(K)
#         K.values[i] = 0.5*(U.x.values[i]^2 + U.y.values[i]^2 + U.z.values[i]^2)
#     end
#     interpolate!(Kf, K, config)
#     for i ∈ eachindex(Kf)
#         Kf.values[i] = 0.5*(Uf.x.values[i]^2 + Uf.y.values[i]^2 + Uf.z.values[i]^2)
#     end
#     # correct_face_interpolation!(Kf, K, mdotf) # This forces KE to be upwind, MIGHT NOT BE WORKING
#     @. Kf.values *= mdotf.values
#     div!(divK, Kf, config)

#     if config.schemes.h.time <: SteadyState
#         @. dKdt.values = 0.0
#     else
#         @. dKdt.values = rho.values*(K.values - prev)/dt
#     end

#     # Set up and solve energy equation
#     @. prev = h.values
#     discretise!(energy_eqn, h, config)
#     apply_boundary_conditions!(energy_eqn, boundaries.h, nothing, time, config)
#     implicit_relaxation_diagdom!(energy_eqn, h.values, solvers.h.relax, nothing, config)
#     update_preconditioner!(energy_eqn.preconditioner, mesh, config)
#     h_res = solve_system!(energy_eqn, solvers.h, h, nothing, config)

#     if !isnothing(solvers.h.limit)
#         Tmin = solvers.h.limit[1]; Tmax = solvers.h.limit[2]
#         thermoClamp!(model, h, Tmin, Tmax)
#     end

#     htoT!(model, h, T)
#     interpolate!(hf, h, config)
#     correct_boundaries!(hf, h, boundaries.h, time, config)

#     residuals = (:h, h_res)
#     converged = h_res <= solvers.h.convergence
#     state.residuals = residuals
#     state.converged = converged

#     return nothing
# end


# """
#     thermo_Psi!(model::Physics{T,F,SO,M,Tu,E,D,BI}, Psi::ScalarField) 
#     where {T,F<:AbstractCompressible,M,Tu,E,D,BI}

# Model updates the value of Psi.

# ### Input
# - `model`  -- Physics model defined by user.
# - `Psi`    -- Compressibility factor ScalarField.

# ### Algorithm
# Weakly compressible currently uses the ideal gas equation for establishing the
# compressibility factor where ``\\rho = p * \\Psi``. ``\\Psi`` is calculated from the sensible 
# enthalpy, reference temperature and fluid model specified ``C_p`` and ``R`` value where 
# ``R`` is calculated from ``C_p`` and ``\\gamma`` specified in the fluid model.
# """
# function thermo_Psi!(
#     model::Physics{T,F,SO,M,Tu,E,D,BI}, Psi::ScalarField
#     ) where {T,F<:AbstractCompressible,SO,M,Tu,E,D,BI}
#     (; coeffs, h) = model.energy
#     (; Tref) = coeffs
#     Cp = model.fluid.cp; R = model.fluid.R
#     @. Psi.values = Cp.values/(R.values*(h.values + Cp.values*Tref))
# end

# """
#     thermo_Psi!(model::Physics{T,F,SO,M,Tu,E,D,BI}, Psif::FaceScalarField) 
#     where {T,F<:AbstractCompressible,M,Tu,E,D,BI}

# Function updates the value of Psi.

# ### Input
# - `model`  -- Physics model defined by user.
# - `Psif`    -- Compressibility factor FaceScalarField.

# ### Algorithm
# Weakly compressible currently uses the ideal gas equation for establishing the
# compressibility factor where ``\\rho = p * \\Psi``. ``\\Psi`` is calculated from the sensible 
# enthalpy, reference temperature and fluid model specified ``C_p`` and ``R`` value where 
# ``R`` is calculated from ``C_p`` and ``\\gamma`` specified in the fluid model.
# """
# function thermo_Psi!(
#     model::Physics{T,F,SO,M,Tu,E,D,BI}, Psif::FaceScalarField, config
#     ) where {T,F<:AbstractCompressible,SO,M,Tu,E,D,BI}
#     (; coeffs, hf, h) = model.energy
#     interpolate!(hf, h, config)
#     correct_boundaries!(hf, h, config.boundaries.h, time, config)
#     (; Tref) = coeffs
#     Cp = model.fluid.cp; R = model.fluid.R
#     @. Psif.values = Cp.values/(R.values*(hf.values + Cp.values*Tref))
# end

# """
#     Ttoh!(model::Physics{T1,F,SO,M,Tu,E,D,BI}, T::ScalarField, h::ScalarField
#     ) where {T1,F<:AbstractCompressible,M,Tu,E,D,BI}

# Function coverts temperature ScalarField to sensible enthalpy ScalarField.

# ### Input
# - `model`  -- Physics model defined by user.
# - `T`      -- Temperature ScalarField.
# - `h`      -- Sensible enthalpy ScalarField.
# """
# function Ttoh!(
#     model::Physics{T1,F,SO,M,Tu,E,D,BI}, T::ScalarField, h::ScalarField
#     ) where {T1,F<:AbstractCompressible,SO,M,Tu,E,D,BI}
#     (; coeffs) = model.energy
#     (; Tref) = coeffs
#     Cp = model.fluid.cp
#     @. h.values = Cp.values*(T.values-Tref)
# end

# """
#     htoT!(model::Physics{T1,F,SO,M,Tu,E,D,BI}, h::ScalarField, T::ScalarField
#     ) where {T1,F<:AbstractCompressible,M,Tu,E,D,BI}

# Function coverts sensible enthalpy ScalarField to temperature ScalarField.

# ### Input
# - `model`  -- Physics model defined by user.
# - `h`      -- Sensible enthalpy ScalarField.
# - `T`      -- Temperature ScalarField.
# """
# function htoT!(
#     model::Physics{T1,F,SO,M,Tu,E,D,BI}, h::ScalarField, T::ScalarField
#     ) where {T1,F<:AbstractCompressible,SO,M,Tu,E,D,BI}
#     (; coeffs) = model.energy
#     (; Tref) = coeffs
#     Cp = model.fluid.cp
#     @. T.values = (h.values/Cp.values) + Tref
# end

# function thermoClamp!(
#     model::Physics{T1,F,SO,M,Tu,E,D,BI}, h::ScalarField, Tmin, Tmax
#     ) where {T1,F<:AbstractCompressible,SO,M,Tu,E,D,BI}
#     (; coeffs) = model.energy
#     (; Tref) = coeffs
#     Cp = model.fluid.cp
#     hmin = Cp.values*(Tmin-Tref)
#     hmax = Cp.values*(Tmax-Tref)
#     clamp!(h.values, hmin, hmax)
# end

# function thermoClamp!(
#     model::Physics{T1,F,SO,M,Tu,E,D,BI}, hf::FaceScalarField, Tmin, Tmax
#     ) where {T1,F<:AbstractCompressible,SO,M,Tu,E,D,BI}
#     (; coeffs) = model.energy
#     (; Tref) = coeffs
#     Cp = model.fluid.cp
#     hmin = Cp.values*(Tmin-Tref)
#     hmax = Cp.values*(Tmax-Tref)
#     clamp!(hf.values, hmin, hmax)
# end

# function correct_face_interpolation!(phif::FaceScalarField, phi, Uf::FaceScalarField)
#     mesh = phif.mesh
#     (; faces, cells) = mesh
#     for fID ∈ eachindex(faces)
#         face = faces[fID]
#         (; ownerCells, area, normal) = face
#         cID1 = ownerCells[1]
#         cID2 = ownerCells[2]
#         phi1 = phi[cID1]
#         phi2 = phi[cID2]
#         flux = Uf[fID]
#         if flux >= 0.0
#             phif.values[fID] = phi1
#         else
#             phif.values[fID] = phi2
#         end
#     end
# end
