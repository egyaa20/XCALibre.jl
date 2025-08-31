export Phase, Fluid, Multiphase
export Gravity
export ConstEos, PerfectGas, HelmholtzEnergy, ConstMu, Sutherland, Andrade
export Phase, physicsProperties

export ConstSurfaceTension, SurfaceTensionModel, LeeModel, NucleateBoilingModel, DriftVelocity, Drag_SchillerNaumann
export AbstractModel, AbstractEosModel, AbstractViscosityModel, AbstractPhysicsProperty

export LeeModelState, DriftVelocityState, GravityState



abstract type AbstractModel end
abstract type AbstractEosModel <: AbstractModel end
abstract type AbstractViscosityModel <: AbstractModel end
abstract type AbstractPhysicsProperty end


abstract type AbstractDrag <: AbstractPhysicsProperty end


abstract type AbstractPhase <: AbstractMultiphase end

# abstract type AbstractMultiphaseModel end


Base.@kwdef struct ConstEos{T<:AbstractFloat} <: AbstractEosModel
    rho::T
end
(eos::ConstEos)(phase, model) = begin
    rho_field = phase.rho
    initialise!(rho_field, eos.rho)
end


Base.@kwdef struct PerfectGas{T<:AbstractFloat} <: AbstractEosModel
    rho::T
    R::T
end
(eos::PerfectGas)(phase, model) = begin
    (; p) = model.momentum

    T_ref = 273.0
    R = phase.eosModel.R
    rho_field = phase.rho

    @. rho_field.values = (p.values) / (R * T_ref) # CAREFUL WITH p=0 initialisation
end


Base.@kwdef struct PengRobinson{T<:AbstractFloat} <: AbstractEosModel # edit
    T_crit::T
    p_crit::T
    omega::T
    M::T # Molar mass in g/mol
end


abstract type HelmholtzEnergyFluid end
struct N2 <: HelmholtzEnergyFluid end
struct H2 <: HelmholtzEnergyFluid end

Base.@kwdef struct HelmholtzEnergy{F<:HelmholtzEnergyFluid, I<:Bool} <: AbstractEosModel
    name::F
    interpolationMode::I
end



Base.@kwdef struct ConstMu{T<:AbstractFloat} <: AbstractViscosityModel
    mu::T
end
(viscosityModel::ConstMu)(phase, T) = begin
    mu_field = phase.mu
    initialise!(mu_field, viscosityModel.mu)
end


Base.@kwdef struct Sutherland{T<:AbstractFloat} <: AbstractViscosityModel
    mu_ref::T
    S::T
end
(viscosityModel::Sutherland)(phase, T) = begin
    mu_ref = viscosityModel.mu_ref
    S = viscosityModel.S

    T_ref = 273.0
    mu_field = phase.mu

    @. mu_field.values = (mu_ref * (T.values/T_ref)^1.5) * ((T_ref + S)/(T.values + S))
end


Base.@kwdef struct Andrade{T<:AbstractFloat} <: AbstractViscosityModel
    B::T
    C::T
end
(viscosityModel::Andrade)(phase, T) = begin
    B = viscosityModel.B
    C = viscosityModel.C
    
    mu_field = phase.mu

    @. mu_field.values = B * exp(C / T.values)
end


Base.@kwdef struct HydrogenViscosity <: AbstractPhysicsProperty end
Base.@kwdef struct NitrogenViscosity <: AbstractPhysicsProperty end


Base.@kwdef struct Gravity{V<:AbstractVector{<:AbstractFloat}} <: AbstractPhysicsProperty
    g::V
end
@kwdef struct GravityState{V<:AbstractVector{<:AbstractFloat}} <: AbstractPhysicsProperty
    g::V
    momentum_source_pointer::Ref{Int}
end
Adapt.@adapt_structure GravityState

function build_gravityModel(setup::Gravity, mesh)
    return GravityState(
        g=setup.g,
        momentum_source_pointer=Ref(0)
    )
end


Base.@kwdef struct ConstSurfaceTension{T<:AbstractFloat} <: AbstractPhysicsProperty
    s::T
end


Base.@kwdef struct LeeModel{T<:AbstractFloat} <: AbstractPhysicsProperty
    evap_coeff::T
    condens_coeff::T
end
@kwdef struct LeeModelState{T<:AbstractFloat, M1,M2,LH} <: AbstractPhysicsProperty
    evap_coeff::T
    condens_coeff::T
    m_qp::M1
    m_pq::M2
    latentHeat::LH
    alpha_source_pointer::Ref{Int}
    energy_source_pointer::Ref{Int}
end
Adapt.@adapt_structure LeeModelState

function build_leeModel(setup::LeeModel, mesh)
    m_qp        = ScalarField(mesh)
    m_pq        = ScalarField(mesh)
    latentHeat  = ScalarField(mesh)

    return LeeModelState(
        evap_coeff=setup.evap_coeff,
        condens_coeff=setup.condens_coeff,
        m_qp=m_qp,
        m_pq=m_pq,
        latentHeat=latentHeat,
        alpha_source_pointer=Ref(0),
        energy_source_pointer=Ref(0)
    )
end



Base.@kwdef struct SurfaceTensionModel <: AbstractPhysicsProperty end
Base.@kwdef struct NucleateBoilingModel <: AbstractPhysicsProperty end






@kwdef struct Phase{E<:AbstractEosModel, V<:AbstractViscosityModel} <: AbstractPhase
    eosModel::E
    viscosityModel::V
end
@kwdef struct PhaseState{E<:AbstractEosModel, V<:AbstractViscosityModel, S1,S2,S3,S4,S5} <: AbstractPhase
    eosModel::E
    viscosityModel::V

    rho::S1
    mu::S2
    k::S3
    cp::S4
    beta::S5
end
Adapt.@adapt_structure PhaseState

function build_phase(phase_setup::Phase, mesh)
    rho   = ScalarField(mesh)
    mu    = ScalarField(mesh)
    k     = ScalarField(mesh)
    cp    = ScalarField(mesh)
    beta  = ScalarField(mesh)

    return PhaseState(
        eosModel=phase_setup.eosModel,
        viscosityModel=phase_setup.viscosityModel,
        rho=rho,
        mu=mu,
        k=k,
        cp=cp,
        beta=beta
    )
end




Base.@kwdef struct Drag_SchillerNaumann <: AbstractDrag end


Base.@kwdef struct DriftVelocity{G<:Gravity, T<:AbstractFloat, D<:AbstractDrag} <: AbstractPhysicsProperty
    gravity::G      # required for acceleration field computation
    d_p::T          # d_p, either computed from subgrid model or defined as const
    drag::D         # Drag_SchillerNaumann for now
end
@kwdef struct DriftVelocityState{G<:Gravity, T<:AbstractFloat, D<:AbstractDrag, V1,V2,V3,V4,V5,V6,V7} <: AbstractPhysicsProperty
    gravity::G 
    d_p::T
    drag::D

    v_dr_p::V1
    v_dr_q::V2
    v_p::V3
    v_q::V4

    v_pq::V5
    v_pq_prev::V6
    U_prev::V7

    momentum_source_pointer::Ref{Int}
    alpha_source_pointer::Ref{Int}
    energy_source_pointer::Ref{Int}
end
Adapt.@adapt_structure DriftVelocityState

function build_driftVelocity(setup::DriftVelocity, mesh)
    v_dr_p  = VectorField(mesh)
    v_dr_q  = VectorField(mesh)
    v_p     = VectorField(mesh)
    v_q     = VectorField(mesh)

    v_pq      = VectorField(mesh)
    v_pq_prev = VectorField(mesh)
    U_prev    = VectorField(mesh)

    return DriftVelocityState(
        gravity=setup.gravity,
        d_p=setup.d_p,
        drag=setup.drag,
        v_dr_p=v_dr_p,
        v_dr_q=v_dr_q,
        v_pq=v_pq,
        v_pq_prev=v_pq_prev,
        U_prev=U_prev,
        v_p=v_p,
        v_q=v_q,
        momentum_source_pointer=Ref(0),
        alpha_source_pointer=Ref(0),
        energy_source_pointer=Ref(0)
    )
end





@kwdef struct Multiphase{P1,P2,S1,F1,S2,F2} <: AbstractMultiphase
    phases::P1
    physics_properties::P2
    alpha::S1
    alphaf::F1
    rho::S2
    rhof::F2
end
Adapt.@adapt_structure Multiphase

Fluid{Multiphase}(; phases::Tuple, kwargs...) = begin
    coeffs = (; phases, kwargs...)
    ARG = typeof(coeffs)
    Fluid{Multiphase, ARG}(coeffs)
end


(fluid::Fluid{Multiphase, ARG})(mesh) where {ARG} = begin
    coeffs = fluid.args

    physics_properties = Base.structdiff(coeffs, (phases = nothing,))

    build_multiphase(coeffs.phases, physics_properties, mesh)
end


build_property(property, mesh) = property
build_property(setup::Gravity, mesh) = build_gravityModel(setup, mesh)
build_property(setup::LeeModel, mesh) = build_leeModel(setup, mesh)
build_property(setup::DriftVelocity, mesh) = build_driftVelocity(setup, mesh)

function build_multiphase(phase_setups::Tuple{<:AbstractPhase, <:AbstractPhase}, physics_properties_setup::NamedTuple, mesh)
    phases = map(setup -> build_phase(setup, mesh), phase_setups)

    built_properties = map(prop_setup -> build_property(prop_setup, mesh), physics_properties_setup)

    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)

    rho  = ScalarField(mesh)
    rhof = FaceScalarField(mesh)
    
    Multiphase(phases=phases, physics_properties=built_properties, alpha=alpha, alphaf=alphaf, rho=rho, rhof=rhof)
end