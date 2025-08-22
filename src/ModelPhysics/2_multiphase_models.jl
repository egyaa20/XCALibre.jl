export Phase, Fluid, Multiphase
export Gravity
export ConstEos, PerfectGas, HelmholtzEnergy, ConstMu, Sutherland, Andrade
export Phase, physicsProperties

export ConstSurfaceTension, SurfaceTensionModel, LeeModel, NucleateBoilingModel, DriftVelocity, Drag_SchillerNaumann


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

Base.@kwdef struct PerfectGas{T<:AbstractFloat} <: AbstractEosModel
    rho::T
    R::T
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

Base.@kwdef struct Sutherland{T<:AbstractFloat} <: AbstractViscosityModel
    mu_ref::T
    S::T
end

Base.@kwdef struct Andrade{T<:AbstractFloat} <: AbstractViscosityModel
    B::T
    C::T
end
Base.@kwdef struct HydrogenViscosity <: AbstractPhysicsProperty end
Base.@kwdef struct NitrogenViscosity <: AbstractPhysicsProperty end

Base.@kwdef struct Gravity{V<:AbstractVector{<:AbstractFloat}} <: AbstractPhysicsProperty
    g::V
end

Base.@kwdef struct ConstSurfaceTension{T<:AbstractFloat} <: AbstractPhysicsProperty
    s::T
end


Base.@kwdef struct LeeModel{T<:AbstractFloat} <: AbstractPhysicsProperty
    evap_coeff::T
    condens_coeff::T
end



Base.@kwdef struct SurfaceTensionModel <: AbstractPhysicsProperty end
Base.@kwdef struct NucleateBoilingModel <: AbstractPhysicsProperty end


Base.@kwdef struct Phase{E<:AbstractEosModel, M<:AbstractViscosityModel} <: AbstractPhase
    eos::E
    mu::M
end

Base.@kwdef struct Drag_SchillerNaumann <: AbstractDrag end

# to be changed later I suppose.... needs to be coupled with subgrid models!
Base.@kwdef struct DriftVelocity{G<:Gravity, T<:AbstractFloat, D<:AbstractDrag} <: AbstractPhysicsProperty
    gravity::G # required for acceleration field computation
    d_p::T #d_p, either computed from subgrid model or defined as const
    drag::D
end





@kwdef struct Multiphase{P1,P2,S1,S2,F1,F2} <: AbstractMultiphase
    phases::P1
    physics_properties::P2
    alpha::S1
    rho::S2
    alphaf::F1
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


function build_multiphase(phases::Tuple{<:AbstractPhase, <:AbstractPhase}, physics_properties::NamedTuple, mesh)
    phase_1, phase_2 = phases
    
    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)
    rho    = ScalarField(mesh)
    rhof   = FaceScalarField(mesh)
    
    Multiphase(phases=phases, physics_properties=physics_properties, alpha=alpha, rho=rho, alphaf=alphaf, rhof=rhof)
end