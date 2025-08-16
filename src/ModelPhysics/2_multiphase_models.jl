export Phase, Fluid, Multiphase
export Gravity#, gravity
export ConstEos, PerfectGas, HelmholtzEnergy, ConstMu, Sutherland, Andrade
# export constEos, perfectGas, constMu, sutherland, andrade
export Phase, physicsProperties

# export constSurfaceTension, surfaceTensionModel, leeModel, nucleateBoilingModel
export ConstSurfaceTension, SurfaceTensionModel, LeeModel, NucleateBoilingModel


abstract type AbstractModel end
abstract type AbstractEosModel <: AbstractModel end
abstract type AbstractViscosityModel <: AbstractModel end
abstract type AbstractPhysicsProperty end


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
    rho_ref::T
    R::T
end

Base.@kwdef struct HelmholtzEnergy <: AbstractEosModel
    name::Symbol
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


# ConstEos(rho::AbstractFloat) = begin
#     println("Const EoS is loaded")
#     ConstEos(rho)
# end

# PerfectGas(; rho::AbstractFloat, R::AbstractFloat) = begin
#     println("PerfectGas is loaded")
#     PerfectGas(rho, R)
# end

# ConstMu(mu::AbstractFloat) = begin
#     println("Const mu is loaded")
#     ConstMu(mu)
# end


# Gravity(g::AbstractVector{<:AbstractFloat}) = begin
#     println("Gravity is loaded")
#     Gravity(g)
# end
# ConstSurfaceTension(s::AbstractFloat) = begin
#     println("Constant Surface Tension is loaded")
#     ConstSurfaceTension(s)
# end
# LeeModel(; evap_coeff::T, condens_coeff::T) = begin
#     println("Lee Model is loaded")
#     LeeModel(evap_coeff, condens_coeff)
# end


# ConstEos(rho::AbstractFloat) = ConstEos(rho)
# PerfectGas(; rho::AbstractFloat, R::AbstractFloat) = PerfectGas(rho, R)

# HelmholtzEnergy(name::Symbol) = HelmholtzEnergy(name)

# ConstMu(mu::AbstractFloat) = ConstMu(mu)
# Sutherland(; mu_ref::AbstractFloat, S::AbstractFloat) = Sutherland(mu_ref, S)
# Andrade(; B::AbstractFloat, C::AbstractFloat) = Andrade(B, C)

# Gravity(g::AbstractVector{<:AbstractFloat}) = Gravity(g)
# ConstSurfaceTension(s::T) where {T<:AbstractFloat} = ConstSurfaceTension(s)
# SurfaceTensionModel() = SurfaceTensionModel()
# LeeModel(; evap_coeff::T, condens_coeff::T) where {T<:AbstractFloat} = LeeModel(evap_coeff, condens_coeff)
# NucleatBoilingModel() = NucleateBoilingModel()


Base.@kwdef struct Phase{E<:AbstractEosModel, M<:AbstractViscosityModel} <: AbstractPhase
    eos::E
    mu::M
end
# Phase(; eos::AbstractEosModel, mu::AbstractViscosityModel) = Phase(eos, mu)


# phase_1(; eos::AbstractEosModel, mu::AbstractViscosityModel) = Phase(eos, mu)
# phase_2(; eos::AbstractEosModel, mu::AbstractViscosityModel) = Phase(eos, mu)







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
    println("Phases Dispatch")
    
    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)
    rho    = ScalarField(mesh)
    rhof   = FaceScalarField(mesh)
    
    Multiphase(phases=phases, physics_properties=physics_properties, alpha=alpha, rho=rho, alphaf=alphaf, rhof=rhof)
end