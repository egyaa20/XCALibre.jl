export liquid, gas, mixture
export Water, Air, H2, N2
export gravity
export constEos, perfectGas, constMu, sutherland, andrade


abstract type AbstractPhase <: AbstractMultiphase end

abstract type AbstractLiquid <: AbstractPhase end
abstract type AbstractGas <: AbstractPhase end
abstract type AbstractMixture <: AbstractPhase end

struct Water <: AbstractLiquid end
struct Air <: AbstractGas end
struct H2 <: AbstractMixture end
struct N2 <: AbstractMixture end


abstract type AbstractEosModel end
abstract type AbstractViscosityModel end

abstract type AbstractPhysicsProperty end

@kwdef struct Gravity{V<:Vector{<:Float64}} <: AbstractPhysicsProperty
    g::V = [0.0, -9.81, 0.0]
end

# EOS models
struct ConstEos{T} <: AbstractEosModel; rho::T; end
struct PerfectGas <: AbstractEosModel end

# VISCOSITY models
struct ConstMu{T} <: AbstractViscosityModel; mu::T; end
struct Sutherland <: AbstractViscosityModel end
struct Andrade <: AbstractViscosityModel end



struct LiquidPhase{N,E,M} <: AbstractLiquid
    name::N
    eos::E
    mu::M
end

struct GasPhase{N,E,M,R} <: AbstractGas
    name::N
    eos::E
    mu::M
    R::R
end

struct MixturePhase{T<:AbstractMixture} <: AbstractMixture
    name::T
end


liquid(; name::N, eos::E, mu::M) where {
    N<:AbstractLiquid,
    E<:ConstEos,
    M<:Union{ConstMu, Andrade}
} = LiquidPhase(name, eos, mu)

gas(; name::N, eos::E, mu::M, R::Float64) where {
    N<:AbstractGas,
    E<:Union{ConstEos, PerfectGas},
    M<:Union{ConstMu, Sutherland}
} = GasPhase(name, eos, mu,R)

mixture(; name::T) where {T<:AbstractMixture} = MixturePhase(name)

gravity(g::Vector{<:Float64} = [0.0, -9.81, 0.0]) = Gravity(g)

# EOS models
constEos(; rho) = ConstEos(rho)
perfectGas() = PerfectGas()

# VISCOSITY models
constMu(; mu) = ConstMu(mu)
sutherland() = Sutherland()
andrade() = Andrade()


@kwdef struct Multiphase{P1,P2,S1,S2,F1,F2} <: AbstractMultiphase
    phases::P1
    physicsProperties::P2
    alpha::S1
    rho::S2
    alphaf::F1
    rhof::F2
end
Adapt.@adapt_structure Multiphase

Fluid{Multiphase}(; phases::Tuple, physicsProperties::Tuple = ()) = begin
    coeffs = (; phases = phases, physicsProperties = physicsProperties)
    ARG = typeof(coeffs)
    Fluid{Multiphase, ARG}(coeffs)
end


(fluid::Fluid{Multiphase, ARG})(mesh) where {ARG} = begin
    coeffs = fluid.args
    build_multiphase(coeffs.phases, coeffs.physicsProperties, mesh)
end


function build_multiphase(phases::Tuple{<:LiquidPhase, <:GasPhase}, physicsProperties::Tuple, mesh)
    println("Dispatching to: L and G")
    liquid_phase, gas_phase = phases
    println(" -> Liquid: $(typeof(liquid_phase.name))")
    println(" -> Gas: $(typeof(gas_phase.name))")
    # L and G
    
    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)
    rho    = ScalarField(mesh)
    rhof   = FaceScalarField(mesh)
    
    Multiphase(phases=phases, physicsProperties=physicsProperties, alpha=alpha, rho=rho, alphaf=alphaf, rhof=rhof)
end

function build_multiphase(phases::Tuple{<:MixturePhase}, physicsProperties::Tuple, mesh)
    println("Dispatching to: Mixture")
    mixture_phase = phases[1]
    println(" -> Mixture: $(typeof(mixture_phase.name))")
    #mixture....

    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)
    rho    = ScalarField(mesh)
    rhof   = FaceScalarField(mesh)
    
    Multiphase(phases=mixture_phase, physicsProperties=physicsProperties, alpha=alpha, rho=rho, alphaf=alphaf, rhof=rhof)
end