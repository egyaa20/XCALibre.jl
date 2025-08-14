# export liquid, gas, mixture
# export Water, Air, H2, N2
# export gravity
# export constEos, perfectGas, constMu, sutherland, andrade
# export Andrade, Sutherland, PerfectGas, ConstEos, ConstMu
# export AbstractPhase_1, AbstractPhase_2



# abstract type AbstractPhase <: AbstractMultiphase end

# abstract type AbstractPhase_1 <: AbstractPhase end
# abstract type AbstractPhase_2 <: AbstractPhase end
# # abstract type AbstractMixture <: AbstractPhase end

# # struct CustomLiquid <: AbstractLiquid end
# # struct Water <: AbstractLiquid end

# # struct CustomGas <: AbstractGas end
# # struct Air <: AbstractGas end

# # struct H2 <: AbstractMixture end
# # struct N2 <: AbstractMixture end


# abstract type AbstractEosModel end
# abstract type AbstractViscosityModel end

# abstract type AbstractPhysicsProperty end

# @kwdef struct Gravity{V<:Vector{<:Float64}} <: AbstractPhysicsProperty
#     g::V
# end


# struct ConstEos <: AbstractEosModel; rho::Float64; end

# struct ConstMu <: AbstractViscosityModel; mu::Float64; end


# struct PerfectGas <: AbstractEosModel
#     R::Union{Float64,Nothing}
# end
# perfectGas() = PerfectGas(nothing)
# perfectGas(; R::Float64) = PerfectGas(R)

# struct Sutherland <: AbstractViscosityModel
#     mu_ref::Union{Float64,Nothing}  # Pa*s
#     S::Union{Float64,Nothing}       # K
# end
# # Sutherland() = Sutherland(nothing, nothing)
# Sutherland(; mu_ref::Float64, S::Float64) = Sutherland(mu_ref, S)
# Sutherland() = Sutherland(nothing, nothing)


# struct Andrade <: AbstractViscosityModel
#     B::Union{Float64,Nothing}
#     C::Union{Float64,Nothing}
# end
# andrade() = Andrade(nothing, nothing)
# andrade(; B::Float64, C::Float64) = Andrade(B, C)



# struct LiquidPhase{N,E,T} <: AbstractLiquid
#     name::N
#     eos::E
#     transport::T
# end

# struct GasPhase{N,E,T,R} <: AbstractGas
#     name::N
#     eos::E
#     transport::T
# end

# # struct MixturePhase{T<:AbstractMixture} <: AbstractMixture
# #     name::T
# # end


# # default_eos(::Water) = ConstEos(1000.0)
# # default_eos(::Air) = ConstEos(1.225)

# # default_transport(::Water) = ConstMu(1.0e-3)
# # default_transport(::Air) = ConstMu(1.8e-5)

# # default_R(::Air) = 287.0
# # default_perfectGas(::Air) = PerfectGas(default_R(Air()))

# # default_sutherland(::Air) = Sutherland(1.8e-5, 110.4)
# # default_andrade(::Water) = Andrade(1.732e-6, 1863.0)


# liquid(; name::Union{AbstractLiquid,Nothing}=nothing,
#          eos::Union{ConstEos,Nothing}=nothing,
#          transport::Union{ConstMu,Andrade,Nothing}=nothing) = begin
#     if name !== nothing
#         if eos === nothing
#             eos_val = default_eos(name)
#         else
#             eos_val = eos
#         end

#         if transport === nothing
#             transport_val = default_transport(name)
#         else
#             if transport isa Andrade && (isnothing(transport.B) || isnothing(transport.C))
#                 transport_val = default_andrade(name)
#             else
#                 transport_val = transport
#             end
#         end

#         return LiquidPhase(name, eos_val, transport_val)
#     else
#         if eos === nothing || transport === nothing
#             println("NOT ENOUGH ARGS PROVIDED")
#         end

#         return LiquidPhase(CustomLiquid(), eos, transport)
#     end
# end


# # liquid(; name::N, ## THIS DOES NOT WORK
# #          eos::ConstEos = default_eos(name),
# #          transport::Union{ConstMu,Andrade} = default_transport(name)
# # ) where {N<:AbstractLiquid} = LiquidPhase(name, eos, transport)

# # liquid(; eos::ConstEos,
# #          transport::Union{ConstMu,Andrade}
# # ) = LiquidPhase(CustomLiquid(), eos, transport)

# gas(; name::Union{AbstractGas,Nothing}=nothing,
#       eos::Union{ConstEos,PerfectGas,Nothing}=nothing,
#       transport::Union{ConstMu,Sutherland,Nothing}=nothing,
#       R::Union{Float64,Nothing}=nothing) = begin
#     if name !== nothing
#         if eos === nothing
#             eos_val = default_eos(name)
#         else
#             if transport isa PerfectGas && (isnothing(eos.R))
#                 eos_val = default_perfectGas(name)
#             else
#                 eos_val = eos
#             end
#         end

#         if transport === nothing
#             transport_val = default_transport(name)
#             # if typeof(transport) <: Sutherland ### THIS ONE DOES NOT REALLY WORK!
#             #     transport_val = default_sutherland(name)
#             # else
#             #     transport_val = default_transport(name)
#             # end
#         else
#             println("else section")
#             if transport isa Sutherland && (isnothing(transport.mu_ref) || isnothing(transport.S))
#                 println("right")
#                 transport_val = default_sutherland(name)
#             else
#                 println("wrong")
#                 transport_val = transport
#             end
#         end

#         if R === nothing
#             R_val = default_R(name)
#         else
#             R_val = R
#         end

#         return GasPhase(name, eos_val, transport_val, R_val)
#     else
#         if eos === nothing || transport === nothing || R === nothing
#             println("NOT ENOUGH ARGS PROVIDED")
#         end

#         return GasPhase(CustomGas(), eos, transport, R)
#     end
# end


# # gas(; name::N,
# #       eos::Union{ConstEos,PerfectGas}=default_eos(name),
# #       transport::Union{ConstMu,Sutherland}=default_transport(name),
# #       R::Float64 = default_R(name)
# # ) where {N<:AbstractGas} = GasPhase(name, eos, transport, R)

# # gas(; eos::Union{ConstEos,PerfectGas},
# #       transport::Union{ConstMu,Sutherland},
# #       R::Float64
# # ) = GasPhase(CustomGas(), eos, transport, R)












# # mixture(; name::T) where {T<:AbstractMixture} = MixturePhase(name)

# gravity(g::Vector{<:Float64}) = Gravity(g)

# # EOS models
# constEos(rho::Float64) = ConstEos(rho)
# # perfectGas() = PerfectGas()

# # VISCOSITY models
# constMu(mu::Float64) = ConstMu(mu)

# # sutherland(mu::Float64, S::Float64) = Sutherland(mu, S)
# # andrade(B::Float64, C::Float64) = Andrade(B, C)




# gas(name=Air(), eos=PerfectGas(1.0), transport=Sutherland()),
















export Phase, Fluid, Multiphase
export Gravity, gravity
export ConstEos, PerfectGas, HelmholtzEnergy, ConstMu, Sutherland, Andrade
export constEos, perfectGas, constMu, sutherland, andrade
export phase_1, phase_2, physicsProperties

export constSurfaceTension, surfaceTensionModel, leeModel, nucleateBoilingModel
export ConstSurfaceTension, SurfaceTensionModel, LeeModel, NucleateBoilingModel


abstract type AbstractModel end
abstract type AbstractEosModel <: AbstractModel end
abstract type AbstractViscosityModel <: AbstractModel end
abstract type AbstractPhysicsProperty end


abstract type AbstractPhase <: AbstractMultiphase end

# abstract type AbstractMultiphaseModel end


struct ConstEos{T<:AbstractFloat} <: AbstractEosModel
    rho::T
end

struct PerfectGas{T<:AbstractFloat} <: AbstractEosModel
    rho_ref::T
    R::T
end

struct HelmholtzEnergy <: AbstractEosModel
    name::Symbol
end

struct ConstMu{T<:AbstractFloat} <: AbstractViscosityModel
    mu::T
end

struct Sutherland{T<:AbstractFloat} <: AbstractViscosityModel
    mu_ref::T
    S::T
end

struct Andrade{T<:AbstractFloat} <: AbstractViscosityModel
    B::T
    C::T
end

struct Gravity{V<:AbstractVector{<:AbstractFloat}} <: AbstractPhysicsProperty
    g::V
end

struct ConstSurfaceTension{T<:AbstractFloat} <: AbstractPhysicsProperty
    s::T
end

struct LeeModel{T<:AbstractFloat} <: AbstractPhysicsProperty
    evap_coeff::T
    condens_coeff::T
end

struct SurfaceTensionModel <: AbstractPhysicsProperty end
struct NucleateBoilingModel <: AbstractPhysicsProperty end



constEos(rho::AbstractFloat) = ConstEos(rho)
perfectGas(; rho::AbstractFloat, R::AbstractFloat) = PerfectGas(rho, R)
helmholtzEnergy(name::Symbol) = HelmholtzEnergy(name)

constMu(mu::AbstractFloat) = ConstMu(mu)
sutherland(; mu_ref::AbstractFloat, S::AbstractFloat) = Sutherland(mu_ref, S)
andrade(; B::AbstractFloat, C::AbstractFloat) = Andrade(B, C)

gravity(g::AbstractVector{<:AbstractFloat}) = Gravity(g)
constSurfaceTension(s::T) where {T<:AbstractFloat} = ConstSurfaceTension(s)
surfaceTensionModel() = SurfaceTensionModel()
leeModel(; evap_coeff::T, condens_coeff::T) where {T<:AbstractFloat} = LeeModel(evap_coeff, condens_coeff)
nucleatBoilingModel() = NucleateBoilingModel()


struct Phase1{E<:AbstractEosModel, M<:AbstractViscosityModel} <: AbstractPhase
    eos::E
    mu::M
end


struct Phase2{E<:AbstractEosModel, M<:AbstractViscosityModel} <: AbstractPhase
    eos::E
    mu::M
end

phase_1(; eos::AbstractEosModel, mu::AbstractViscosityModel) = Phase1(eos, mu)
phase_2(; eos::AbstractEosModel, mu::AbstractViscosityModel) = Phase2(eos, mu)







@kwdef struct Multiphase{P1,P2,S1,S2,F1,F2} <: AbstractMultiphase
    phases::P1
    physics_properties::P2 # QUESTION!!!!!!!!
    alpha::S1
    rho::S2
    alphaf::F1
    rhof::F2
end
Adapt.@adapt_structure Multiphase

# Fluid{Multiphase}(; phases::Tuple, physicsProperties::Tuple = ()) = begin
#     coeffs = (; phases = phases, physicsProperties = physicsProperties)
#     ARG = typeof(coeffs)
#     Fluid{Multiphase, ARG}(coeffs)
# end
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
    # println(" -> Liquid: $(typeof(liquid_phase.name))")
    # println(" -> Gas: $(typeof(gas_phase.name))")
    # L and G
    
    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)
    rho    = ScalarField(mesh)
    rhof   = FaceScalarField(mesh)
    
    Multiphase(phases=phases, physics_properties=physics_properties, alpha=alpha, rho=rho, alphaf=alphaf, rhof=rhof)
end

# function build_multiphase(phases::Tuple{<:MixturePhase}, physicsProperties::Tuple, mesh)
#     println("Dispatching to: Mixture")
#     mixture_phase = phases[1]
#     println(" -> Mixture: $(typeof(mixture_phase.name))")
#     #mixture....

#     alpha  = ScalarField(mesh)
#     alphaf = FaceScalarField(mesh)
#     rho    = ScalarField(mesh)
#     rhof   = FaceScalarField(mesh)
    
#     Multiphase(phases=mixture_phase, physicsProperties=physicsProperties, alpha=alpha, rho=rho, alphaf=alphaf, rhof=rhof)
# end