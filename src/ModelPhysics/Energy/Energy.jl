# module Energy

using Atomix
using KernelAbstractions
using Accessors
using StaticArrays
using LinearAlgebra
using Adapt
using XCALibre.Mesh
using XCALibre.Fields
using XCALibre.ModelFramework
using XCALibre.Discretise
using XCALibre.ModelPhysics
using XCALibre.Solve
using XCALibre.Calculate
using XCALibre.IOFormats

include("energy_types.jl")

# Energy models
include("Sensible_Enthalpy.jl")
include("Multiphase_Energy.jl")
include("Conduction.jl")



# Property Models
include("PropertyModels/Cryogenic_metal_properties.jl")

# Equations Of State
include("PropertyModels/EquationsOfState/multiparameter_H2.jl")
include("PropertyModels/EquationsOfState/multiparameter_N2.jl")

# Viscosity Models
include("PropertyModels/Viscosity/high_fidelity_mu_H2.jl")
include("PropertyModels/Viscosity/high_fidelity_mu_N2.jl")

# Thermal Conductivity
include("PropertyModels/ThermalConductivity/thermal_conductivity_H2.jl")
include("PropertyModels/ThermalConductivity/thermal_conductivity_N2.jl")
# Surface Tension
include("PropertyModels/surface_tension.jl")

# Top level functor
include("PropertyModels/EOS_closures.jl")


export initialise, energy!

# end # end module