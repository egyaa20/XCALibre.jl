module MED

using Revise
using StaticArrays
using LinearAlgebra
using Accessors
using Adapt
using Printf
using Statistics
using Setfield
using HDF5

using XCALibre.Mesh

include("MED_0_types.jl")
include("MED_1_reader.jl")
include("MED_2_builder.jl")
include("MED_check_connectivity.jl")

end