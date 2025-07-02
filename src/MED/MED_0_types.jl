# MED_0_types.jl - Data Structures for MED Mesh Parser (Intermediate)

# using StaticArrays ## is it needed?

# Point structure (Consistent with UNV)
struct Point{F<:AbstractFloat, SV3<:SVector{3,F}}
    xyz::SV3
end
Point(z::TF) where TF<:AbstractFloat = Point(SVector{3, TF}(zero(TF), zero(TF), zero(TF)))

# Edge (Consistent - Commented out)
# mutable struct Edge...

# Face structure (Consistent with UNV, nnodes used)
mutable struct Face{I<:Integer, VI<:AbstractArray{I}}
    index::I      # Unique identifier FOR THIS PARSING STEP (e.g., sequential reader index)
    nodeCount::I     # Number of nodes defining the face
    nodesID::VI   # Array of node indices (referencing Point structures)
end
Face(z::TI) where TI<:Integer = Face(zero(TI), zero(TI), TI[])

# Cell_MED structure (Consistent with UNV, nnodes used)
mutable struct Cell_MED{I<:Integer,VI<:AbstractArray{I}}
    index::I      # Unique identifier FOR THIS PARSING STEP (e.g., sequential reader index)
    nodeCount::I     # Number of nodes defining the cell
    nodesID::VI   # Array of node indices (referencing Point structures)
end
Cell_MED(z::TI) where TI<:Integer = Cell_MED(zero(TI), zero(TI), TI[])

# BoundaryElement structure
mutable struct BoundaryElement{S<:String,I<:Integer,VI<:AbstractArray{I}}
    name::S       # Name of the boundary group
    index::I    
    facesID::VI   # these are nodes IDs - should probably just call them that
end
BoundaryElement(z::TI) where TI<:Integer = BoundaryElement("default", zero(TI), TI[])

# Element (Consistent - Commented out)
# mutable struct Element...

println("MED Mesh Parser Intermediate Types Loaded.") # Optional print