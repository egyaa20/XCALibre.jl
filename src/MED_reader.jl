module MED_types

using StaticArrays

# Abstract supertype for all cells
abstract type MedCell end

# Concrete cell types
struct TetraCell  <: MedCell
  index::Int
  node_ids::NTuple{4,Int}
end

struct PyramidCell <: MedCell
  index::Int
  node_ids::NTuple{5,Int}
end

struct WedgeCell   <: MedCell
  index::Int
  node_ids::NTuple{6,Int}
end

struct HexaCell   <: MedCell
  index::Int
  node_ids::NTuple{8,Int}
end

# Container for the entire mesh
struct MedMesh{F}
  coords  :: Vector{SVector{3,F}}          # list of node coordinates
  cells   :: Vector{MedCell}               # all volume elements
  groups  :: Dict{String, Vector{Int}}     # named groups (e.g. boundaries)
  # Optional: submeshes, result fields, etc.
end

end # module






function read_MED(path::String)
    # Open HDF5 file
    h5 = h5open(path, "r")
    try
      # Navigate to the mesh group (often /Mesh/0)
      mesh_grp = h5["/Mesh/0"]
  
      # 1) Read node coordinates (Nx3 array)
      raw_coords = read(mesh_grp["Nodes"])            # Array{Float64,2}
      coords = [SVector{3,Float64}(raw_coords[:,i]) for i in 1:size(raw_coords, 2)]
  
      # 2) Read cell connectivity datasets
      types = read(mesh_grp["CellTypes"])             # Vector{Int}
      idxs  = read(mesh_grp["CellIndex"])             # Vector{Int}
      conn  = read(mesh_grp["CellNodeConnectivity"])  # Vector{Int}
  
      # 3) Assemble MedCell instances
      cells = MedCell[]
      for i in 1:length(types)
        start = idxs[i] + 1
        stop  = idxs[i+1]
        node_ids = conn[start:stop]
        c = create_cell(types[i], node_ids, i)
        push!(cells, c)
      end
  
      # 4) (Optional) Read groups under /Mesh/0/Groups
      groups = Dict{String,Vector{Int}}()
      if haskey(mesh_grp, "Groups")
        for name in keys(mesh_grp["Groups"])
          ds = mesh_grp["Groups"][name]
          members = read(ds)
          groups[name] = members
        end
      end
  
      return MedMesh(coords, cells, groups)
    finally
      close(h5)
    end
  end