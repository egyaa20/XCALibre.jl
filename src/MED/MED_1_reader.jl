using HDF5
using StaticArrays

# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

"""
    read_MED_HDF5(path; scale = 1.0, integer = Int64, float = Float64)

Parse a **MED‑HDF5** file (single mesh, single step) and return four vectors
that can be fed straight into *XCALibre.jl*’s builder:

```julia
points      :: Vector{Point}
faces       :: Vector{Face}
cells       :: Vector{Cell_MED}
boundaries  :: Vector{BoundaryElement}
```

The routine:
1. Loads coordinates *as written*, with optional uniform `scale`.
2. Maps solver‑original node IDs (`NOE/NUM`) to coordinate row indices.
3. Walks every `/MAI/<elemType>` group, remapping & (optionally) re‑ordering
   connectivity so that bricks/tets/quads have positive orientation.
4. Collects faces vs. cells according to `_ELEM_KIND`. 2‑node segments (`SE2`,
   `SEG2`) are treated as faces so that 2‑D meshes work out of the box.
5. Groups faces by *family* using `/FAS/<mesh>/ELEME`, yielding one
   `BoundaryElement` per family ID (or `"Unnamed_Family_$id"`).

Throws a descriptive `error()` on any size mismatch, unknown element type, or
missing dataset.
"""
function read_MED_HDF5(path; scale = 1.0, integer = Int64, float = Float64)
    fid = h5open(path, "r")
    try
        mesh_name  = first(keys(fid["/ENS_MAA"]))
        n_meshes    = length(keys(fid["/ENS_MAA"]))
        n_meshes > 1 && @warn "MEDReader: only first mesh is read (found $n_meshes)."

        step_name  = first(keys(fid["/ENS_MAA/$mesh_name"]))

        # ──────────── 1. Coordinates ────────────
        coords = read(fid["/ENS_MAA/$mesh_name/$step_name/NOE/COO"])
        coords = ndims(coords) == 1 ? reshape(coords, 3, :)' :
                 (size(coords, 2) == 3 ? coords : coords')
        xyz = scale .* coords
        points = [Point(SVector{3, float}(xyz[i,1], xyz[i,2], xyz[i,3])) for i in 1:size(xyz,1)]

        # ──────────── 2. Node‑ID ↔ row index map ────────────
        has_num = haskey(fid, "/ENS_MAA/$mesh_name/$step_name/NOE/NUM")
        orig_ids = has_num ? vec(read(fid["/ENS_MAA/$mesh_name/$step_name/NOE/NUM"])) : collect(1:length(points))
        node_id_to_idx = Dict{integer, Int}(orig_ids[i] => i for i in eachindex(orig_ids))

        # ──────────── 3. Connectivity pass ────────────
        faces  = Vector{Face}(undef,0)
        cells  = Vector{Cell_MED}(undef,0)
        fam_map_faces = Dict{integer, Vector{Int}}()

        mai_root = "/ENS_MAA/$mesh_name/$step_name/MAI"
        for elem_type in keys(fid[mai_root])
            g = fid["$mai_root/$elem_type"]
            F  = _ELEM_NODES[elem_type]  # may throw if unsupported
            raw = read(g["NOD"])
            # read attribute "NBR" via HDF5.jl API
            nbr = read_attribute(g["NOD"], "NBR")
            length(raw) == F * nbr || error("Corrupt NOD length for $elem_type: expected $(F*nbr), got $(length(raw))")
            conn = reshape(raw, F, nbr)'
            # Remap node IDs
            conn = map(x -> node_id_to_idx[x], conn)
            # Optional MED → UNV permutation
            if haskey(_PERM_MED2UNV, elem_type)
                perm = _PERM_MED2UNV[elem_type]
                conn = conn[:, perm]
            end
            fam = haskey(g, "FAM") ? vec(read(g["FAM"])) : fill(0, nbr)

            kind = _ELEM_KIND[elem_type]  # :cell or :face
            if kind === :cell
                start_idx = length(cells) + 1
                for i in 1:nbr
                    push!(cells, Cell_MED(start_idx + i - 1, F, conn[i,:]))
                end
            else # :face
                start_idx = length(faces) + 1
                for i in 1:nbr
                    idx = start_idx + i - 1
                    push!(faces, Face(idx, F, conn[i,:]))
                    push!(get!(fam_map_faces, fam[i], Int[]), idx)
                end
            end
        end

        # ──────────── 4. Boundary families ────────────
        fas_root = "/FAS/$mesh_name/ELEME"
        fam_to_name = Dict{integer,String}()
        if haskey(fid, fas_root)
            for grp in keys(fid[fas_root])
                if occursin(r"^FAM_", grp)
                    m = match(r"^FAM_(-?\d+)_([^/]+)$", grp)
                    if m !== nothing
                        fam_to_name[parse(Int, m.captures[1])] = m.captures[2]
                    end
                end
            end
        end

        boundaries = Vector{BoundaryElement}(undef,0)
        for (fam_id, face_ids) in fam_map_faces
            name = get(fam_to_name, fam_id, "Unnamed_Family_$(fam_id)")
            push!(boundaries, BoundaryElement(name, length(boundaries) + 1, face_ids))
        end
        sort!(boundaries; by = x -> x.name)

        return points, faces, cells, boundaries
    finally
        close(fid)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Internals & constants
# ──────────────────────────────────────────────────────────────────────────────

const _ELEM_NODES = Dict(
    "HE8"  => 8,
    "TE4"  => 4,
    "QU4"  => 4,
    "TR3"  => 3,
    "SE2"  => 2,
    "SEG2" => 2
)

const _PERM_MED2UNV = Dict(
    "HE8"  => [1,5,7,3,2,6,8,4],
    "QU4"  => [1,3,4,2],
    "TE4"  => [1,3,2,4],
    "TR3"  => [1,3,2],
    "SE2"  => [1,2],
    "SEG2" => [1,2]
)

const _ELEM_KIND = Dict(
    "HE8"  => :cell,
    "TE4"  => :cell,
    "QU4"  => :face,
    "TR3"  => :face,
    "SE2"  => :face,
    "SEG2" => :face
)



















# using HDF5
# using StaticArrays

# # Helper function to determine expected nodes per element based on MED type string
# function expectedNodesPerElement(elemType::AbstractString)::Union{Int, Nothing}
#     # Based on common MED element types and Appendix B
#     # Extend this list as needed for other element types.
#     if elemType == "TR3" return 3
#     elseif elemType == "QU4" return 4
#     elseif elemType == "TE4" return 4
#     elseif elemType == "HE8" || elemType == "HEXA8" return 8
#     # Add other types like segments (SE2), pyramids (PY5), prisms (PE6) if required
#     elseif elemType == "SE2" return 2 # Segment
#     elseif elemType == "PY5" return 5 # Pyramid
#     elseif elemType == "PE6" return 6 # Prism (Wedge)
#     # ... add more mappings ...
#     else
#         @warn "Unknown element type encountered: $elemType. Cannot determine nodes per element."
#         return nothing
#     end
# end

# # Helper function to check if an element type typically represents a surface element (Face)
# function isFaceElementType(elemType::AbstractString)::Bool
#     # Add more face types if necessary (e.g., TR6, QU8)
#     return elemType in ("TR3", "QU4", "SE2") # Assuming SE2 is boundary line in 2D/3D
# end

# # Helper function to get node permutation from MED to UNV convention (Appendix B)
# function getNodePermutation(elemType::AbstractString)::Union{Vector{Int}, Nothing}
#     if elemType == "HE8" || elemType == "HEXA8" return [1, 5, 7, 3, 2, 6, 8, 4]
#     elseif elemType == "QU4" return [1, 3, 4, 2]
#     elseif elemType == "TE4" return [1, 3, 2, 4]
#     elseif elemType == "TR3" return [1, 3, 2]
#     # Add other permutations if known and needed
#     else
#         return nothing # No permutation needed or known
#     end
# end

# # Helper function to read connectivity and family data, handling Fortran order
# function read_connectivity_and_families(group::HDF5.Group, elemType::String, F::Int, I::Type{<:Integer})
#     local raw_conn_vec, fam_vec, num_elements

#     # --- Read Connectivity (NOD) ---
#     if !haskey(group, "NOD")
#         error("Element block '$elemType' is missing the mandatory 'NOD' dataset.")
#     end
#     nod_dset = group["NOD"]
#     raw_conn_vec = read(nod_dset) # Reads as a flat vector (Fortran order)

#     # Check element count consistency if NBR attribute exists
#     expected_elements = -1
#     if haskey(HDF5.attributes(nod_dset), "NBR")
#         expected_elements = read(HDF5.attributes(nod_dset)["NBR"])
#         if length(raw_conn_vec) != F * expected_elements
#              @warn "Mismatch between NOD vector length ($(length(raw_conn_vec))) and expected size ($F * $expected_elements) based on NBR attribute for '$elemType'. Using vector length."
#              # Decide recovery strategy: error or proceed based on vector length?
#              # Proceeding based on actual length, but requires integer number of elements.
#              if length(raw_conn_vec) % F != 0
#                  error("NOD vector length ($(length(raw_conn_vec))) for '$elemType' is not a multiple of nodes per element ($F). File is likely corrupt.")
#              end
#              num_elements = div(length(raw_conn_vec), F)
#         else
#              num_elements = expected_elements
#         end
#     else
#         # Infer element count from vector length
#         if length(raw_conn_vec) % F != 0
#             error("NOD vector length ($(length(raw_conn_vec))) for '$elemType' is not a multiple of nodes per element ($F). File is likely corrupt.")
#         end
#         num_elements = div(length(raw_conn_vec), F)
#     end

#     if num_elements == 0
#         # Return empty arrays if the block has no elements
#         return Matrix{I}(undef, 0, F), Vector{I}(undef, 0), 0
#     end

#     # Reshape and transpose to Julia's row-major format (E x F)
#     # reshape column-major -> permutedims -> collect for Matrix
#     conn_matrix = collect(permutedims(reshape(raw_conn_vec, F, num_elements)))

#     # --- Read Families (FAM) ---
#     if haskey(group, "FAM")
#         fam_dset = group["FAM"]
#         fam_vec_raw = read(fam_dset)
#         # Ensure family IDs are read as the specified integer type I
#         fam_vec = convert(Vector{I}, vec(fam_vec_raw)) # Ensure it's a vector of type I
#         if length(fam_vec) != num_elements
#             error("Family ID ('FAM') vector length ($(length(fam_vec))) does not match element count ($num_elements) for '$elemType'.")
#         end
#     else
#         # If FAM dataset is missing, assign a default value (e.g., 0) or handle as needed.
#         # Using 0 as a placeholder for "no specific family".
#         fam_vec = zeros(I, num_elements)
#         @debug "Element block '$elemType' is missing the 'FAM' dataset. Assigning family ID 0 to all elements."
#     end

#     # --- Read Original Numbers (NUM) - Optional, not strictly needed by spec, but good practice ---
#     # original_element_ids = nothing
#     # if haskey(group, "NUM")
#     #     num_dset = group["NUM"]
#     #     original_element_ids = convert(Vector{I}, vec(read(num_dset)))
#     #     if length(original_element_ids) != num_elements
#     #         @warn "Original element ID ('NUM') vector length ($(length(original_element_ids))) does not match element count ($num_elements) for '$elemType'. Ignoring NUM."
#     #         original_element_ids = nothing
#     #     end
#     # end

#     return conn_matrix, fam_vec, num_elements
# end


# """
#     read_MED_HDF5(path::AbstractString;
#                   scale::Real=1.0,
#                   integer::Type{<:Integer}=Int64,
#                   float::Type{<:AbstractFloat}=Float64)
#         → points :: Vector{Point{float, SVector{3, float}}},
#           faces  :: Vector{Face{integer, Vector{integer}}},
#           cells  :: Vector{Cell_MED{integer, Vector{integer}}},
#           boundaries :: Vector{BoundaryElement{String, integer, Vector{integer}}}

# Reads mesh data from a MED-HDF5 file format specified by `path`.

# This function implements the reading logic based on common MED conventions,
# handling Fortran-ordered arrays, node ID remapping via NOE/NUM, and family/boundary groups.
# It aims to provide clean, raw lists of points, faces, cells, and boundary definitions
# suitable for downstream processing (e.g., by XCALibre.jl's builder).

# **Arguments:**

# *   `path`: Path to the MED-HDF5 file.
# *   `scale`: Optional scaling factor applied to all point coordinates (default: 1.0).
# *   `integer`: Optional integer type to use for indices and counts (default: `Int64`).
# *   `float`: Optional float type for coordinates (default: `Float64`).

# **Returns:**

# A tuple containing four vectors:
# 1.  `points`: A vector of `Point` structs, preserving the order from the file.
# 2.  `faces`: A vector of `Face` structs, containing remapped and potentially reordered node indices.
# 3.  `cells`: A vector of `Cell_MED` structs, containing remapped and potentially reordered node indices.
# 4.  `boundaries`: A vector of `BoundaryElement` structs, grouping faces by family ID and name, sorted alphabetically by name.

# **MED Conventions Handled:**

# *   Locates the first mesh and first step/iteration within `/ENS_MAA`. Issues warnings if multiple are found.
# *   Reads coordinates (`/NOE/COO`), handling both N×3 and 3×N layouts.
# *   Reads original node IDs (`/NOE/NUM`) if present and builds a mapping to the 1-based point indices. If `NUM` is absent, assumes identity mapping.
# *   Iterates through element blocks (`/MAI/<elemType>`).
# *   Reads connectivity (`NOD`), handling Fortran order (reshape F×E -> transpose E×F).
# *   Reads family IDs (`FAM`).
# *   Remaps node IDs in connectivity arrays using the `NOE/NUM` mapping. Reports errors if a node ID is not found in the map.
# *   Optionally reorders nodes within elements (e.g., HE8, QU4) to a consistent convention (UNV) if specified in `getNodePermutation`.
# *   Reads boundary names from `/FAS/<meshName>/ELEME/FAM_<id>_<name>`.
# *   Groups faces into `BoundaryElement` structs based on their family IDs. Handles missing family names.
# *   Does **not** perform coordinate deduplication.
# *   Does **not** delete or modify elements after reading; flags or errors are preferred during the read phase.

# **Potential Warnings/Errors:**

# *   Warning if multiple meshes or steps are found in the file (only the first is used).
# *   Warning if `/NOE/NUM` dataset is missing (identity mapping assumed).
# *   Warning if `/MAI/<elemType>/FAM` dataset is missing (family ID 0 assumed).
# *   Warning if `/FAS/<meshName>/ELEME` group or specific family names are missing.
# *   Warning for unknown element types.
# *   Error if mandatory datasets (`NOE/COO`, `MAI/*/NOD`) are missing.
# *   Error if array dimensions mismatch (e.g., `NOD` length vs. `F * NBR`, `FAM` length vs. element count).
# *   Error if a node ID found in `NOD` is not present in the `nodeIdToIdx` map.
# """
# function read_MED_HDF5(path::AbstractString;
#                        scale::Real=1.0,
#                        integer::Type{<:Integer}=Int64,
#                        float::Type{<:AbstractFloat}=Float64)

#     # Type aliases for convenience
#     I = integer
#     F = float
#     SV3F = SVector{3, F}
#     PointF = Point{F, SV3F}
#     FaceI = Face{I, Vector{I}}
#     CellI = Cell_MED{I, Vector{I}}
#     BoundaryElementSI = BoundaryElement{String, I, Vector{I}}

#     # Output vectors initialization
#     points = Vector{PointF}()
#     faces = Vector{FaceI}()
#     cells = Vector{CellI}()
#     boundaries = Vector{BoundaryElementSI}()

#     # Internal tracking
#     faces_by_family = Dict{I, Vector{I}}() # Maps family ID -> list of face indices
#     face_counter = 0
#     cell_counter = 0
#     boundary_counter = 0 # Counter for BoundaryElement index

#     h5open(path, "r") do fid
#         # --- 3.1 Open file & locate mesh ---
#         if !haskey(fid, "ENS_MAA")
#             error("HDF5 file is missing the '/ENS_MAA' group. Is this a valid MED file?")
#         end
#         ens_maa_group = fid["ENS_MAA"]
#         mesh_names = keys(ens_maa_group)
#         if isempty(mesh_names)
#             error("No meshes found under '/ENS_MAA'.")
#         end
#         if length(mesh_names) > 1
#             @warn "Multiple meshes found: $(join(mesh_names, ", ")). Using the first one: '$(first(mesh_names))'."
#         end
#         meshName = first(mesh_names)
#         mesh_group = ens_maa_group[meshName]

#         step_names = keys(mesh_group)
#         if isempty(step_names)
#             error("No steps/iterations found under '/ENS_MAA/$meshName'.")
#         end
#         if length(step_names) > 1
#             # MED steps often look like "NOEUD____00000000000000000001" or similar iteration counters
#             # Or sometimes time steps like "0.000000e+00"
#             # Choosing the first one is usually safe for single-state meshes.
#             @warn "Multiple steps/iterations found: $(join(step_names, ", ")). Using the first one: '$(first(step_names))'."
#         end
#         stepName = first(step_names) # Often something like 'NOEUD...' or 'INFO_MED' or time step
        
#         # Base path for mesh data for this step
#         basePath = "/ENS_MAA/$meshName/$stepName"
#         if !haskey(fid, basePath)
#              error("Could not access mesh step group at path '$basePath'")
#         end
#         step_group = fid[basePath]

#         # Check for essential groups: NOE (nodes) and MAI (elements/maillage)
#         if !haskey(step_group, "NOE")
#             error("Mesh step group '$basePath' is missing the 'NOE' (nodes) subgroup.")
#         end
#         if !haskey(step_group, "MAI")
#             error("Mesh step group '$basePath' is missing the 'MAI' (elements) subgroup.")
#         end
#         noe_group = step_group["NOE"]
#         mai_group = step_group["MAI"]

#         # --- 3.2 Load coordinates (no dedup!) ---
#         if !haskey(noe_group, "COO")
#             error("Node group '$basePath/NOE' is missing the 'COO' (coordinates) dataset.")
#         end
#         coords_raw = read(noe_group["COO"])

#         # Handle potential Fortran vs C array ordering and 2D/3D
#         local coords_matrix::Matrix{F}
#         if ndims(coords_raw) == 1
#             # Assume flattened 3xN or Nx3 - check length
#             if length(coords_raw) % 3 == 0
#                  num_nodes_inferred = div(length(coords_raw), 3)
#                  coords_matrix = convert(Matrix{F}, permutedims(reshape(coords_raw, 3, num_nodes_inferred)))
#                  @debug "Read 1D coordinate array, reshaped to ($num_nodes_inferred, 3)."
#              else
#                  error("1D Coordinate array length ($(length(coords_raw))) is not divisible by 3.")
#              end
#         elseif ndims(coords_raw) == 2
#             if size(coords_raw, 2) == 3 # N x 3 (C order or already transposed Fortran)
#                  coords_matrix = convert(Matrix{F}, coords_raw)
#                  @debug "Read 2D coordinate array ($(size(coords_matrix,1)) x 3)."
#             elseif size(coords_raw, 1) == 3 # 3 x N (Fortran order)
#                  coords_matrix = convert(Matrix{F}, permutedims(coords_raw))
#                  @debug "Read 2D coordinate array (3 x $(size(coords_matrix,1))), transposed."
#             elseif size(coords_raw, 2) == 2 # N x 2 (2D mesh)
#                  coords_2d = convert(Matrix{F}, coords_raw)
#                  num_nodes_2d = size(coords_2d, 1)
#                  coords_matrix = zeros(F, num_nodes_2d, 3)
#                  coords_matrix[:, 1:2] .= coords_2d
#                  @debug "Read 2D coordinate array ($(num_nodes_2d) x 2), padded Z to 0."
#             elseif size(coords_raw, 1) == 2 # 2 x N (2D mesh, Fortran order)
#                  coords_2d_t = convert(Matrix{F}, permutedims(coords_raw))
#                  num_nodes_2d = size(coords_2d_t, 1)
#                  coords_matrix = zeros(F, num_nodes_2d, 3)
#                  coords_matrix[:, 1:2] .= coords_2d_t
#                  @debug "Read 2D coordinate array (2 x $(num_nodes_2d)), transposed and padded Z to 0."
#             else
#                  error("Unsupported coordinate array shape: $(size(coords_raw)). Expected N x 3, 3 x N, N x 2, or 2 x N.")
#             end
#         else
#             error("Unsupported coordinate array dimensionality: $(ndims(coords_raw)). Expected 1 or 2.")
#         end

#         num_points = size(coords_matrix, 1)
#         sizehint!(points, num_points) # Pre-allocate memory

#         # Create Point objects
#         for i in 1:num_points
#             # Apply scaling during Point creation
#             push!(points, PointF(SV3F(coords_matrix[i, 1] * scale,
#                                       coords_matrix[i, 2] * scale,
#                                       coords_matrix[i, 3] * scale))) # Z is already padded if needed
#         end
#         @info "Read $num_points points."

#         # --- 3.3 Build the node-ID ↔ point-index map ---
#         local nodeIdToIdx::Dict{I, Int}
#         if haskey(noe_group, "NUM")
#             origIds_raw = read(noe_group["NUM"])
#             origIds = convert(Vector{I}, vec(origIds_raw)) # Ensure vector of type I
#             if length(origIds) != num_points
#                 error("'/NOE/NUM' dataset length ($(length(origIds))) does not match number of points ($num_points).")
#             end
#             # Check for duplicate original IDs, which would be invalid
#             if length(Set(origIds)) != num_points
#                  @warn "Duplicate node IDs found in '/NOE/NUM'. This might indicate an issue with the MED file."
#                  # Keep going for now, but the Dict will arbitrarily pick one index per duplicate ID
#             end
#             nodeIdToIdx = Dict{I, Int}(id => i for (i, id) in enumerate(origIds))
#             @info "Built node ID map from '/NOE/NUM'."
#         else
#             @warn "'/NOE/NUM' dataset not found. Assuming identity mapping (original node ID == point index)."
#             # Assume 1-based indexing matches the row index in COO
#             nodeIdToIdx = Dict{I, Int}(I(i) => i for i in 1:num_points)
#         end

#         # --- 3.4 Read every element block ---
#         element_types = keys(mai_group)
#         if isempty(element_types)
#             @warn "No element blocks found under '$basePath/MAI'. Mesh will be empty."
#         end

#         total_faces_read = 0
#         total_cells_read = 0

#         for elemType in element_types
#             elem_group = mai_group[elemType]
#             @info "Processing element block: $elemType"

#             nodes_per_elem = expectedNodesPerElement(elemType)
#             if isnothing(nodes_per_elem)
#                  @warn "Skipping element block '$elemType' due to unknown nodes per element."
#                  continue
#             end

#             # Read raw connectivity (original node IDs) and family IDs
#             conn_matrix_orig, fam_vec, num_elements = read_connectivity_and_families(elem_group, elemType, nodes_per_elem, I)

#             if num_elements == 0
#                 @info "Element block '$elemType' contains 0 elements. Skipping."
#                 continue
#             end

#             # Get node permutation if applicable
#             permutation = getNodePermutation(elemType)

#             is_face_type = isFaceElementType(elemType)

#             # Process each element in the block
#             for i in 1:num_elements
#                 elem_nodes_orig = conn_matrix_orig[i, :]
#                 elem_nodes_idx = Vector{I}(undef, nodes_per_elem)

#                 # Remap original node IDs to point indices
#                 valid_element = true
#                 for j in 1:nodes_per_elem
#                     orig_node_id = elem_nodes_orig[j]
#                     if !haskey(nodeIdToIdx, orig_node_id)
#                         @error "Node ID '$orig_node_id' (from '$elemType' element #$i, local node #$j) not found in the node ID map derived from '/NOE/NUM'. Cannot process this element."
#                         valid_element = false
#                         break # Stop processing this element
#                     end
#                     elem_nodes_idx[j] = nodeIdToIdx[orig_node_id]
#                 end

#                 if !valid_element
#                     # Decide action: skip element with warning, or error out?
#                     # Following prompt: "if a node ID is missing, abort or mark..." - Aborting is safer.
#                     error("Aborting due to missing node ID reference. Check element $i of type $elemType.")
#                     # Or, to skip:
#                     # @warn "Skipping element $i of type $elemType due to missing node ID reference(s)."
#                     # continue # Skip to the next element
#                 end

#                 # Apply node reordering if permutation exists
#                 if !isnothing(permutation)
#                     if length(permutation) != nodes_per_elem
#                          @warn "Permutation length mismatch for $elemType. Skipping reordering."
#                     else
#                          permuted_nodes = elem_nodes_idx[permutation]
#                          elem_nodes_idx = permuted_nodes
#                     end
#                 end

#                 family_id = fam_vec[i]

#                 # Add to faces or cells list
#                 if is_face_type
#                     face_counter += 1
#                     push!(faces, FaceI(I(face_counter), I(nodes_per_elem), elem_nodes_idx))
#                     total_faces_read += 1

#                     # Store face index by family ID for boundary creation
#                     if family_id != 0 # Assuming 0 means no specific family or internal face
#                          if !haskey(faces_by_family, family_id)
#                              faces_by_family[family_id] = Vector{I}()
#                          end
#                          push!(faces_by_family[family_id], I(face_counter))
#                      end
#                 else # Assume it's a cell type
#                     cell_counter += 1
#                     push!(cells, CellI(I(cell_counter), I(nodes_per_elem), elem_nodes_idx))
#                     total_cells_read += 1
#                     # Note: Cell family IDs are read but not typically used for boundaries in the same way.
#                     # Could store them if needed later: e.g., cell_family_map[cell_counter] = family_id
#                 end
#             end # end loop over elements in block
#             @info "Processed $num_elements elements of type $elemType."
#         end # end loop over element types

#         @info "Finished reading elements: $total_faces_read faces, $total_cells_read cells."

#         # --- 3.5 Group boundary faces ---
#         family_id_to_name = Dict{I, String}()
#         fas_path = "/FAS/$meshName/ELEME" # Standard path for element family names

#         if haskey(fid, fas_path)
#             fas_group = fid[fas_path]
#             family_names_count = 0
#             for group_name in keys(fas_group)
#                 # Expect names like FAM_<id>_<name> or potentially just FAM_<id>
#                 parts = split(group_name, '_')
#                 if length(parts) >= 2 && parts[1] == "FAM"
#                     try
#                         # Family IDs can be negative
#                         family_id = parse(I, parts[2])
#                         # Reconstruct name if it contained underscores
#                         family_name = length(parts) > 2 ? join(parts[3:end], "_") : "Family_$(family_id)"
#                         family_id_to_name[family_id] = family_name
#                         family_names_count += 1
#                     catch e
#                         if e isa ArgumentError
#                             @warn "Could not parse family ID from group name '$group_name' under '$fas_path'. Skipping."
#                         else
#                             rethrow(e)
#                         end
#                     end
#                 else
#                      @warn "Unexpected group format '$group_name' under '$fas_path'. Expected 'FAM_<id>_<name>'. Skipping."
#                 end
#             end
#             @info "Read $family_names_count family names from '$fas_path'."
#         else
#             @warn "Family names group '$fas_path' not found. Boundary elements will use generated names."
#         end

#         # Create BoundaryElement objects
#         sizehint!(boundaries, length(faces_by_family))
#         for (family_id, face_indices) in faces_by_family
#             if isempty(face_indices)
#                 # This shouldn't happen based on the logic above, but check just in case
#                  @warn "Family ID $family_id has an empty list of faces. Skipping boundary creation for it."
#                  continue
#             end

#             boundary_name = get(family_id_to_name, family_id, "Unnamed_Family_$(family_id)")
#             boundary_counter += 1
#             push!(boundaries, BoundaryElementSI(boundary_name, I(boundary_counter), face_indices))
#         end
#         @info "Created $boundary_counter boundary elements from $(length(faces_by_family)) unique non-zero family IDs."

#         # --- 3.6 Return (Sorting happens outside HDF5 block) ---

#     end # h5open

#     # Sort boundary elements by name for deterministic output
#     sort!(boundaries, by = x -> x.name)

#     return points, faces, cells, boundaries
# end
