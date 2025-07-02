# connectivity_check.jl - Mesh Connectivity Verification Functions

export check_boundary_face_count
export check_internal_face_owners
export check_cell_face_owners # New: Inverse of check_internal_face_owners
export check_node_cell_consistency
export check_boundary_consistency # New: Check mesh.boundary_cellsID
export check_range_consistency # New: Check ranges vs array lengths

using StaticArrays # If needed by Mesh types

"""
    check_boundary_face_count(mesh) -> Int

Calculates the total number of boundary faces based on the ranges defined
in `mesh.boundaries`. Assumes boundary faces are indexed contiguously from 1.
"""
function check_boundary_face_count(mesh)
    nbfaces = 0
    if isempty(mesh.boundaries)
        return 0
    end
    # Find the maximum endpoint of all boundary ranges
    try
        nbfaces = maximum(maximum(boundary.IDs_range) for boundary in mesh.boundaries if !isempty(boundary.IDs_range))
    catch e # Handles cases where all ranges might be empty or boundaries is empty
        if e isa ArgumentError # Thrown by maximum([])
             return 0
        else
             rethrow(e)
        end
    end
    return nbfaces
end

"""
    check_internal_face_owners(mesh) -> Tuple{Vector{Int}, Vector{Bool}}

Checks if the owner cells listed for each internal face correctly reference
that face in their `cell_faces` list.

Returns:
    - `faces_checked`: Array of internal face IDs that were checked.
    - `results`: Boolean array indicating pass/fail for each checked face.
"""
function check_internal_face_owners(mesh)
    (; cells, faces, cell_faces) = mesh
    IntType = eltype(cell_faces) # Assuming IntType can be inferred

    nbfaces = check_boundary_face_count(mesh)
    num_total_faces = length(faces)
    num_internal_faces = num_total_faces - nbfaces

    if num_internal_faces <= 0
        println("Internal Face Owners check: No internal faces found. Skipping.")
        return IntType[], Bool[]
    end

    test_range = (nbfaces + 1):num_total_faces # Test only internal faces
    results = falses(num_internal_faces) # Preallocate results
    faces_checked = Vector{IntType}(undef, num_internal_faces) # IDs of internal faces checked

    @inbounds for (i, face_ID) in enumerate(test_range)
        faces_checked[i] = face_ID
        face = faces[face_ID]
        ownerCells = face.ownerCells # Should be [owner1, owner2]

        # Basic check for valid owner cell IDs
        if ownerCells[1] <= 0 || ownerCells[1] > length(cells) ||
           ownerCells[2] <= 0 || ownerCells[2] > length(cells)
            # println("Face $face_ID has invalid owner cell IDs: $ownerCells")
            continue # Fails (results[i] is already false)
        end

        owner1 = ownerCells[1] # Cell ID of owner 1
        owner2 = ownerCells[2] # Cell ID of owner 2

        # Check owner 1
        owner1_face_range = cells[owner1].faces_range
        owner1_fIDs = view(cell_faces, owner1_face_range) # Use view for efficiency
        owner1_check = face_ID in owner1_fIDs

        # Check owner 2
        owner2_face_range = cells[owner2].faces_range
        owner2_fIDs = view(cell_faces, owner2_face_range)
        owner2_check = face_ID in owner2_fIDs

        if owner1_check && owner2_check
            results[i] = true
        else
            # Optional: Print detailed failure info
            # println("Face owners not consistent for face $face_ID: Owners $ownerCells. Check Owner1: $owner1_check, Check Owner2: $owner2_check")
        end
    end

    passed = sum(results)
    failed = num_internal_faces - passed
    println("Internal Face Owners check: Passed $passed, Failed $failed (out of $num_internal_faces internal faces)")
    return faces_checked, results
end


"""
    check_cell_face_owners(mesh) -> Tuple{Vector{Int}, Vector{Bool}}

Checks if the faces listed for each cell in `cell_faces` correctly
reference that cell as one of their owners. This is the inverse check
of `check_internal_face_owners`.

Returns:
    - `cells_checked`: Array of cell IDs that were checked.
    - `results`: Boolean array indicating pass/fail for each checked cell.
"""
function check_cell_face_owners(mesh)
    (; cells, faces, cell_faces) = mesh
    IntType = eltype(cell_faces)
    num_cells = length(cells)

    if num_cells == 0
        println("Cell Face Owners check: No cells found. Skipping.")
        return IntType[], Bool[]
    end

    results = trues(num_cells) # Assume pass initially
    cells_checked = Vector{IntType}(1:num_cells)

    @inbounds for cell_ID in 1:num_cells
        cell_face_range = cells[cell_ID].faces_range
        # Skip cells with no internal faces listed
        if isempty(cell_face_range)
            continue
        end

        cell_fIDs = view(cell_faces, cell_face_range)

        for face_ID in cell_fIDs
            # Check if face_ID is valid
            if face_ID <= 0 || face_ID > length(faces)
                # println("Cell $cell_ID lists invalid face ID $face_ID")
                results[cell_ID] = false
                break # Fail this cell
            end

            face = faces[face_ID]
            ownerCells = face.ownerCells

            # Check if the current cell_ID is one of the owners
            if !(cell_ID in ownerCells)
                # println("Cell $cell_ID lists face $face_ID, but face owners are $ownerCells")
                results[cell_ID] = false
                break # Fail this cell
            end
        end
    end

    passed = sum(results)
    failed = num_cells - passed
    println("Cell Face Owners check: Passed $passed, Failed $failed (out of $num_cells cells)")
    return cells_checked, results
end


"""
    check_node_cell_consistency(mesh) -> Tuple{Vector{Int}, Vector{Bool}}

Checks if the cells listed for each node (`node_cells` via `node.cells_range`)
correctly contain that node in their `cell_nodes` list.

Returns:
    - `nodes_checked`: Array of node IDs that were checked.
    - `results`: Boolean array indicating pass/fail for each checked node.
"""
function check_node_cell_consistency(mesh)
    (; nodes, cells, node_cells, cell_nodes) = mesh
    IntType = eltype(node_cells)
    num_nodes = length(nodes)

    if num_nodes == 0
        println("Node-Cell Consistency check: No nodes found. Skipping.")
        return IntType[], Bool[]
    end

    results = trues(num_nodes) # Assume pass initially
    nodes_checked = Vector{IntType}(1:num_nodes)

    @inbounds for node_ID in 1:num_nodes
        node = nodes[node_ID]
        node_cell_range = node.cells_range

        # Skip nodes not connected to any cells
        if isempty(node_cell_range)
            continue
        end

        # Get the list of cells connected to this node
        cells_for_node = view(node_cells, node_cell_range)

        for cell_ID in cells_for_node
            # Check if cell_ID is valid
            if cell_ID <= 0 || cell_ID > length(cells)
                # println("Node $node_ID lists invalid cell ID $cell_ID")
                results[node_ID] = false
                break # Fail this node
            end

            cell = cells[cell_ID]
            cell_node_range = cell.nodes_range
            cell_nodes_for_cell = view(cell_nodes, cell_node_range)

            # Check if the current node_ID is listed in the cell's nodes
            if !(node_ID in cell_nodes_for_cell)
                # println("Node $node_ID lists cell $cell_ID, but cell nodes are $cell_nodes_for_cell")
                results[node_ID] = false
                break # Fail this node
            end
        end
    end

    passed = sum(results)
    failed = num_nodes - passed
    println("Node-Cell Consistency check: Passed $passed, Failed $failed (out of $num_nodes nodes)")
    return nodes_checked, results
end


"""
    check_boundary_consistency(mesh) -> Bool

Checks if the owner cell listed for each boundary face in `mesh.faces`
matches the corresponding entry in `mesh.boundary_cellsID`.
"""
function check_boundary_consistency(mesh)
    (; faces, boundary_cellsID) = mesh
    nbfaces = check_boundary_face_count(mesh)

    if nbfaces == 0
        println("Boundary Consistency check: No boundary faces found. Skipping.")
        return true
    end

    if length(boundary_cellsID) != nbfaces
        println("Failed: Length of boundary_cellsID ($(length(boundary_cellsID))) does not match calculated number of boundary faces ($nbfaces).")
        return false
    end

    passed = true
    @inbounds for bface_ID in 1:nbfaces
        face = faces[bface_ID]
        # Boundary faces should have owner1 == owner2
        if face.ownerCells[1] != face.ownerCells[2]
            println("Failed: Boundary face $bface_ID has inconsistent owners: $(face.ownerCells).")
            passed = false
            # continue # Continue checking other faces
        end
        # Check if owner matches boundary_cellsID
        if face.ownerCells[1] != boundary_cellsID[bface_ID]
            println("Failed: Boundary face $bface_ID owner ($(face.ownerCells[1])) does not match boundary_cellsID ($(boundary_cellsID[bface_ID])).")
            passed = false
        end
    end

    if passed
        println("Boundary Consistency check: Passed.")
    else
        println("Boundary Consistency check: Failed.")
    end
    return passed
end


"""
    check_range_consistency(mesh) -> Bool

Checks if the sum of lengths of ranges matches the length of the corresponding flat array
for nodes, faces, and cells connectivity.
"""
function check_range_consistency(mesh)
    (; cells, faces, nodes, cell_nodes, cell_faces, node_cells, face_nodes) = mesh
    all_passed = true

    # Cell Nodes
    total_nodes_in_ranges = isempty(cells) ? 0 : sum(length(c.nodes_range) for c in cells)
    if total_nodes_in_ranges != length(cell_nodes)
        println("Failed: Cell Nodes range lengths ($(total_nodes_in_ranges)) do not sum to length of cell_nodes ($(length(cell_nodes))).")
        all_passed = false
    end

    # Cell Faces
    total_faces_in_ranges = isempty(cells) ? 0 : sum(length(c.faces_range) for c in cells)
    if total_faces_in_ranges != length(cell_faces)
        println("Failed: Cell Faces range lengths ($(total_faces_in_ranges)) do not sum to length of cell_faces ($(length(cell_faces))).")
        all_passed = false
    end

    # Face Nodes
    total_nodes_in_face_ranges = isempty(faces) ? 0 : sum(length(f.nodes_range) for f in faces)
    if total_nodes_in_face_ranges != length(face_nodes)
        println("Failed: Face Nodes range lengths ($(total_nodes_in_face_ranges)) do not sum to length of face_nodes ($(length(face_nodes))).")
        all_passed = false
    end

    # Node Cells
    total_cells_in_node_ranges = isempty(nodes) ? 0 : sum(length(n.cells_range) for n in nodes)
    if total_cells_in_node_ranges != length(node_cells)
        println("Failed: Node Cells range lengths ($(total_cells_in_node_ranges)) do not sum to length of node_cells ($(length(node_cells))).")
        all_passed = false
    end

    if all_passed
        println("Range Consistency check: Passed.")
    else
        println("Range Consistency check: Failed.")
    end
    return all_passed
end


"""
    run_all_connectivity_checks(mesh)

Runs all standard connectivity checks on the mesh object.
"""
function run_all_connectivity_checks(mesh)
    println("\n--- Running Mesh Connectivity Checks ---")
    passed_all = true
    check_boundary_face_count(mesh) # Just prints count implicitly via other checks
    _, _ = check_internal_face_owners(mesh)
    _, _ = check_cell_face_owners(mesh)
    _, _ = check_node_cell_consistency(mesh)
    passed_all &= check_boundary_consistency(mesh)
    passed_all &= check_range_consistency(mesh)
    println("----------------------------------------")
    if passed_all
         println("All critical connectivity checks passed.")
    else
         println("One or more critical connectivity checks failed.")
    end
    return passed_all
end

# --- Obsolete Functions (Commented out or removed) ---
# check_cell_face_nodes: Relied on intermediate data and mesh uniformity
# check_all_cell_faces: Relied on intermediate data and mesh uniformity
# check_boundary_faces: Relied on intermediate data from UNV builder