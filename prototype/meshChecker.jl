module meshCheck

using XCALibre
using LinearAlgebra
using Statistics #to compute mean values

export detect_element_type, detect_mesh_element_types, count_element_types, mesh_identifier, print_element_counts, 
        tetra_checker, volume_check_simple, volume_check_complex, nsign_check, e_value_check,
        compute_tetrahedron_volume, is_tetrahedron, general_mesh_check


function detect_element_type(cell, mesh)
    nodes_pointer = mesh.cell_nodes[cell.nodes_range]
    num_nodes = length(nodes_pointer)
    if num_nodes == 4
        return "tetrahedron"
    elseif num_nodes == 5
        return "pyramid"
    elseif num_nodes == 6
        return "prism"
    elseif num_nodes == 8
        return "hexahedron"
    else
        return "polyhedron"
    end
    #Q: maybe worth returning code-numbers rather than strings for efficiency? (applies to other functions too)
end

function detect_mesh_element_types(mesh)
    num_cells = length(mesh.cells)
    types = Vector{String}(undef, num_cells)
    for i in 1:num_cells
        types[i] = detect_element_type(mesh.cells[i], mesh)
    end
    return types
end

function count_element_types(mesh)
    types = detect_mesh_element_types(mesh)
    counts = Dict("tetrahedron" => 0,
                  "pyramid"     => 0,
                  "prism"       => 0,
                  "hexahedron"  => 0,
                  "polyhedron"  => 0)
    for t in types
        if haskey(counts, t)
            counts[t] += 1
        else
            counts[t] = 1
        end
    end
    return counts
end

function mesh_identifier(tetra, hexa, prisms, pyramids, polyhedra)
    if tetra > 0 && hexa == 0 && prisms == 0 && pyramids == 0 && polyhedra == 0
        return "[Mesh type] : Tetra"
    end
    if tetra > 0 && prisms > 0 && hexa == 0 && pyramids == 0 && polyhedra == 0
        return "[Mesh type] : Tetra with prismatic boundary layers"
    end
    
    if hexa > 0 && tetra == 0 && prisms == 0 && pyramids == 0 && polyhedra == 0
        return "[Mesh type] : Hexahedron"
    end
    
    if hexa > 0 && (tetra > 0 || prisms > 0 || pyramids > 0 || polyhedra > 0)
        return "[Mesh type] : Poly-hexcore"
    end
    
    return "mixed/unidentified mesh type"
end

function print_element_counts(mesh)
    counts = count_element_types(mesh)
    println("Element type counts:")
    for (element_type, count) in counts
        println("> $(element_type): $(count)")
    end

    tetra = get(counts, "tetrahedron", 0)
    hexa = get(counts, "hexahedron", 0)
    prisms = get(counts, "prism", 0)
    pyramids = get(counts, "pyramid", 0)
    polyhedra = get(counts, "polyhedron", 0)

    mesh_type = mesh_identifier(tetra, hexa, prisms, pyramids, polyhedra)
    println(mesh_type)
end























# PREVIOUS UPDATES:

# SIMPLE TETRA CHECKER:

function tetra_checker(mesh)
    cell = mesh.cells[1]

    faces_pointer = cell.faces_range
    nodes_array = mesh.cell_nodes[faces_pointer]
    
    vertex_ids = nodes_array
    cell_coords = [mesh.nodes[ID].coords for ID in vertex_ids]

    result = is_tetrahedron(cell_coords, mesh) # if true then it is tetrahedra element
    return result

    # would not work if it is a mixed mesh (e.g. tetra with prisms)
    # improvement: loop over all cells and sort tetra elements and other types of elements?
end



function volume_check_simple(expected_volume, mesh) #sums individual volumes by taking the volume value from mesh data
    sum_of_volumes = 0.0

    for cell ∈ mesh.cells
        sum_of_volumes += cell.volume
    end

    error = 100*(expected_volume-sum_of_volumes)/expected_volume
    return sum_of_volumes, error
end




function volume_check_complex(expected_volume::Real, mesh) # works for either polyhedral or arbitrary tetra meshes
    total_volume = 0.0
    is_tetra     = tetra_checker(mesh)     # your existing test

    ########################################################################
    # Fast path ─ tetrahedral mesh
    ########################################################################
    if is_tetra
        for cell in mesh.cells
            node_ids = mesh.cell_nodes[cell.nodes_range]              # 4 nodes
            coords   = [mesh.nodes[id].coords for id in node_ids]
            total_volume += compute_tetrahedron_volume(coords, mesh)
        end

    ########################################################################
    # General path ─ arbitrary polyhedral mesh
    ########################################################################
    else
        for cell in mesh.cells
            centroid   = cell.centre
            cell_vol   = 0.0
            face_ids   = mesh.cell_faces[cell.faces_range]

            for fID in face_ids
                node_ids = mesh.face_nodes[mesh.faces[fID].nodes_range]
                verts    = [mesh.nodes[id].coords for id in node_ids]
                n        = length(verts)
                n ≥ 3 || error("Face $fID has only $n vertices")

                # fan‑triangulate the face (centroid stays the same)
                v1 = verts[1]
                for i in 2:(n-1)
                    v2 = verts[i]
                    v3 = verts[i+1]
                    tetra = (centroid, v1, v2, v3)
                    cell_vol += compute_tetrahedron_volume(tetra, mesh)
                end
            end

            total_volume += cell_vol
        end
    end

    rel_error_pct = 100 * (expected_volume - total_volume) / expected_volume
    return total_volume, rel_error_pct
end













# function volume_check_complex(expected_volume, mesh) 
# #decompises non-tetrahedra mesh into tetrahedra if required => computes cells' volume using faces area data from mesh => sums it up

# # Loop over all cells
#     # Extract centroid coordinates
#     # Loop over all faces
#         # Split quadrilateral face into triangles
#         # Connect newly created triangular surface's vertices to the centroid
#         # Compute the volume for newly created tetrahedra and save it (using compute_tetrahedron_volume function)
#     # Compare the volume of tetrahedras to the volume of polyhedra element


#     total_volume = 0.0

#     for cell ∈ mesh.cells

#         centroid = cell.centre
#         cell_volume = 0.0
#         faces_pointer = cell.faces_range
#         faces_array = mesh.cell_faces[faces_pointer]

#         is_mesh_tetrahedral = tetra_checker(mesh)

#         #if mesh is non tetra, then decompose it into tetra
#         #if it is tetra, compute volumes straight away

#         if (is_mesh_tetrahedral)
#             for face ∈ faces_array
#                 nodes_pointer = mesh.faces[face].nodes_range
#                 nodes_array = mesh.face_nodes[nodes_pointer]
    
                
#                 nodes = [mesh.nodes[ID] for ID in nodes_array]

#                 tetra = [nodes[1], nodes[2], nodes[3], nodes[4]]
#                 vol = compute_tetrahedron_volume(tetra, mesh)

#                 total_volume += vol
#             end
#         else
#             for face ∈ faces_array
#                 nodes_pointer = mesh.faces[face].nodes_range
#                 nodes_array = mesh.face_nodes[nodes_pointer]
    

#                 nodes = [mesh.nodes[ID] for ID in nodes_array]
#                 n = length(nodes)
    
#                 if n < 3
#                     error("Face has less than 3 nodes")
#                 end
    
#                 v1 = nodes[1].coords #define vertex
#                 for i in 1:(n-2)
#                     v2 = nodes[i+1].coords
#                     v3 = nodes[i+2].coords
                    
#                     tetra = [centroid, v1, v2, v3]
                    
#                     vol = compute_tetrahedron_volume(tetra, mesh)
#                     cell_volume += vol
#                 end
#             end
    
#             total_volume += cell_volume # add tetra element volume to the total calculation
#         end
#     end

#     error = 100*(expected_volume-total_volume)/expected_volume
#     return total_volume, error
# end



function nsign_check(mesh)

    for (face_ID, face) in enumerate(mesh.faces)

        cell1_ID = face.ownerCells[1]
        cell2_ID = face.ownerCells[2]


        boundary_face = false

        if cell1_ID == cell2_ID
            boundary_face = true
            cell = mesh.cells[cell1_ID]

            # the difference is that for boundary faces, mesh.cell_faces will only return internal faces (so findfirst() is useless)
            # do we even need to check those?
            continue
        else
            cell1_faces_range = mesh.cells[cell1_ID].faces_range
            cell2_faces_range = mesh.cells[cell2_ID].faces_range #pointer
            
            face_IDs_1 = mesh.cell_faces[cell1_faces_range]
            face_IDs_2 = mesh.cell_faces[cell2_faces_range]

            face_index_1 = findfirst(==(face_ID), face_IDs_1) #looks for index in array
            face_index_2 = findfirst(==(face_ID), face_IDs_2) #looks for index in array

            actual_nsign_1 = mesh.cell_nsign[cell1_faces_range][face_index_1]
            actual_nsign_2 = mesh.cell_nsign[cell2_faces_range][face_index_2]

            if cell1_ID > cell2_ID
                if actual_nsign_1 != -1
                    error("Expected nsign to be -1 since cell1_ID ($cell1_ID) > cell2_ID ($cell2_ID), but got $actual_nsign_1")
                end
            elseif cell2_ID > cell1_ID
                if actual_nsign_2 != -1
                    error("Expected nsign to be -1 since cell2_ID ($cell2_ID) > cell1_ID ($cell1_ID), but got $actual_nsign_2")
                end
            end
        end
    end

    return true
end


function e_value_check(mesh) 
#takes coordinates of two centroids => computes unit vector between them => compares with e value from the mesh

    for (face_id, face) in enumerate(mesh.faces)
        face_id = 40
        face = mesh.faces[face_id]
        cell1_ID = face.ownerCells[1]
        cell2_ID = face.ownerCells[2]

        boundary_face = false

        if cell1_ID == cell2_ID
            boundary_face = true
            cell = mesh.cells[cell1_ID]
            
            # e is expected to point outwards of the domain e.g. be the same as normal
            # check if dot product equals to one?

            result = dot(face.e, face.normal)

            if result != 1.0
                error("Unit vector e at boundary face is not in the direction of face's normal")
            end
            
        else
            cell1 = mesh.cells[cell1_ID]
            cell2 = mesh.cells[cell2_ID]

            centroid1 = cell1.centre
            centroid2 = cell2.centre

            if cell1_ID < cell2_ID
                direction = centroid2 .- centroid1
            else
                direction = centroid1 .- centroid2
            end

            magnitude = norm(direction)
            unit_vector = direction / magnitude
            if (unit_vector != face.e)
                error("Unit vector e is not properly calculated")
            end
        end
    end

    return true
end


function compute_tetrahedron_volume(cell_vertices, mesh) #self-explanitory
    A, B, C, D = cell_vertices

    AB = B .- A
    AC = C .- A
    AD = D .- A

    cross_product = cross(AC, AD)
    tripple_product = dot(AB, cross_product)
    volume = abs(tripple_product) / 6.0

    return volume
end

function is_tetrahedron(cell_vertices, mesh) #logic: check if there are 4 vertices and that the volume can be computed (found this method online)
    if length(cell_vertices) != 4
        return false
    else
        volume = compute_tetrahedron_volume(cell_vertices, mesh)

        epsilon = 1e-10 #floating-point error threshold
        return volume > epsilon #make sure we avoid floating-point errors
    end
    
end



function general_mesh_check(expected_volume, mesh)

    # #Stage 1: tetra check

    # println("[Checking if mesh is tetrahedral....]")
    # is_mesh_tetrahedral = tetra_checker(mesh)
    # if is_mesh_tetrahedral
    #     println("Mesh is tetrahedral")
    # else
    #     println("Mesh is not tetrahedral... Decomposition might be required for some checks")
    # end

    # #Stage 2: simple volume check based on cell volume data recorded in mesh

    # println("\n[Running simple volume inspection...]")
    # simple_volume_result, simple_volume_error = volume_check_simple(expected_volume, mesh)
    # println("Simple volume inspection completed.\nExpected volume: $expected_volume\nComputed volume: $simple_volume_result\nError: $simple_volume_error %")
    
    # #Stage 3: complex volume check (split mesh into tetra if required and calculate individual volumes based on geometry)
    # println("\n[Running complex volume inspection...]")
    # complex_volume_result, complex_volume_error = volume_check_complex(expected_volume, mesh)
    # println("Complex volume inspection completed.\nExpected volume: $expected_volume\nComputed volume: $complex_volume_result\nError: $complex_volume_error %")
   
    # #Stage 4: n_signs check
    # println("\n[Running n signs inspection...]")
    # nsign_inspection_state = nsign_check(mesh)
    # if nsign_inspection_state
    #     println("n sign inspection completed. Everything is correct.")
    # end

    # #Stage 5: e vectors check
    # println("\n[Running e vectors inspection...]")
    # e_inspection_state = e_value_check(mesh)
    # if e_inspection_state
    #     println("e vectors inspection completed. Everything is correct.")
    # end
end


end
