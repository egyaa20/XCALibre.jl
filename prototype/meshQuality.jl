module meshQuality

include("meshSkewenessTest.jl")

using XCALibre
using LinearAlgebra
using Statistics #to compute mean values
using .meshQualitySkew

export compute_mesh_quality_stats, calculate_aspect_ratio, calculate_nonorthogonality,
            compute_aspect_ratio_stats, compute_nonorthogonality_stats,
            compute_skewness_stats, face_area_vector, is_degenerate_face, detect_degenerate_faces,
            cell_volume, detect_degenerate_volumes, report_mesh_quality



function calculate_aspect_ratio(cell, mesh)
#We compute a simple aspect ratio as the ratio of the longest edge length to the shortest

    nodes_pointer = mesh.cell_nodes[cell.nodes_range]
    coords = mesh.nodes[nodes_pointer]

    max_length = 0.0
    min_length = Inf
    n = length(coords)
    
    for i in 1:(n-1)
        for j in (i+1):n
#we check only unique pairs: if i=1 then j=2,3,4..... if i=2, then j=3,4 etc... 
#i=n is skipped because the pair "1<->n" is covered when i=1... no duplicates
            edge_length = norm(coords[j].coords - coords[i].coords)
            max_length = max(max_length, edge_length) #if there is longer edge => update
            min_length = min(min_length, edge_length) #if there is shorter edge => update
        end
    end
    return max_length / min_length
end

function calculate_nonorthogonality(cell, cell_id, mesh)
    max_angle = 0.0
    
    for pointer in cell.faces_range
        face_id = mesh.cell_faces[pointer]
        face = mesh.faces[face_id]

        if face.ownerCells[1] == face.ownerCells[2] #skip boundary faces
            continue
        end

        #Determine the neighbour
        if cell_id == face.ownerCells[1]
            neighbor_id = face.ownerCells[2]
        else
            neighbor_id = face.ownerCells[1]
        end
        neighbor_cell = mesh.cells[neighbor_id]


        if cell_id > neighbor_id #account for nsign
            vector_between_cells = cell.centre - neighbor_cell.centre
        elseif cell_id < neighbor_id
            vector_between_cells = neighbor_cell.centre - cell.centre
        end


        dot_product = dot(normalize(vector_between_cells), normalize(face.normal))
        #we compare the vector between cell centres and the vector normal to the face between cells
        #dot product result would be the cosine of the angle between them

        angle_rad = acos(clamp(dot_product, -1, 1)) #take the inverse of cosine
        angle_deg = rad2deg(angle_rad)
        max_angle = max(max_angle, angle_deg) #if there is a greater angle => update
    end
    return max_angle

    #comment : I considered to only compute angle for faces that have +ve nsign in order to not compute the same face two times...
    #... but I think it would break the logic because to find max angle we actually need to find max_angle for each individual cell, right?
end



function compute_aspect_ratio_stats(mesh)
    num_cells = length(mesh.cells)

    cell_AR_values = zeros(Float64, num_cells)

    for (cell_id, cell) in enumerate(mesh.cells)
        cell_AR_values[cell_id] = calculate_aspect_ratio(cell, mesh)
    end


    min = minimum(cell_AR_values)
    max = maximum(cell_AR_values)
    avg = mean(cell_AR_values)

    return (min, max, avg)
end

function compute_nonorthogonality_stats(mesh)
    num_cells = length(mesh.cells)

    cell_nonortho_values = zeros(Float64, num_cells)

    for (cell_id, cell) in enumerate(mesh.cells)
        cell_nonortho_values[cell_id] = calculate_nonorthogonality(cell, cell_id, mesh)
    end


    min = minimum(cell_nonortho_values)
    max = maximum(cell_nonortho_values)
    avg = mean(cell_nonortho_values)

    return (min, max, avg)
end



function compute_mesh_quality_stats(mesh)
    (ar_min, ar_max, ar_avg) = compute_aspect_ratio_stats(mesh)
    
    (nonortho_min, nonortho_max, nonortho_avg) = compute_nonorthogonality_stats(mesh)

    (sk_min, sk_max, sk_avg) = compute_skewness_stats(mesh)

    
    println("\nMesh Quality Statistics:")
    println("------------------------")
    println("Aspect Ratio:")
    println("   Minimum: $(round(ar_min, digits=3))")
    println("   Maximum: $(round(ar_max, digits=3))")
    println("   Average: $(round(ar_avg, digits=3))")

    println("\nNon-orthogonality:")
    println("   Minimum: $(round(nonortho_min, digits=3))°")
    println("   Maximum: $(round(nonortho_max, digits=3))°") 
    println("   Average: $(round(nonortho_avg, digits=3))°")

    println("\nSkewness:")
    println("   Minimum: $(round(sk_min, digits=3))")
    println("   Maximum: $(round(sk_max, digits=3))")
    println("   Average: $(round(sk_avg, digits=3))")
    
    # return (
    #     aspect_ratio = (min=ar_min, max=ar_max, avg=ar_avg),
    #     nonorthogonality = (min=nonortho_min, max=nonortho_max, avg=nonortho_avg)
    # )
end






function face_area_vector(mesh::Mesh3, fID::Int)
    # gather vertex coordinates
    node_ids = mesh.face_nodes[mesh.faces[fID].nodes_range]
    coords = [mesh.nodes[i].coords for i in node_ids]
    A = zero(typeof(coords[1]))
    n = length(coords)
    # sum cross products around the loop
    for j in 1:n
        vj = coords[j]
        vk = coords[ j % n + 1 ]
        A += cross(vj, vk)
    end
    return A / 2
end


function is_degenerate_face(mesh::Mesh3, fID::Int; length_tol=1e-12, area_tol=1e-24)
    node_ids = mesh.face_nodes[mesh.faces[fID].nodes_range]
    coords = [mesh.nodes[i].coords for i in node_ids]
    # coincident vertices
    if length(unique(coords)) < length(coords)
        return true
    end
    # tiny area
    area_vec = face_area_vector(mesh, fID)
    if norm(area_vec) < area_tol
        return true
    end
    return false
end

function detect_degenerate_faces(mesh::Mesh3; length_tol=1e-12, area_tol=1e-24)
    return [fID for fID in eachindex(mesh.faces) if is_degenerate_face(mesh, fID; length_tol=length_tol, area_tol=area_tol)]
end


function cell_volume(mesh::Mesh3, cID::Int)
    cell     = mesh.cells[cID]
    face_ids = mesh.cell_faces[cell.faces_range]
    signs    = mesh.cell_nsign[cell.faces_range]

    V = 0.0
    c = cell.centre                       # centroid inside the cell
    for (j, fID) in enumerate(face_ids)
        face  = mesh.faces[fID]
        A_vec = face.normal * face.area * signs[j]   # outward for this cell
        r_rel = face.centre - c                      # local position vector
        V    += dot(A_vec, r_rel)
    end
    return V / 3            # always positive; no abs needed
end



function detect_degenerate_volumes(mesh::Mesh3; vol_tol=1e-18)
    degenerate = Int[]               # empty result vector
    for cID in eachindex(mesh.cells) # loop over 1:length(mesh.cells)
        V = cell_volume(mesh, cID)
        if V < vol_tol
            push!(degenerate, cID)
        end
    end
    return degenerate
end


function report_mesh_quality(mesh::Mesh3; vol_tol=1e-18, length_tol=1e-12, area_tol=1e-24)
    # 1) Degenerate faces
    degen_faces = detect_degenerate_faces(mesh; length_tol=1e-12, area_tol=1e-24)
    # println("Degenerate faces (IDs): ", degen_faces)
    println("Number of degenerate faces: ", length(degen_faces))

    # 2) Degenerate cells
    degen_cells = detect_degenerate_volumes(mesh; vol_tol=1e-18)
    # println("Degenerate cells (IDs): ", degen_cells)
    println("Number of degenerate cells: ", length(degen_cells))

    vols = [ cell_volume(mesh, cID) for cID in eachindex(mesh.cells) ]
    vmin = minimum(vols)
    vmax = maximum(vols)
    vavg = mean(vols)

    println("Cell volume stats:")
    println("  min = ", vmin)
    println("  max = ", vmax)
    println("  avg = ", vavg)
    println("  sum = ", sum(vols))

end

end
