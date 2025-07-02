module meshQualitySkew

using XCALibre
using LinearAlgebra
using Statistics

export calculate_skewness, compute_skewness_stats

function calculate_skewness(cell, cell_id, mesh)
    cell_centre = mesh.cells[cell_id].centre
    max_skeweness = 0.0

    for pointer in cell.faces_range
        face_id = mesh.cell_faces[pointer]
        face = mesh.faces[face_id]
        face_centre = face.centre
        
        n = face.normal

        if face.ownerCells[1] != face.ownerCells[2] # if the face is not a boundary face
            if cell_id == face.ownerCells[1]
                neighbour_id = face.ownerCells[2]
            else
                neighbour_id = face.ownerCells[1]
            end
            neighbour_cell_centre = mesh.cells[neighbour_id].centre
        else
            continue #skip boundary faces
        end

        if cell_id > neighbour_id #account for nsign
            vector_between_cells = cell_centre - neighbour_cell_centre
        elseif cell_id < neighbour_id
            vector_between_cells = neighbour_cell_centre - cell_centre
        end
        
        denominator = dot(n, vector_between_cells)
        if abs(denominator) < eps() #basically check if denominator is zero
            continue
        end
        t = dot(n, face_centre - cell_centre) / denominator
        intersection_point = cell_centre + t * vector_between_cells

        skew = norm(face_centre - intersection_point) / norm(vector_between_cells)
        max_skeweness = max(max_skeweness, skew)
    end

    return max_skeweness
end


function compute_skewness_stats(mesh)
    num_cells = length(mesh.cells)
    values = zeros(Float64, num_cells)

    for (cell_id, cell) in enumerate(mesh.cells)
        values[cell_id] = calculate_skewness(cell, cell_id, mesh)
    end

    return (minimum(values), maximum(values), mean(values))
end

end
