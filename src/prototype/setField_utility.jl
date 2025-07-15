function setField!(mesh, field, value::Real, region_box::String)
    coords = map(s -> parse(Float64, s), split(replace(region_box, r"[()]" => ""))) #split input region into array of coords
    p1 = coords[1:3] #min corner of the box
    p2 = coords[4:6] #max corner of the box
    
    min_corner = min.(p1, p2) #account for any order of input coords
    max_corner = max.(p1, p2) #account for any order of input coords

    cells_in_region = Int[]
    
    for (id, cell) in enumerate(mesh.cells) #check that X, Y, Z coords are within the box
        center = cell.centre
        if (min_corner[1] <= center[1] <= max_corner[1] &&
            min_corner[2] <= center[2] <= max_corner[2] &&
            min_corner[3] <= center[3] <= max_corner[3])
            
            push!(cells_in_region, id)
            # Warning: if outer cell boundary is outside the box but its centre is within the box, this cell will count too
        end
    end

    if !isempty(cells_in_region)
        field.values[cells_in_region] .= value # Reassign values inside the field
    end
    
    return length(cells_in_region)
end