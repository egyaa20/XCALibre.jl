export Div
export div! 
export div_tensor!

# Define Divergence type and functionality

struct Div{VF<:VectorField,FVF<:FaceVectorField,F,M}
    vector::VF
    face_vector::FVF
    values::Vector{F}
    mesh::M
end
Adapt.@adapt_structure Div
Div(vector::VectorField) = begin
    mesh = vector.mesh
    face_vector = FaceVectorField(mesh)
    values = zeros(F, length(mesh.cells))
    Div(vector, face_vector, values, mesh)
end

# Divergence function definition

function div!(phi::ScalarField, psif::FaceVectorField, config)
    # Extract variables for function
    mesh = phi.mesh
    # backend = _get_backend(mesh)
    (; cells, cell_nsign, cell_faces, faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    # Retrieve user-selected float type
    F = _get_float(mesh)

    # Launch main calculation kernel
    ndrange = length(cells)
    kernel! = div_kernel!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, F, cell_faces, cell_nsign, faces, phi, psif)
    # KernelAbstractions.synchronize(backend)

    # Retrieve number of boundary faces
    nbfaces = length(mesh.boundary_cellsID)

    # Launch boundary faces contribution kernel
    ndrange = nbfaces
    kernel! = div_boundary_faces_contribution_kernel!(_setup(backend, workgroup, ndrange)...)
    kernel!(faces, cells, phi, psif)
    # KernelAbstractions.synchronize(backend)
end

# Divergence calculation kernel

@kernel inbounds=true function div_kernel!(cells::AbstractArray{Cell{TF,SV,UR}}, F, cell_faces, cell_nsign, faces, phi, psif) where {TF,SV,UR}
    i = @index(Global)
    
    @inbounds begin
        # Extract required fields from cells structure
        (; volume, faces_range) = cells[i]
        
        # Set work item scalar field value as zero
        # phi.values[i] = 0.0 #zero(TF)
        reduction = zero(TF)
        # Loop over faces to iterate work item scalar field value 
        for fi ∈ faces_range
            # Extract face ID and corresponding normal direction
            fID = cell_faces[fi]
            nsign = cell_nsign[fi]

            # Extract required fields from faces structure
            (; area, normal) = faces[fID]

            # Scalar field values calculation
            Sf = area*normal
            # Atomix.@atomic phi.values[i] += psif[fID]⋅Sf*nsign/volume
            reduction += psif[fID]⋅Sf*nsign
        end
        phi.values[i] = reduction/volume # divide only once
    end
end

# Boundary faces contribution kernel

@kernel function div_boundary_faces_contribution_kernel!(faces, cells, phi, psif)
    i = @index(Global)
    
    @inbounds begin
        # Retreive variables from work item boundary face
        cID = faces[i].ownerCells[1]
        volume = cells[cID].volume
        (; area, normal) = faces[i]

        # Boundary contribution calculation (boundary normals are correct by definition)
        Sf = area*normal
        Atomix.@atomic phi.values[cID] += psif[i]⋅Sf/volume
        # phi.values[cID] += psif[i]⋅Sf/volume
    end
end

# Divergence function definition - FaceScalarField

function div!(phi::ScalarField, psif::FaceScalarField, config)
    # Extract variables for function
    mesh = phi.mesh
    # backend = _get_backend(mesh)
    (; cells, cell_nsign, cell_faces, faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    # Retrieve user-selected float type
    F = _get_float(mesh)

    # Launch main calculation kernel
    ndrange = length(cells)
    kernel! = div_noS_kernel!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, F, cell_faces, cell_nsign, faces, phi, psif)
    # KernelAbstractions.synchronize(backend)

    # Retrieve number of boundary faces
    nbfaces = length(mesh.boundary_cellsID)

    # Launch boundary faces contribution kernel
    ndrange = nbfaces
    kernel! = div_noS_boundary_faces_contribution_kernel!(
        _setup(backend, workgroup, ndrange)...)
    kernel!(faces, cells, phi, psif)
    # KernelAbstractions.synchronize(backend)
end

# Divergence calculation kernel - FaceScalarField

@kernel function div_noS_kernel!(cells::AbstractArray{Cell{TF,SV,UR}}, F, cell_faces, cell_nsign, faces, phi, psif) where {TF,SV,UR}
    i = @index(Global)
    
    @inbounds begin
        # Extract required fields from cells structure
        (; volume, faces_range) = cells[i]
        
        # Set work item scalar field value as zero
        # phi.values[i] = 0.0 #zero(TF)
        reduction = zero(TF)
        # Loop over faces to iterate work item scalar field value 
        for fi ∈ faces_range
            # Extract face ID and corresponding normal direction
            fID = cell_faces[fi]
            nsign = cell_nsign[fi]

            # Extract required fields from faces structure
            (; area, normal) = faces[fID]

            # Atomix.@atomic phi.values[i] += psif[fID]⋅Sf*nsign/volume
            reduction += psif[fID]*nsign
        end
        phi.values[i] = reduction/volume # divide only once
    end
end

# Boundary faces contribution kernel

@kernel function div_noS_boundary_faces_contribution_kernel!(faces, cells, phi, psif)
    i = @index(Global)
    
    @inbounds begin
        # Retreive variables from work item boundary face
        cID = faces[i].ownerCells[1]
        volume = cells[cID].volume
        (; area, normal) = faces[i]

        # Boundary contribution calculation (boundary normals are correct by definition)
        # Sf = area*normal
        Atomix.@atomic phi.values[cID] += psif[i]/volume
        # phi.values[cID] += psif[i]⋅Sf/volume
    end
end








# Divergence of Tensor Field

function div_tensor!(vector::VectorField, face_tensor::FaceTensorField, config)
    mesh = vector.mesh
    (; cells, cell_nsign, cell_faces, faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    F = _get_float(mesh)

    ndrange = length(cells)
    kernel! = div_tensor_kernel!(_setup(backend, workgroup, ndrange)...)
    kernel!(cells, F, cell_faces, cell_nsign, faces, vector, face_tensor)

    nbfaces = length(mesh.boundary_cellsID)

    ndrange = nbfaces
    kernel! = div_tensor_boundary_faces_contribution_kernel!(_setup(backend, workgroup, ndrange)...)
    kernel!(faces, cells, vector, face_tensor)
end


@kernel inbounds=true function div_tensor_kernel!(cells::AbstractArray{Cell{TF,SV,UR}}, F, cell_faces, cell_nsign, faces, vector, face_tensor) where {TF,SV,UR}
    i = @index(Global)
    
    @inbounds begin
        (; volume, faces_range) = cells[i]
        
        reduction_x = zero(TF)
        reduction_y = zero(TF)
        reduction_z = zero(TF)
        
        for fi ∈ faces_range
            fID = cell_faces[fi]
            nsign = cell_nsign[fi]

            (; area, normal) = faces[fID]
            
            Tx_dot_n = face_tensor.xx[fID]*normal[1] + face_tensor.xy[fID]*normal[2] + face_tensor.xz[fID]*normal[3]
            Ty_dot_n = face_tensor.yx[fID]*normal[1] + face_tensor.yy[fID]*normal[2] + face_tensor.yz[fID]*normal[3]
            Tz_dot_n = face_tensor.zx[fID]*normal[1] + face_tensor.zy[fID]*normal[2] + face_tensor.zz[fID]*normal[3]

            reduction_x += (Tx_dot_n * area) * nsign
            reduction_y += (Ty_dot_n * area) * nsign
            reduction_z += (Tz_dot_n * area) * nsign
        end
        
        vector.x[i] = reduction_x / volume
        vector.y[i] = reduction_y / volume
        vector.z[i] = reduction_z / volume
    end
end

@kernel function div_tensor_boundary_faces_contribution_kernel!(faces, cells, vector, face_tensor)
    i = @index(Global)
    
    @inbounds begin
        cID = faces[i].ownerCells[1]
        volume = cells[cID].volume
        (; area, normal) = faces[i]

        Tx_dot_n = face_tensor.xx[i]*normal[1] + face_tensor.xy[i]*normal[2] + face_tensor.xz[i]*normal[3]
        Ty_dot_n = face_tensor.yx[i]*normal[1] + face_tensor.yy[i]*normal[2] + face_tensor.yz[i]*normal[3]
        Tz_dot_n = face_tensor.zx[i]*normal[1] + face_tensor.zy[i]*normal[2] + face_tensor.zz[i]*normal[3]

        Atomix.@atomic vector.x.values[cID] += (Tx_dot_n * area) / volume
        Atomix.@atomic vector.y.values[cID] += (Ty_dot_n * area) / volume
        Atomix.@atomic vector.z.values[cID] += (Tz_dot_n * area) / volume
    end
end