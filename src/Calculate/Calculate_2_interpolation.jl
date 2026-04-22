# export correct_boundaries!
export interpolate!
export interpolate_harmonic!
export interpolate_upwind!
export interpolate_vanleer!
export interpolate_superbee!
export interpolate_vanalbada!
export interpolate_minmod!


# Temporary functions to extract boundary array
function to_cpu(boundaries::AbstractArray)
    return boundaries
end

# Function to copy from GPU to CPU
function to_cpu(boundaries::AbstractGPUArray)
    # Copy boundaries to CPU
    boundaries_cpu = Array{eltype(boundaries)}(undef, length(boundaries))
    KernelAbstractions.copyto!(CPU(), boundaries_cpu, boundaries)
    return boundaries_cpu
end

# Function to correct interpolation at boundaries (expands loop to reduce allocations)


# @generated function correct_boundaries!(phif, phi, BCs, time, config)
#     unpacked_BCs = []
#     for i ∈ 1:length(BCs.parameters)
#         unpack = quote
#             #KERNEL LAUNCH
#             adjust_boundary!(BCs[$i], phif, phi, boundaries, boundary_cellsID, time, backend, workgroup)
#         end
#         push!(unpacked_BCs, unpack)
#     end
#     quote
#     (; mesh) = phif
#     (; boundary_cellsID, boundaries) = mesh 
#     (; hardware) = config
#     (; backend, workgroup) = hardware
#     $(unpacked_BCs...) 
#     end
# end

## SCALAR INTERPOLATION

function interpolate!(phif::FaceScalarField, phi::Union{ScalarField, ConstantScalar}, config)
    # Extract values arrays from scalar fields 
    vals = phi.values
    fvals = phif.values

    # Extract faces from mesh
    mesh = phif.mesh
    (; cells, faces) = mesh

    # Launch interpolate kernel
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = interpolate_Scalar!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, cells, faces)
    # # KernelAbstractions.synchronize(backend)
end

@kernel function interpolate_Scalar!(fvals, vals, cells, faces)
    # Define index for thread
    i = @index(Global)

    @inbounds begin
        # Deconstruct faces to use weight and ownerCells in calculations
        face = faces[i]
        (; weight, ownerCells) = face

        # Calculate initial values based on index queried from ownerCells
        owner1 = ownerCells[1]
        owner2 = ownerCells[2]
        phi1 = vals[owner1]
        phi2 = vals[owner2]

        one_minus_weight = one(weight) - weight
        fvals[i] = weight*phi1 + one_minus_weight*phi2 # check weight is used correctly!
    end
end


## VAN LEER HELPERS

@inline function _vanleer_component(phi_U, dphi, gp)
    if abs(dphi) < eps(typeof(dphi))
        return phi_U
    else
        r = gp / dphi
        psi = (r + abs(r)) / (one(r) + abs(r))
        return phi_U + psi / 2 * dphi
    end
end

@inline function _vanleer_vector(psi_U, dpsi, grad_proj)
    px = _vanleer_component(psi_U[1], dpsi[1], grad_proj[1])
    py = _vanleer_component(psi_U[2], dpsi[2], grad_proj[2])
    pz = _vanleer_component(psi_U[3], dpsi[3], grad_proj[3])
    return SVector(px, py, pz)
end

## VAN LEER SCALAR INTERPOLATION

function interpolate_vanleer!(phif::FaceScalarField, phi::ScalarField, grad::Grad, mdotf, config)
    vals = phi.values
    fvals = phif.values

    mesh = phif.mesh
    (; faces) = mesh

    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_vanleer!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, grad, mdotf, faces)
end

@kernel function _interpolate_vanleer!(fvals, vals, grad, mdotf, faces)
    i = @index(Global)

    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face

        owner1 = ownerCells[1]  # [o] - normal owner
        owner2 = ownerCells[2]  # [n] - neighbour

        flux = mdotf[i]

        if owner1 == owner2
            fvals[i] = vals[owner1]
        else
            # Select upwind (U) and downwind (D) based on flux direction
            # e points from owner1 to owner2
            if flux >= 0.0
                uID = owner1
                dID = owner2
                grad_U = grad[uID]
            else
                uID = owner2
                dID = owner1
                grad_U = -grad[uID]  # flip so projection is along owner1→owner2
            end

            phi_U = vals[uID]
            phi_D = vals[dID]
            dphi = phi_D - phi_U

            grad_proj = 2 * (grad_U ⋅ e) * delta

            fvals[i] = _vanleer_component(phi_U, dphi, grad_proj)
        end
    end
end

## SUPERBEE HELPERS

# SuperBee limiter: ψ(r) = max(0, min(2r,1), min(r,2)); zero for r ≤ 0.
@inline function _superbee_component(phi_U, dphi, gp)
    if abs(dphi) < eps(typeof(dphi))
        return phi_U
    else
        r = gp / dphi
        if r <= zero(r)
            psi = zero(r)
        else
            psi = max(min(2*r, one(r)), min(r, 2*one(r)))
        end
        return phi_U + psi / 2 * dphi
    end
end

@inline function _superbee_vector(psi_U, dpsi, grad_proj)
    px = _superbee_component(psi_U[1], dpsi[1], grad_proj[1])
    py = _superbee_component(psi_U[2], dpsi[2], grad_proj[2])
    pz = _superbee_component(psi_U[3], dpsi[3], grad_proj[3])
    return SVector(px, py, pz)
end

## SUPERBEE SCALAR INTERPOLATION

function interpolate_superbee!(phif::FaceScalarField, phi::ScalarField, grad::Grad, mdotf, config)
    vals = phi.values
    fvals = phif.values

    mesh = phif.mesh
    (; faces) = mesh

    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_superbee!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, grad, mdotf, faces)
end

@kernel function _interpolate_superbee!(fvals, vals, grad, mdotf, faces)
    i = @index(Global)

    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face

        owner1 = ownerCells[1]
        owner2 = ownerCells[2]

        flux = mdotf[i]

        if owner1 == owner2
            fvals[i] = vals[owner1]
        else
            if flux >= 0.0
                uID = owner1
                dID = owner2
                grad_U = grad[uID]
            else
                uID = owner2
                dID = owner1
                grad_U = -grad[uID]
            end

            phi_U = vals[uID]
            phi_D = vals[dID]
            dphi = phi_D - phi_U

            grad_proj = 2 * (grad_U ⋅ e) * delta

            fvals[i] = _superbee_component(phi_U, dphi, grad_proj)
        end
    end
end

## VAN ALBADA HELPERS

# Van Albada limiter: ψ(r) = (r² + r) / (r² + 1); smooth, differentiable.
# Approaches van Leer for large r but has continuous derivatives → better for
# adjoint / implicit methods. Zero for r ≤ 0.
@inline function _vanalbada_component(phi_U, dphi, gp)
    if abs(dphi) < eps(typeof(dphi))
        return phi_U
    else
        r = gp / dphi
        psi = r <= zero(r) ? zero(r) : (r*r + r) / (r*r + one(r))
        return phi_U + psi / 2 * dphi
    end
end

@inline function _vanalbada_vector(psi_U, dpsi, grad_proj)
    px = _vanalbada_component(psi_U[1], dpsi[1], grad_proj[1])
    py = _vanalbada_component(psi_U[2], dpsi[2], grad_proj[2])
    pz = _vanalbada_component(psi_U[3], dpsi[3], grad_proj[3])
    return SVector(px, py, pz)
end

## VAN ALBADA SCALAR INTERPOLATION

function interpolate_vanalbada!(phif::FaceScalarField, phi::ScalarField, grad::Grad, mdotf, config)
    vals  = phi.values
    fvals = phif.values
    mesh  = phif.mesh
    (; faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_vanalbada!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, grad, mdotf, faces)
end

@kernel function _interpolate_vanalbada!(fvals, vals, grad, mdotf, faces)
    i = @index(Global)
    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face
        owner1 = ownerCells[1]
        owner2 = ownerCells[2]
        flux   = mdotf[i]
        if owner1 == owner2
            fvals[i] = vals[owner1]
        else
            if flux >= 0.0
                uID    = owner1
                dID    = owner2
                grad_U = grad[uID]
            else
                uID    = owner2
                dID    = owner1
                grad_U = -grad[uID]
            end
            phi_U     = vals[uID]
            phi_D     = vals[dID]
            dphi      = phi_D - phi_U
            grad_proj = 2 * (grad_U ⋅ e) * delta
            fvals[i]  = _vanalbada_component(phi_U, dphi, grad_proj)
        end
    end
end

## VAN ALBADA VECTOR INTERPOLATION

function interpolate_vanalbada!(psif::FaceVectorField, psi::VectorField, mdotf, config)
    (; mesh) = psif
    (; faces) = mesh
    interpolate!(psif, psi, config)
    ∇psi = Grad{Midpoint}(psi)
    green_gauss!(∇psi, psif, config)
    (; xx, xy, xz, yx, yy, yz, zx, zy, zz) = ∇psi.result
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_vanalbada_vector!(_setup(backend, workgroup, ndrange)...)
    kernel!(psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
end

@kernel function _interpolate_vanalbada_vector!(
    psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
    i = @index(Global)
    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face
        owner1 = ownerCells[1]
        owner2 = ownerCells[2]
        flux   = mdotf[i]
        if owner1 == owner2
            psif[i] = psi[owner1]
        else
            if flux >= 0.0
                uID = owner1; dID = owner2; se = e
            else
                uID = owner2; dID = owner1; se = -e
            end
            psi_U  = psi[uID]
            psi_D  = psi[dID]
            grad_x = SVector(xx[uID], xy[uID], xz[uID])
            grad_y = SVector(yx[uID], yy[uID], yz[uID])
            grad_z = SVector(zx[uID], zy[uID], zz[uID])
            gp_x   = 2 * (grad_x ⋅ se) * delta
            gp_y   = 2 * (grad_y ⋅ se) * delta
            gp_z   = 2 * (grad_z ⋅ se) * delta
            px = _vanalbada_component(psi_U[1], psi_D[1] - psi_U[1], gp_x)
            py = _vanalbada_component(psi_U[2], psi_D[2] - psi_U[2], gp_y)
            pz = _vanalbada_component(psi_U[3], psi_D[3] - psi_U[3], gp_z)
            psif[i] = SVector(px, py, pz)
        end
    end
end


## MINMOD HELPERS

# MinMod limiter: ψ(r) = max(0, min(1, r)); most diffusive TVD limiter.
# Gives clean monotone results at the cost of accuracy — good baseline / debug.
@inline function _minmod_component(phi_U, dphi, gp)
    if abs(dphi) < eps(typeof(dphi))
        return phi_U
    else
        r    = gp / dphi
        psi  = max(zero(r), min(one(r), r))
        return phi_U + psi / 2 * dphi
    end
end

@inline function _minmod_vector(psi_U, dpsi, grad_proj)
    px = _minmod_component(psi_U[1], dpsi[1], grad_proj[1])
    py = _minmod_component(psi_U[2], dpsi[2], grad_proj[2])
    pz = _minmod_component(psi_U[3], dpsi[3], grad_proj[3])
    return SVector(px, py, pz)
end

## MINMOD SCALAR INTERPOLATION

function interpolate_minmod!(phif::FaceScalarField, phi::ScalarField, grad::Grad, mdotf, config)
    vals  = phi.values
    fvals = phif.values
    mesh  = phif.mesh
    (; faces) = mesh
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_minmod!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, grad, mdotf, faces)
end

@kernel function _interpolate_minmod!(fvals, vals, grad, mdotf, faces)
    i = @index(Global)
    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face
        owner1 = ownerCells[1]
        owner2 = ownerCells[2]
        flux   = mdotf[i]
        if owner1 == owner2
            fvals[i] = vals[owner1]
        else
            if flux >= 0.0
                uID    = owner1
                dID    = owner2
                grad_U = grad[uID]
            else
                uID    = owner2
                dID    = owner1
                grad_U = -grad[uID]
            end
            phi_U     = vals[uID]
            phi_D     = vals[dID]
            dphi      = phi_D - phi_U
            grad_proj = 2 * (grad_U ⋅ e) * delta
            fvals[i]  = _minmod_component(phi_U, dphi, grad_proj)
        end
    end
end

## MINMOD VECTOR INTERPOLATION

function interpolate_minmod!(psif::FaceVectorField, psi::VectorField, mdotf, config)
    (; mesh) = psif
    (; faces) = mesh
    interpolate!(psif, psi, config)
    ∇psi = Grad{Midpoint}(psi)
    green_gauss!(∇psi, psif, config)
    (; xx, xy, xz, yx, yy, yz, zx, zy, zz) = ∇psi.result
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_minmod_vector!(_setup(backend, workgroup, ndrange)...)
    kernel!(psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
end

@kernel function _interpolate_minmod_vector!(
    psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
    i = @index(Global)
    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face
        owner1 = ownerCells[1]
        owner2 = ownerCells[2]
        flux   = mdotf[i]
        if owner1 == owner2
            psif[i] = psi[owner1]
        else
            if flux >= 0.0
                uID = owner1; dID = owner2; se = e
            else
                uID = owner2; dID = owner1; se = -e
            end
            psi_U  = psi[uID]
            psi_D  = psi[dID]
            grad_x = SVector(xx[uID], xy[uID], xz[uID])
            grad_y = SVector(yx[uID], yy[uID], yz[uID])
            grad_z = SVector(zx[uID], zy[uID], zz[uID])
            gp_x   = 2 * (grad_x ⋅ se) * delta
            gp_y   = 2 * (grad_y ⋅ se) * delta
            gp_z   = 2 * (grad_z ⋅ se) * delta
            px = _minmod_component(psi_U[1], psi_D[1] - psi_U[1], gp_x)
            py = _minmod_component(psi_U[2], psi_D[2] - psi_U[2], gp_y)
            pz = _minmod_component(psi_U[3], psi_D[3] - psi_U[3], gp_z)
            psif[i] = SVector(px, py, pz)
        end
    end
end


## UPWIND SCALAR INTERPOLATION

function interpolate_upwind!(phif::FaceScalarField, phi::ScalarField, mdotf, config)
    # Extract values arrays from scalar fields 
    vals = phi.values
    fvals = phif.values

    # Extract faces from mesh
    mesh = phif.mesh
    (; cells, faces) = mesh

    # Launch interpolate kernel
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_upwind!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, mdotf, faces)
    # # KernelAbstractions.synchronize(backend)
end

@kernel function _interpolate_upwind!(fvals, vals, mdotf, faces)
    i = @index(Global)

    @inbounds begin
        face = faces[i]
        (; weight, ownerCells, normal) = face
        F = face.centre

        owner1 = ownerCells[1] # [o]
        owner2 = ownerCells[2] # [n]

        flux = mdotf[i]

        if owner1 == owner2
            fvals[i] = vals[owner1]
        else
            if flux >= 0.0
                fvals[i] = vals[owner1]
            else
                fvals[i] = vals[owner2]
            end
        end
    end
end


## UPWIND VECTOR INTERPOLATION

function interpolate_upwind!(phif::FaceVectorField, phi::VectorField, mdotf, config)
    # Extract faces from mesh
    mesh = phif.mesh
    (; cells, faces) = mesh

    # Launch interpolate kernel
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_upwind_vec!(_setup(backend, workgroup, ndrange)...)
    kernel!(phif, phi, mdotf, faces)
    # # KernelAbstractions.synchronize(backend)
end

@kernel function _interpolate_upwind_vec!(phif, phi, mdotf, faces)
    i = @index(Global)

    @inbounds begin
        face = faces[i]
        (; weight, ownerCells, normal) = face
        # F = face.centre

        owner1 = ownerCells[1] # [o]
        owner2 = ownerCells[2] # [n]

        flux = mdotf[i]

        if owner1 == owner2 #unnecesary line!
            phif[i] = phi[owner1]
        else
            if flux >= 0.0
                phif[i] = phi[owner1]
            else
                phif[i] = phi[owner2]
            end
        end
    end
end

## HARMONIC SCALAR INTERPOLATION

function interpolate_harmonic!(phif::FaceScalarField, phi::ScalarField, config)
    # Extract values arrays from scalar fields 
    vals = phi.values
    fvals = phif.values

    # Extract faces from mesh
    mesh = phif.mesh
    (; cells, faces) = mesh

    # Launch interpolate kernel
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = interpolate_harmonic_Scalar!(_setup(backend, workgroup, ndrange)...)
    kernel!(fvals, vals, cells, faces)
    # # KernelAbstractions.synchronize(backend)
end

@kernel function interpolate_harmonic_Scalar!(fvals, vals, cells, faces)
    # Define index for thread
    i = @index(Global)

    @inbounds begin
        # Deconstruct faces to use weight and ownerCells in calculations
        face = faces[i]
        (; ownerCells) = face

        # Calculate initial values based on index queried from ownerCells
        owner1 = ownerCells[1]
        owner2 = ownerCells[2]
        phi1 = vals[owner1]
        phi2 = vals[owner2]

        fvals[i] = 2*((phi1*phi2)/(phi1+phi2))
    end
end

# VECTOR INTERPOLATION
function interpolate!(psif::FaceVectorField, psi::VectorField, config)
    # Extract x, y, z, values from FaceVectorField
    (; mesh) = psif

    #Redefine x, y, z values to be used in kernel
    xf = psif.x
    yf = psif.y
    zf = psif.z

    # Extract x, y, z, values from VectorField
    xv = psi.x
    yv = psi.y
    zv = psi.z

    #Extract faces array from mesh
    faces = mesh.faces

    # Launch interpolate kernel
    # backend = _get_backend(mesh)
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = interpolate_Vector!(_setup(backend, workgroup, ndrange)...)
    kernel!(xv, yv, zv, xf, yf, zf, faces)
    # # KernelAbstractions.synchronize(backend)
end

@kernel function interpolate_Vector!(@Const(xv), @Const(yv), @Const(zv), xf, yf, zf, @Const(faces))
    # Define index for thread
    i = @index(Global)

    @inbounds begin
        # Deconstruct faces to use weight and ownerCells in calculations
        # @synchronize # commented out on 2024/10/24
        (; weight, ownerCells) = faces[i]

        # Define indices for initial x and y values from psi struct
        cID1 = ownerCells[1]; cID2 = ownerCells[2]
        x1 = xv[cID1]; x2 = xv[cID2]
        y1 = yv[cID1]; y2 = yv[cID2]
        z1 = zv[cID1]; z2 = zv[cID2]

        # Calculate one minus weight
        one_minus_weight = one(weight) - weight

        # Update psif x and y arrays for interpolation (IMPLEMENT 3D)
        xf[i] = weight*x1 + one_minus_weight*x2 # check weight is used correctly!
        yf[i] = weight*y1 + one_minus_weight*y2 # check weight is used correctly!
        zf[i] = weight*z1 + one_minus_weight*z2 # check weight is used correctly!
    end
end

## VAN LEER VECTOR INTERPOLATION

function interpolate_vanleer!(psif::FaceVectorField, psi::VectorField, mdotf, config)
    (; mesh) = psif
    (; faces) = mesh

    # Interpolate psi onto faces (linear) for green gauss input
    interpolate!(psif, psi, config)

    # Compute gradient of vector field via green gauss → tensor field
    ∇psi = Grad{Midpoint}(psi)
    green_gauss!(∇psi, psif, config)

    # Extract tensor field components
    (; xx, xy, xz, yx, yy, yz, zx, zy, zz) = ∇psi.result

    # Launch van Leer kernel component by component
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_vanleer_vector!(_setup(backend, workgroup, ndrange)...)
    kernel!(psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
end

@kernel function _interpolate_vanleer_vector!(
    psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
    i = @index(Global)

    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face

        owner1 = ownerCells[1]
        owner2 = ownerCells[2]

        flux = mdotf[i]

        if owner1 == owner2
            psif[i] = psi[owner1]
        else
            if flux >= 0.0
                uID = owner1
                dID = owner2
                se = e
            else
                uID = owner2
                dID = owner1
                se = -e
            end

            psi_U = psi[uID]
            psi_D = psi[dID]

            # Gradient rows for upwind cell (∇ψ_x, ∇ψ_y, ∇ψ_z)
            grad_x = SVector(xx[uID], xy[uID], xz[uID])
            grad_y = SVector(yx[uID], yy[uID], yz[uID])
            grad_z = SVector(zx[uID], zy[uID], zz[uID])

            # Project each gradient row onto signed e and apply van Leer per component
            gp_x = 2 * (grad_x ⋅ se) * delta
            gp_y = 2 * (grad_y ⋅ se) * delta
            gp_z = 2 * (grad_z ⋅ se) * delta

            px = _vanleer_component(psi_U[1], psi_D[1] - psi_U[1], gp_x)
            py = _vanleer_component(psi_U[2], psi_D[2] - psi_U[2], gp_y)
            pz = _vanleer_component(psi_U[3], psi_D[3] - psi_U[3], gp_z)
            psif[i] = SVector(px, py, pz)
        end
    end
end

## SUPERBEE VECTOR INTERPOLATION

function interpolate_superbee!(psif::FaceVectorField, psi::VectorField, mdotf, config)
    (; mesh) = psif
    (; faces) = mesh

    interpolate!(psif, psi, config)

    ∇psi = Grad{Midpoint}(psi)
    green_gauss!(∇psi, psif, config)

    (; xx, xy, xz, yx, yy, yz, zx, zy, zz) = ∇psi.result

    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(faces)
    kernel! = _interpolate_superbee_vector!(_setup(backend, workgroup, ndrange)...)
    kernel!(psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
end

@kernel function _interpolate_superbee_vector!(
    psif, psi, xx, xy, xz, yx, yy, yz, zx, zy, zz, mdotf, faces)
    i = @index(Global)

    @inbounds begin
        face = faces[i]
        (; ownerCells, e, delta) = face

        owner1 = ownerCells[1]
        owner2 = ownerCells[2]

        flux = mdotf[i]

        if owner1 == owner2
            psif[i] = psi[owner1]
        else
            if flux >= 0.0
                uID = owner1
                dID = owner2
                se = e
            else
                uID = owner2
                dID = owner1
                se = -e
            end

            psi_U = psi[uID]
            psi_D = psi[dID]

            grad_x = SVector(xx[uID], xy[uID], xz[uID])
            grad_y = SVector(yx[uID], yy[uID], yz[uID])
            grad_z = SVector(zx[uID], zy[uID], zz[uID])

            gp_x = 2 * (grad_x ⋅ se) * delta
            gp_y = 2 * (grad_y ⋅ se) * delta
            gp_z = 2 * (grad_z ⋅ se) * delta

            px = _superbee_component(psi_U[1], psi_D[1] - psi_U[1], gp_x)
            py = _superbee_component(psi_U[2], psi_D[2] - psi_U[2], gp_y)
            pz = _superbee_component(psi_U[3], psi_D[3] - psi_U[3], gp_z)
            psif[i] = SVector(px, py, pz)
        end
    end
end


# GRADIENT INTERPOLATION

function interpolate!(
    gradf::FaceVectorField, grad::Grad, phi
    )
    (; mesh, x, y, z) = gradf
    (; cells, faces) = mesh
    (; values) = phi
    nbfaces = total_boundary_faces(mesh)
    start = nbfaces + 1
    @inbounds for fID ∈ start:length(faces)
        face = faces[fID]
        (; delta, ownerCells, e) = face
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        grad1 = grad(cID1)
        grad2 = grad(cID2)
        # get weight for current scheme
        w, df = weight(get_scheme(grad), cells, faces, fID)
        one_minus_weight = 1.0 - w
        # calculate interpolated value
        grad_ave = w*grad1 + one_minus_weight*grad2
        # correct interpolation
        grad_corr = grad_ave + ((values[cID2] - values[cID1])/delta - (grad_ave⋅e))*e
        x[fID] = grad_corr[1]
        y[fID] = grad_corr[2]
        z[fID] = grad_corr[3]
    end
end
