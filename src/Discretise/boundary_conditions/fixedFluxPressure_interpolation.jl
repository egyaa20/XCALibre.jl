@inline function boundary_interpolation!(
    BC::FixedFluxPressure, phif::FaceScalarField, phi, boundary_cellsID, time, fID)
    @inbounds begin
        (; faces) = phi.mesh
        face = faces[fID]
        (; delta) = face
        cID = boundary_cellsID[fID]
        local_i = fID - BC.IDs_range.start + 1
        phif[fID] = phi[cID] + delta*BC.value[local_i]
    end
    nothing
end
