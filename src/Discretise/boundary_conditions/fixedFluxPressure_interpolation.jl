@inline function boundary_interpolation!(
    BC::fixedFluxPressure, phif::FaceScalarField, phi, boundary_cellsID, time, fID)
    @inbounds begin
        (; faces) = phi.mesh
        face = faces[fID]
        (; delta) = face
        cID = boundary_cellsID[fID]
        phif[fID] = phi[cID] + delta*BC.value 
    end
    nothing
end


@inline function boundary_interpolation!(
    BC::fixedFluxPressure, psif::FaceVectorField, psi, boundary_cellsID, time, fID)
    @inbounds begin
        error("fixedFluxPressure boundary condition for vector fields is not implemented yet.")
    end
    nothing
end