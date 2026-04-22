@inline function boundary_interpolation!(
    BC::InletOutlet, phif::FaceScalarField, phi, boundary_cellsID, time, fID)
    @inbounds begin
        cID = boundary_cellsID[fID]
        phif[fID] = phi[cID]
    end
    nothing
end

@inline function boundary_interpolation!(
    BC::InletOutlet, psif::FaceVectorField, psi, boundary_cellsID, time, fID)
    @inbounds begin
        cID = boundary_cellsID[fID]
        psif[fID] = psi[cID]
    end
    nothing
end
