@inline function boundary_interpolation!(
    BC::totalPressure, phif::FaceScalarField, phi, boundary_cellsID, time, fID, rho_field, U_field)

    rhof = rho_field[fID]
    Uf = U_field[fID]
    Uf_mag = norm(Uf)

    @inbounds phif[fID] = BC.value - 0.5*rhof*Uf_mag^2
    nothing
end

@inline function boundary_interpolation!(
    BC::totalPressure, psif::FaceVectorField, psi, boundary_cellsID, time, fID, rho_field, U_field)

    rhof = rho_field[fID]
    Uf = U_field[fID]
    Uf_mag = norm(Uf)

    @inbounds psif[fID] = BC.value - 0.5*rhof*Uf_mag^2
    nothing
end







# Need to change this too.