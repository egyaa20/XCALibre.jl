export totalPressure

"""
    totalPressure <: AbstractDirichlet

totalPressure boundary condition model.

# Inputs
- `ID` Name of the boundary given as a symbol (e.g. :inlet). Internally it gets replaced with the boundary index ID
- `value` Scalar p0 value for totalPressure boundary condition
"""
struct totalPressure{I,V,R<:UnitRange} <: AbstractDirichlet
    ID::I 
    value::V
    IDs_range::R
end
Adapt.@adapt_structure totalPressure

@define_boundary totalPressure Laplacian{Linear} begin
    J = term.flux[fID]
    (; area, delta) = face 
    flux = J*area/delta
    ap = term.sign*(-flux)
    
    U = U_field[i]
    rho = rho_field[i]
    gh = gh_field[i]

    U_mag = norm(U)
    # println(bc.value - 0.5 * (U_mag^2)*ap)
    ap, ap*(bc.value - rho*gh - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{Linear} begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    0.0, ap*(bc.value - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{Upwind} begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    0.0, ap*(bc.value - 0.5 * (U_mag^2)*rho)
end

@define_boundary Dirichlet Divergence{LUST} begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    0.0, ap*(bc.value - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{BoundedUpwind} begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    flux, ap*(bc.value - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Laplacian{Linear} VectorField begin
    J = term.flux[fID]
    (; area, delta) = face 
    flux = J*area/delta
    ap = term.sign*(-flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    ap, ap*(bc.value[component.value] - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{Linear} VectorField begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    0.0, ap*(bc.value[component.value] - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{Upwind} VectorField begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    0.0, ap*(bc.value[component.value] - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{LUST} VectorField begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    0.0, ap*(bc.value[component.value] - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Divergence{BoundedUpwind} VectorField begin
    flux = -term.flux[fID]
    ap = term.sign*(flux)

    U = U_field[i]
    rho = rho_field[i]

    U_mag = norm(U)
    flux, ap*(bc.value[component.value] - 0.5 * (U_mag^2)*rho)
end

@define_boundary totalPressure Si begin
    0.0, 0.0
end
