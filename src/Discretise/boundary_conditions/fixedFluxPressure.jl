export fixedFluxPressure

"""
    fixedFluxPressure <: AbstractNeumann

BC for p_rgh walls

# Example
    fixedFluxPressure(:walls)
"""
struct fixedFluxPressure{I,V,R<:UnitRange} <: AbstractNeumann
    ID::I 
    value::V
    IDs_range::R
end
Adapt.@adapt_structure fixedFluxPressure


@define_boundary fixedFluxPressure Laplacian{Linear} ScalarField begin
    # For now this is hard-coded as zero-gradient. To-do extension to any input gradient
    phi = term.phi 
    values = get_values(phi, component)
    J = term.flux[fID]
    (; area, delta) = face 
    flux = J*area
    0.0, flux*bc.value # draft implementation to test!
end

@define_boundary fixedFluxPressure Divergence{Linear} ScalarField begin
    flux = term.flux[fID]
    (; area, delta) = face 
    ap = term.sign*(flux) 
    ap, -bc.value*ap*delta
end

@define_boundary fixedFluxPressure Divergence{Upwind} ScalarField begin
    flux = term.flux[fID]
    (; area, delta) = face 
    ap = term.sign*(flux) 
    ap, -bc.value*ap*delta
end

@define_boundary fixedFluxPressure Divergence{LUST} ScalarField begin
    flux = term.flux[fID]
    (; area, delta) = face 
    ap = term.sign*(flux) 
    ap, -bc.value*ap*delta
end

@define_boundary fixedFluxPressure Divergence{BoundedUpwind} ScalarField begin
    flux = term.flux[fID]
    (; area, delta) = face 
    ap = term.sign*(flux) 
    ap, -bc.value*ap*delta
end

@define_boundary fixedFluxPressure Si ScalarField begin
    0.0, 0.0
end