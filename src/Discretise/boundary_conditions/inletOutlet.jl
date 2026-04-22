export InletOutlet

"""
    InletOutlet <: AbstractDirichlet

OpenFOAM-style `inletOutlet` boundary condition. On faces where flow is
leaving the domain (outflow, `mdotf >= 0`) it behaves as `Zerogradient`; on
faces where flow is entering (inflow, `mdotf < 0`) it behaves as a
`Dirichlet` with the user-supplied `value`. Used to stabilise pressure-driven
outlets when local backflow would otherwise re-inject dispersed phase or
spurious momentum into the domain.

# Inputs
- `ID` Name of the boundary given as a symbol (e.g. :outlet).
- `value` Scalar or Vector value applied on reversing (inflow) faces.

# Example
    InletOutlet(:outlet, 0.0)            # scalar (e.g. alpha)
    InletOutlet(:outlet, [0.0,0.0,0.0])  # vector (e.g. U)
"""
struct InletOutlet{I,V,R<:UnitRange} <: AbstractDirichlet
    ID::I
    value::V
    IDs_range::R
end
Adapt.@adapt_structure InletOutlet


@define_boundary InletOutlet Laplacian{Linear} begin
    0.0, 0.0
end

@define_boundary InletOutlet Divergence{Linear} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap, 0.0
    else
        0.0, -ap*bc.value
    end
end

@define_boundary InletOutlet Divergence{Upwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap, 0.0
    else
        0.0, -ap*bc.value
    end
end

@define_boundary InletOutlet Divergence{LUST} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap, 0.0
    else
        0.0, -ap*bc.value
    end
end

@define_boundary InletOutlet Divergence{BoundedUpwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap-flux, 0.0
    else
        flux, -ap*bc.value
    end
end


@define_boundary InletOutlet Laplacian{Linear} VectorField begin
    0.0, 0.0
end

@define_boundary InletOutlet Divergence{Linear} VectorField begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap, 0.0
    else
        0.0, -ap*bc.value[component.value]
    end
end

@define_boundary InletOutlet Divergence{Upwind} VectorField begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap, 0.0
    else
        0.0, -ap*bc.value[component.value]
    end
end

@define_boundary InletOutlet Divergence{LUST} VectorField begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap, 0.0
    else
        0.0, -ap*bc.value[component.value]
    end
end

@define_boundary InletOutlet Divergence{BoundedUpwind} VectorField begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    if flux >= zero(flux)
        ap-flux, 0.0
    else
        flux, -ap*bc.value[component.value]
    end
end

@define_boundary InletOutlet Si begin
    0.0, 0.0
end
