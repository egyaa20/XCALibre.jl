export HeatFlux

"""
    HeatFlux <: AbstractNeumann

Fixed heat-flux boundary condition for the enthalpy equation. Injects exactly
`q*A` [W] through each boundary face into the adjacent cell, independent of the
local diffusion coefficient: the Laplacian boundary contribution is
`alpha_eff*A*(dh/dn)` and the required gradient is `q/alpha_eff`, so the
coefficient cancels. This avoids the classic error of prescribing a fixed
gradient computed from reference-state properties, which mis-applies the flux
by the ratio of local to reference conductivity (an order of magnitude across
the pseudo-critical line).

# Inputs
- `ID` Boundary name as a symbol (e.g. :wall)
- `value` Heat flux q [W/m²], positive into the domain

# Example
    HeatFlux(:bottom, 5.0e4)
"""
struct HeatFlux{I,V,R<:UnitRange} <: AbstractNeumann
    ID::I
    value::V
    IDs_range::R
end
Adapt.@adapt_structure HeatFlux

@define_boundary HeatFlux Laplacian{Linear} begin
    (; area) = face
    0.0, -term.sign[1]*bc.value*area
end

@define_boundary HeatFlux Divergence{Linear} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

@define_boundary HeatFlux Divergence{Upwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

@define_boundary HeatFlux Divergence{LUST} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

@define_boundary HeatFlux Divergence{BoundedUpwind} begin
    0.0, 0.0
end

@define_boundary HeatFlux Si begin
    0.0, 0.0
end
