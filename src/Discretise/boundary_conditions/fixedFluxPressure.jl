export FixedFluxPressure


"""
    FixedFluxPressure <: AbstractNeumann

Neumann-type pressure boundary condition that mirrors OpenFOAM's
`fixedFluxPressure`. The face-normal gradient is stored per face in the
`value` array and is populated by the solver before the pressure solve so
that the pressure correction exactly cancels the buoyancy (or any other
explicit) face flux at the patch. This is the standard treatment for p_rgh
at a pressure outlet under variable density (multiphase + gravity), where a
fixed value BC would be inconsistent with the lateral ∇ρ·gh contribution.

# Input
- `ID` Name of the boundary given as a symbol (e.g. :outlet). Internally
  replaced with the boundary index.

# Example
    FixedFluxPressure(:outlet)

The per-face gradient is allocated as a zero-filled array during `assign`
(sized to the patch) and then overwritten by the solver.
"""
struct FixedFluxPressure{I,V,R<:UnitRange} <: AbstractNeumann
    ID::I
    value::V
    IDs_range::R
end
Adapt.@adapt_structure FixedFluxPressure

FixedFluxPressure(name::Symbol) = FixedFluxPressure(name, Float64[], 0:0)

# Called from assign_patches instead of adapt_value so the gradient array is
# sized to the patch length (not treated as a 3-component velocity vector).
adapt_value_for_bc(BC::FixedFluxPressure, region, IDs_range) = begin
    F = _get_float(region)
    zeros(F, length(IDs_range))
end


@define_boundary FixedFluxPressure Laplacian{Linear} begin
    J = term.flux[fID]
    (; area) = face
    flux = J*area
    local_i = fID - bc.IDs_range.start + 1
    grad = bc.value[local_i]
    0.0, flux*grad
end

@define_boundary FixedFluxPressure Divergence{Linear} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

@define_boundary FixedFluxPressure Divergence{Upwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

@define_boundary FixedFluxPressure Divergence{LUST} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

@define_boundary FixedFluxPressure Divergence{BoundedUpwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap-flux, 0.0
end

@define_boundary FixedFluxPressure Si begin
    0.0, 0.0
end
