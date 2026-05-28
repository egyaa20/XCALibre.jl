export NeumannFunction

# abstract type XCALibreUserFunctor end

"""
    NeumannFunction(ID, value) <: AbstractNeumann

Fixed-flux-density boundary condition driven by a user-supplied function,
mirroring `DirichletFunction` for fixed-value BCs. The function returns
the diffusive flux density at the face (e.g. heat flux q [W/m²] when the
equation being solved is a temperature transport equation), and the
discrete contribution to the Laplacian boundary integral is

    ∫_face flux · n dA  ≈  value(face.centre, time, i) · area

The diffusion coefficient at the face does not enter the BC — the user
asserts the physical flux density directly.

# Inputs
- `ID`    Name of the boundary patch (e.g. `:left_heated`).
- `value` Function `(coords, time, faceID) -> q` returning the flux density.

# Example
    # Heat flux of 10 W/m² into the domain through `:left_heated`
    q_left(coords, time, faceID) = 10.0
    NeumannFunction(:left_heated, q_left)
"""
struct NeumannFunction{I,V,R<:UnitRange} <: AbstractNeumann
    ID::I
    value::V
    IDs_range::R
end
Adapt.@adapt_structure NeumannFunction

@define_boundary NeumannFunction Laplacian{Linear} ScalarField begin
    # Boundary contribution to b: + q(face,time) · area
    # a_p = 0 (no coupling to the cell value); no dependence on the face
    # diffusion coefficient.
    (; centre, area) = face
    q = bc.value(centre, time, i)
    0.0, q * area
end

@define_boundary NeumannFunction Divergence{Linear} begin
    flux = term.flux[fID]
    ap = term.sign*(flux) 
    ap, 0.0 # original

    # phi = term.phi 
    # values = get_values(phi, component)
    # 0.0, -ap*values[cellID] # try this
end

@define_boundary NeumannFunction Divergence{Upwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux) 
    ap, 0.0 # original

    # phi = term.phi 
    # values = get_values(phi, component)
    # 0.0, -ap*values[cellID] # try this
end

@define_boundary NeumannFunction Divergence{LUST} begin
    flux = term.flux[fID]
    ap = term.sign*(flux) 
    ap, 0.0 # original

    # phi = term.phi 
    # values = get_values(phi, component)
    # 0.0, -ap*values[cellID] # try this
end

@define_boundary NeumannFunction Divergence{BoundedUpwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap-flux, 0.0
end

@define_boundary NeumannFunction Si begin
    0.0, 0.0
end