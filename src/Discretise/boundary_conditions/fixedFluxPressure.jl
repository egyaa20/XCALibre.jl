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
    (; area, delta) = face

    mdotf = extra_eqn.model.terms[2].flux[fID]
    rDf = term.phi[fID]

    U_d = nothing

    for bc in extra_BCs
        if fID in bc.IDs_range
            U_d = bc.value
            break
        end
    end

    phi_target = rhof[fID] * norm(U_d) * area
    g = (mdotf - phi_target) / ((rDf * area) + eps())

    diff = mdotf - phi_target

    # println("phi_target: $phi_target, mdotf: $mdotf, diff: $diff")
    0.0, g
end

 # No need for those:
# @define_boundary fixedFluxPressure Divergence{Linear} ScalarField begin
#     flux = term.flux[fID]
#     (; area, delta) = face 
#     ap = term.sign*(flux) 
#     ap, -bc.value*ap*delta # No need for those
# end

# @define_boundary fixedFluxPressure Divergence{Upwind} ScalarField begin
#     flux = term.flux[fID]
#     (; area, delta) = face 
#     ap = term.sign*(flux) 
#     ap, -bc.value*ap*delta
# end

# @define_boundary fixedFluxPressure Divergence{LUST} ScalarField begin
#     flux = term.flux[fID]
#     (; area, delta) = face 
#     ap = term.sign*(flux) 
#     ap, -bc.value*ap*delta
# end

# @define_boundary fixedFluxPressure Divergence{BoundedUpwind} ScalarField begin
#     flux = term.flux[fID]
#     (; area, delta) = face 
#     ap = term.sign*(flux) 
#     ap, -bc.value*ap*delta
# end

@define_boundary fixedFluxPressure Si ScalarField begin
    0.0, 0.0
end