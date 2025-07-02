export CryogenicConduction

struct CryogenicConduction{S1,F1} <: AbstractEnergyModel
    T::S1
    Tf::F1
end
Adapt.@adapt_structure CryogenicConduction

# Energy API constructor: allow `Energy{LaplaceEnergy}()`

# Energy{LaplaceEnergy}() = Energy{LaplaceEnergy,Nothing}(nothing)


Energy{CryogenicConduction}() = begin
    args = nothing
    ARGS = typeof(args)
    Energy{CryogenicConduction,ARGS}(args)
end

(energy::Energy{EnergyModel, ARG})(mesh, medium) where {EnergyModel<:CryogenicConduction,ARG} = begin
    T  = ScalarField(mesh)
    Tf = FaceScalarField(mesh)
    CryogenicConduction(T, Tf)
end


function initialise(
    energy::CryogenicConduction, model::Physics{T1,ME,M,Tu,E,D,BI}, _config
) where {T1,ME,M,Tu,E,D,BI} #T?
    # nothing special to set up for pure Laplace
    return energy
end




function energy!(
    energy::CryogenicConduction, model::Physics{T1,ME,M,Tu,E,D,BI}, config
) where {T1,ME,M,Tu,E,D,BI}
    return nothing
end
