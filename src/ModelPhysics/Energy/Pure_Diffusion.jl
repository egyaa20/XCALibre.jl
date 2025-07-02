export PureDiffusion

struct PureDiffusion{S1,F1} <: AbstractEnergyModel
    T::S1
    Tf::F1
end
Adapt.@adapt_structure PureDiffusion

# Energy API constructor: allow `Energy{LaplaceEnergy}()`

# Energy{LaplaceEnergy}() = Energy{LaplaceEnergy,Nothing}(nothing)


Energy{PureDiffusion}() = begin
    args = nothing
    ARGS = typeof(args)
    Energy{PureDiffusion,ARGS}(args)
end

(energy::Energy{EnergyModel, ARG})(mesh, medium) where {EnergyModel<:PureDiffusion,ARG} = begin
    T  = ScalarField(mesh)
    Tf = FaceScalarField(mesh)
    PureDiffusion(T, Tf)
end


function initialise(
    energy::PureDiffusion, model::Physics{T1,ME,M,Tu,E,D,BI}, _config
) where {T1,ME,M,Tu,E,D,BI} #T?
    # nothing special to set up for pure Laplace
    return energy
end




function energy!(
    energy::PureDiffusion, model::Physics{T1,ME,M,Tu,E,D,BI}, config
) where {T1,ME,M,Tu,E,D,BI}
    return nothing
end
