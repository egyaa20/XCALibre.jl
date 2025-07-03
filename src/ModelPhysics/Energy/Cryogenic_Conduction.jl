export CryogenicConduction

struct CryogenicConduction{S1,F1,S2,F2,S3} <: AbstractEnergyModel
    T::S1
    Tf::F1
    k::S2            
    kf::F2            
    cp::S3  
end
Adapt.@adapt_structure CryogenicConduction


Energy{CryogenicConduction}(; material::Symbol, rho::Float64) = begin # maybe assign rho based on the material
    args = (material, rho)
    ARGS = typeof(args)
    Energy{CryogenicConduction,ARGS}(args)
end

(energy::Energy{EnergyModel, ARG})(mesh, medium) where {EnergyModel<:CryogenicConduction,ARG} = begin
    T  = ScalarField(mesh)
    Tf = FaceScalarField(mesh)

    k  = ScalarField(mesh)
    kf = FaceScalarField(mesh)

    cp  = ScalarField(mesh)

    CryogenicConduction(T, Tf, k, kf, cp)

end


function initialise(
    energy::CryogenicConduction, model::Physics{T1,ME,M,Tu,E,D,BI}, T_field, rDf, rhocp_field, rho, material, config
) where {T1,ME,M,Tu,E,D,BI} #T?

    mesh = model.domain

    k  = ScalarField(mesh)
    kf = FaceScalarField(mesh)
    cp  = ScalarField(mesh)


    initialise!(k, 10.0) #for testing
    initialise!(cp, 500.0) #for testing
    # compute k (at cells)
    # compute cp (at cells)


    interpolate_harmonic!(kf, k) # should use harmonic interpolation instead of linear (better for k)

    initialise!(rhocp_field, rho)
    @. rhocp_field.values *= cp.values


    @. rDf.values *= (1.0/kf.values)

    return nothing
end




function energy!(
    energy::CryogenicConduction, model::Physics{T1,ME,M,Tu,E,D,BI}, T_field, rDf, rhocp_field, rho, material, config
) where {T1,ME,M,Tu,E,D,BI}

    # THINGS HAPPEN HERE TOO!!!

    return nothing
end
