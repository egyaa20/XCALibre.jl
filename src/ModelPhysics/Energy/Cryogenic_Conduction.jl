export CryogenicConduction

struct CryogenicConduction{S1,F1,S2,F2,S3,F3} <: AbstractEnergyModel
    T::S1
    Tf::F1
    k::S2            
    kf::F2            
    cp::S3  
    cpf:F3          
end
Adapt.@adapt_structure CryogenicConduction


Energy{CryogenicConduction}() = begin
    args = nothing
    ARGS = typeof(args)
    Energy{CryogenicConduction,ARGS}(args)
end

(energy::Energy{EnergyModel, ARG})(mesh, medium) where {EnergyModel<:CryogenicConduction,ARG} = begin
    T  = ScalarField(mesh)
    Tf = FaceScalarField(mesh)

    k  = ScalarField(mesh)
    kf = FaceScalarField(mesh)

    cp  = ScalarField(mesh)
    cpf = FaceScalarField(mesh)


    CryogenicConduction(T, Tf)
end


function initialise(
    energy::CryogenicConduction, model::Physics{T1,ME,M,Tu,E,D,BI}, _config
) where {T1,ME,M,Tu,E,D,BI} #T?
    
    # THINGS HAPPEN HERE!!!

    return energy
end




function energy!(
    energy::CryogenicConduction, model::Physics{T1,ME,M,Tu,E,D,BI}, config
) where {T1,ME,M,Tu,E,D,BI}

    # THINGS HAPPEN HERE TOO!!!

    return nothing
end
