# # abstract type AbstractFluid end

# #maybe sub types for liquid and gas, and more complex is the fluid

# struct Air <: AbstractFluid end
# struct Water <: AbstractFluid end

# struct H2 <: AbstractFluid end
# struct N2 <: AbstractFluid end

# # Air has rho=1.225, nu=1.48e-5, R=287
# # Water has rho=1000.0, nu=1.0e-6, R=?, not a gas...



# function EoS_default(p::Float64, T::Float64, R::Float64) # ONLY FOR GAS;    liquid ~ const

#     return p/(R*T) #Returns density
# end


# function Viscosity_Const_Mix!(nueff::FaceScalarField, alphaf::FaceScalarField, nu_0::Float64, nu_1::Float64)

#     @. nueff.values = (nu_1 * alphaf.values) + (nu_0 * (1 - alphaf.values))
#     #maybe clamp alpha in here?

# end

# function Density_Const_Mix!(rho::ScalarField, alpha::ScalarField, rho_0::Float64, rho_1::Float64)

#     @. rho.values = (rho_1 * alpha.values) + (rho_0 * (1 - alpha.values))
#     #maybe clamp alpha in here?

# end
