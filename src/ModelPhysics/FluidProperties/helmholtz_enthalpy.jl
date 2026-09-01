export HelmholtzEnthalpy

"""
    Energy{HelmholtzEnthalpy}(; Pr_t=0.85)

Enthalpy energy model driven by a `HelmholtzTable` equation of state:

    d(rho*h)/dt + div(rho*U*h) = laplacian(lambda/cp + mu_t/Pr_t, h)

with the T <-> h inversion provided by the table. Requires a phase with
`rho=HelmholtzTable(...)`. Initialise `model.energy.T` before `run!`.
"""
struct HelmholtzEnthalpy{S,FS,C} <: AbstractEnergyModel
    h::S
    T::S
    hf::FS
    Tf::FS
    coeffs::C
end
Adapt.@adapt_structure HelmholtzEnthalpy

Energy{HelmholtzEnthalpy}(; Pr_t=0.85) = begin
    coeffs = (Pr_t=Pr_t,)
    ARG = typeof(coeffs)
    Energy{HelmholtzEnthalpy,ARG}(coeffs)
end

(energy::Energy{HelmholtzEnthalpy,ARG})(mesh, fluid) where {ARG} = begin
    h  = ScalarField(mesh)
    T  = ScalarField(mesh)
    hf = FaceScalarField(mesh)
    Tf = FaceScalarField(mesh)
    HelmholtzEnthalpy(h, T, hf, Tf, energy.args)
end
