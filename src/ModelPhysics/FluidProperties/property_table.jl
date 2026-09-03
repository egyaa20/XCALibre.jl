export HelmholtzTable, HelmholtzMu, HeatFluxGradient
export update_phase_properties!, enthalpy_from_temperature!, temperature_from_enthalpy!
export table_enthalpy, table_temperature

struct HelmholtzMu <: AbstractViscosityModel end   # properties supplied by the rho-slot HelmholtzTable

_T_crit(::N2) = 126.192;      _p_crit(::N2) = 33.958e5
_T_crit(::H2) = 33.145;       _p_crit(::H2) = 12.964e5
_T_crit(::H2_para) = 32.938;  _p_crit(::H2_para) = 12.858e5

"""
    HelmholtzTable(; fluid, p_ref, T_min, T_max, n_points=500) <: AbstractEosModel

Uniform-in-T lookup table of rho, mu, k, cp, h built at construction by sampling the
Helmholtz EOS at fixed reference pressure `p_ref` [Pa] over `T_min..T_max` [K] with
`n_points` samples. Includes the inverse T(h) table for enthalpy-based energy models.
Use as `Phase(rho=HelmholtzTable(...), mu=HelmholtzMu())`.
"""
struct HelmholtzTable{F<:AbstractFloat, V<:AbstractVector{F}, FL<:HelmholtzEnergyFluid} <: AbstractEosModel
    fluid::FL
    p_ref::F
    T_min::F
    T_max::F
    dT::F
    rho::V
    mu::V
    k::V
    cp::V
    h::V
    h_min::F
    h_max::F
    dh::F
    T_of_h::V
end
Adapt.@adapt_structure HelmholtzTable

function _branch_index(T, T_sat, Tc, p_ref, pc)
    T >= Tc && return 1        # supercritical: both entries identical
    p_ref >= pc && return 1    # compressed liquid above critical pressure
    return T >= T_sat ? 2 : 1  # vapour : liquid
end

function HelmholtzTable(; fluid, p_ref, T_min, T_max, n_points::Int=500, h_ref_T=nothing)
    @assert T_max > T_min "T_max must exceed T_min"
    @assert n_points >= 2 "n_points must be at least 2"
    F = Float64

    eos = HelmholtzEnergy(name=fluid)
    Tc = F(_T_crit(fluid)); pc = F(_p_crit(fluid))

    T_grid = collect(range(F(T_min), F(T_max), length=n_points))
    dT = T_grid[2] - T_grid[1]

    rho = zeros(F, n_points); mu = zeros(F, n_points); k = zeros(F, n_points)
    cp  = zeros(F, n_points); h  = zeros(F, n_points)

    for (i, T) in enumerate(T_grid)
        res = eos(T, F(p_ref))
        idx = _branch_index(T, res.T_sat, Tc, F(p_ref), pc)
        rho[i] = res.rho[idx]
        mu[i]  = res.mu[idx]
        k[i]   = res.k[idx]
        cp[i]  = res.cp[idx] * F(1e3)   # J/g -> J/kg
        h[i]   = res.h[idx]  * F(1e3)
    end

    all(isfinite, rho) && all(>(zero(F)), rho) || error("HelmholtzTable: non-physical density in T range")
    all(>(zero(F)), cp) || error("HelmholtzTable: non-positive cp in T range")
    issorted(h) || error("HelmholtzTable: h(T) not monotonically increasing — cannot build T(h) inverse")

    # Reference-shift so h = 0 at the bulk/operating temperature. The conservative
    # form d(rho h)/dt + div(rho U h) carries a spurious source h*[drho/dt +
    # div(rho U)] proportional to |h|, since the bracket vanishes only to discrete
    # tolerance. Anchoring at the operating point (not the table floor) is what
    # makes it small where most of the field lives.
    if h_ref_T !== nothing
        h .-= _table_lerp(h, F(h_ref_T), F(T_min), dT, n_points)
    end

    h_min = h[1]; h_max = h[end]
    dh = (h_max - h_min) / (n_points - 1)
    T_of_h = zeros(F, n_points)
    j = 1
    for i in 1:n_points
        hi = h_min + (i - 1) * dh
        while j < n_points - 1 && h[j+1] < hi
            j += 1
        end
        w = (hi - h[j]) / (h[j+1] - h[j] + eps(F))
        T_of_h[i] = T_grid[j] * (1 - w) + T_grid[j+1] * w
    end

    HelmholtzTable(fluid, F(p_ref), F(T_min), F(T_max), dT,
                   rho, mu, k, cp, h, h_min, h_max, dh, T_of_h)
end

@inline function _table_lerp(vals, x, x_min, dx, n)
    ξ = clamp((x - x_min) / dx, zero(x), oftype(x, n - 1))
    i = min(unsafe_trunc(Int, ξ), n - 2)
    w = ξ - i
    @inbounds vals[i+1] * (1 - w) + vals[i+2] * w
end

table_enthalpy(table::HelmholtzTable, T) =
    _table_lerp(table.h, float(T), table.T_min, table.dT, length(table.h))
table_temperature(table::HelmholtzTable, h) =
    _table_lerp(table.T_of_h, float(h), table.h_min, table.dh, length(table.T_of_h))

update_phase_properties!(phase, T_field, config) =
    update_phase_properties!(phase.rho_model, phase, T_field, config)
update_phase_properties!(model, phase, T_field, config) = nothing
update_phase_properties!(model, phase, T_field::Nothing, config) = nothing

function update_phase_properties!(table::HelmholtzTable, phase, T_field, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(T_field.values)
    kernel! = _update_phase_properties!(_setup(backend, workgroup, ndrange)...)
    kernel!(phase.rho.values, phase.mu.values, phase.k.values, phase.cp.values,
            T_field.values, table)
end

@kernel inbounds=true function _update_phase_properties!(rho, mu, k, cp, T, table)
    i = @index(Global)
    n = length(table.rho)
    Ti = T[i]
    rho[i] = _table_lerp(table.rho, Ti, table.T_min, table.dT, n)
    mu[i]  = _table_lerp(table.mu,  Ti, table.T_min, table.dT, n)
    k[i]   = _table_lerp(table.k,   Ti, table.T_min, table.dT, n)
    cp[i]  = _table_lerp(table.cp,  Ti, table.T_min, table.dT, n)
end

function enthalpy_from_temperature!(h_field, T_field, table::HelmholtzTable, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(h_field.values)
    kernel! = _enthalpy_from_temperature!(_setup(backend, workgroup, ndrange)...)
    kernel!(h_field.values, T_field.values, table)
end

@kernel inbounds=true function _enthalpy_from_temperature!(h, T, table)
    i = @index(Global)
    h[i] = _table_lerp(table.h, T[i], table.T_min, table.dT, length(table.h))
end

function temperature_from_enthalpy!(T_field, h_field, table::HelmholtzTable, config)
    (; hardware) = config
    (; backend, workgroup) = hardware
    ndrange = length(T_field.values)
    kernel! = _temperature_from_enthalpy!(_setup(backend, workgroup, ndrange)...)
    kernel!(T_field.values, h_field.values, table)
end

@kernel inbounds=true function _temperature_from_enthalpy!(T, h, table)
    i = @index(Global)
    T[i] = _table_lerp(table.T_of_h, h[i], table.h_min, table.dh, length(table.T_of_h))
end

"""
    HeatFluxGradient(name::Symbol, q, table::HelmholtzTable)

Fixed heat-flux wall for the enthalpy equation applied as a gradient condition:
`dh/dn = q * cp(T_w) / lambda(T_w)`, with the wall temperature re-evaluated from
the local boundary-cell enthalpy through `table` at every assembly. Never uses
frozen reference-state properties, so the applied flux tracks the collapse of
`lambda` across the pseudo-critical line. `q` in W/m², positive into the domain.
See also [`HeatFlux`](@ref), which injects `q*A` exactly by construction.
"""
struct HeatFluxGradient{I,V,R<:UnitRange} <: AbstractNeumann
    ID::I
    value::V
    IDs_range::R
end
Adapt.@adapt_structure HeatFluxGradient

struct HeatFluxGradientValue{F<:AbstractFloat,T}
    q::F
    table::T
end
Adapt.@adapt_structure HeatFluxGradientValue

HeatFluxGradient(name::Symbol, q::Number, table::HelmholtzTable) =
    HeatFluxGradient(name, HeatFluxGradientValue(Float64(q), table), 0:0)

Discretise.adapt_value(value::HeatFluxGradientValue, mesh) =
    HeatFluxGradientValue(_get_float(mesh)(value.q), value.table)

Discretise.@define_boundary HeatFluxGradient Laplacian{Linear} begin
    (; q, table) = bc.value
    J = term.flux[fID]
    (; area) = face
    h_c = get_values(term.phi, component)[cellID]
    n = length(table.T_of_h)
    Tw = _table_lerp(table.T_of_h, h_c, table.h_min, table.dh, n)
    lam = _table_lerp(table.k, Tw, table.T_min, table.dT, n)
    cpw = _table_lerp(table.cp, Tw, table.T_min, table.dT, n)
    g = q * cpw / (lam + eps(typeof(lam)))
    0.0, -term.sign[1]*J*g*area
end

Discretise.@define_boundary HeatFluxGradient Divergence{Linear} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

Discretise.@define_boundary HeatFluxGradient Divergence{Upwind} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

Discretise.@define_boundary HeatFluxGradient Divergence{LUST} begin
    flux = term.flux[fID]
    ap = term.sign*(flux)
    ap, 0.0
end

Discretise.@define_boundary HeatFluxGradient Divergence{BoundedUpwind} begin
    0.0, 0.0
end

Discretise.@define_boundary HeatFluxGradient Si begin
    0.0, 0.0
end

@inline function Discretise.boundary_interpolation!(
    BC::HeatFluxGradient, phif::FaceScalarField, phi, boundary_cellsID, time, fID)
    @inbounds begin
        cID = boundary_cellsID[fID]
        phif[fID] = phi[cID]
    end
    nothing
end

@inline function Discretise.boundary_interpolation!(
    BC::HeatFluxGradient, psif::FaceVectorField, psi, boundary_cellsID, time, fID)
    @inbounds begin
        cID = boundary_cellsID[fID]
        psif[fID] = psi[cID]
    end
    nothing
end
