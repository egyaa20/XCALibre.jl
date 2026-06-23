export HelmholtzTable, lookup_helmholtz, snapshot_phase_setups
export update_phase_properties_from_table!
export Ttoh_helmholtz, htoT_helmholtz
export interp_h_from_T!, invert_T_from_h!

"""
    HelmholtzTable(fluid; pressure, T_range, n_points=100, T_snapshot=nothing)

Pre-tabulated thermophysical properties from the Helmholtz-energy EOS,
sampled at a fixed pressure across `n_points` evenly-spaced temperatures
in `T_range`.

The table caches per-temperature liquid and vapour values for density
(kg/m³), dynamic viscosity (Pa·s), specific heat at constant pressure
(J/(kg·K)), and thermal conductivity (W/(m·K)).

# Phase-1 integration
At solver initialisation, properties are *snapshotted* at `T_snapshot`
(default: midpoint of `T_range`) and the resulting liquid / vapour
values overwrite the corresponding `Phase` constants. The solver then
runs unmodified, but with NIST-grade properties at one operating
point. Live per-cell interpolation requires solver-kernel changes that
are not done in this phase.

# Inputs
- `fluid`      One of `H2()`, `H2_para()`, `N2()`.
- `pressure`   Pressure in Pa at which to tabulate.
- `T_range`    `(T_min, T_max)` in K.
- `n_points`   Number of sample temperatures (default 100; pass more
               near the saturation line for steeper transitions).
- `T_snapshot` Temperature used to extract phase constants for the
               Phase-1 integration. Defaults to the midpoint of
               `T_range`.

# Example
    fluid = Fluid{Multiphase}(
        model = Mixture(diameter=1e-6),
        phases = (Phase(rho=70.0, mu=1.4e-5),
                  Phase(rho=1.4,  mu=1.0e-6)),
        gravity = Gravity([0, -9.81, 0]),
        fluid_properties = HelmholtzTable(
            H2_para();
            pressure = 1.0e6,
            T_range  = (20.0, 35.0),
            n_points = 200,
        ),
    )

The `Phase(...)` constants serve as a fallback if the table is later
disabled; their values get overwritten by the snapshot at solver init.
"""
struct HelmholtzTable{F<:HelmholtzEnergyFluid, V<:AbstractVector{Float64}}
    fluid::F
    pressure::Float64
    T_range::Tuple{Float64, Float64}
    n_points::Int
    T_snapshot::Float64
    T_grid::V
    rho_l::V
    rho_v::V
    mu_l::V
    mu_v::V
    cp_l::V
    cp_v::V
    k_l::V
    k_v::V
    h_l::V          # sensible enthalpy ∫cp dT, liquid branch [J/kg]
    h_v::V          # sensible enthalpy ∫cp dT, vapour branch [J/kg]
end
Adapt.@adapt_structure HelmholtzTable

# Cumulative sensible enthalpy by the trapezoidal rule, referenced to 0 at
# the table's lower edge:  h[1] = 0,  h[i] = h[i-1] + ½(cp[i-1]+cp[i])·ΔT.
# Because cp > 0 the result is *strictly increasing*, so the h-grid is
# invertible — the cp spike integrates to a smooth, monotone ramp rather
# than a peak, which is exactly what lets the enthalpy energy formulation
# step through the pseudo-critical line without the divide-by-cp runaway.
function _cumulative_enthalpy(T::AbstractVector, cp::AbstractVector)
    n = length(T)
    h = zeros(Float64, n)
    @inbounds for i in 2:n
        h[i] = h[i-1] + 0.5 * (cp[i-1] + cp[i]) * (T[i] - T[i-1])
    end
    return h
end

function HelmholtzTable(fluid::HelmholtzEnergyFluid;
                        pressure::Real,
                        T_range::Tuple{<:Real, <:Real},
                        n_points::Integer = 100,
                        T_snapshot::Union{Real, Nothing} = nothing)

    T_lo, T_hi = float(T_range[1]), float(T_range[2])
    T_lo < T_hi || error("HelmholtzTable: T_range[1] must be < T_range[2], got $T_range")
    n_points >= 2 || error("HelmholtzTable: n_points must be >= 2, got $n_points")
    pressure_f = float(pressure)
    T_snap = isnothing(T_snapshot) ? 0.5 * (T_lo + T_hi) : float(T_snapshot)
    T_lo <= T_snap <= T_hi || error("HelmholtzTable: T_snapshot=$T_snap outside T_range=$T_range")

    eos = HelmholtzEnergy(name = fluid)
    Ts = collect(range(T_lo, T_hi, length = n_points))

    rho_l = zeros(Float64, n_points); rho_v = zeros(Float64, n_points)
    mu_l  = zeros(Float64, n_points); mu_v  = zeros(Float64, n_points)
    cp_l  = zeros(Float64, n_points); cp_v  = zeros(Float64, n_points)
    k_l   = zeros(Float64, n_points); k_v   = zeros(Float64, n_points)

    @info "HelmholtzTable: building $(n_points)-point table for $(typeof(fluid)) at $(pressure_f) Pa, T ∈ $(T_range) K"
    for (i, T) in enumerate(Ts)
        result = eos(T, pressure_f)
        rho_l[i], rho_v[i] = result.rho[1], result.rho[2]
        mu_l[i],  mu_v[i]  = result.mu[1],  result.mu[2]
        cp_l[i],  cp_v[i]  = result.cp[1],  result.cp[2]
        k_l[i],   k_v[i]   = result.k[1],   result.k[2]
    end

    # Sensible enthalpy h(T) = ∫_{T_lo}^{T} cp dT' (exact at constant
    # pressure since cp ≡ (∂h/∂T)_p), consistent with the tabulated cp.
    h_l = _cumulative_enthalpy(Ts, cp_l)
    h_v = _cumulative_enthalpy(Ts, cp_v)

    return HelmholtzTable(fluid, pressure_f, (T_lo, T_hi), n_points, T_snap,
                          Ts, rho_l, rho_v, mu_l, mu_v, cp_l, cp_v, k_l, k_v, h_l, h_v)
end

"""
    lookup_helmholtz(table::HelmholtzTable, T)

Linearly interpolate every tabulated property at temperature `T`.
Returns a `NamedTuple` with two-element tuples `(liquid, vapour)` for
`rho`, `mu`, `cp`, `k`. Out-of-range temperatures are clamped to the
table edges (no extrapolation).
"""
function lookup_helmholtz(table::HelmholtzTable, T::Real)
    T_lo, T_hi = table.T_range
    Tf = float(T)
    T_clamped = clamp(Tf, T_lo, T_hi)

    Δ = (T_hi - T_lo) / (table.n_points - 1)
    pos = (T_clamped - T_lo) / Δ
    i_lo = clamp(Int(floor(pos)) + 1, 1, table.n_points - 1)
    i_hi = i_lo + 1
    w = pos - (i_lo - 1)

    @inline interp(a) = (1.0 - w) * a[i_lo] + w * a[i_hi]

    return (
        rho = (interp(table.rho_l), interp(table.rho_v)),
        mu  = (interp(table.mu_l),  interp(table.mu_v)),
        cp  = (interp(table.cp_l),  interp(table.cp_v)),
        k   = (interp(table.k_l),   interp(table.k_v)),
    )
end

"""
    Ttoh_helmholtz(table::HelmholtzTable, T) -> h

Interpolate sensible enthalpy [J/kg] at temperature `T` from the precomputed
monotone `h_l(T)` grid (liquid branch). Out-of-range `T` clamps to the edges.
"""
function Ttoh_helmholtz(table::HelmholtzTable, T::Real)
    T_lo, T_hi = table.T_range
    Tf = clamp(float(T), T_lo, T_hi)
    Δ = (T_hi - T_lo) / (table.n_points - 1)
    pos = (Tf - T_lo) / Δ
    i_lo = clamp(Int(floor(pos)) + 1, 1, table.n_points - 1)
    i_hi = i_lo + 1
    w = pos - (i_lo - 1)
    return (1.0 - w) * table.h_l[i_lo] + w * table.h_l[i_hi]
end

"""
    htoT_helmholtz(table::HelmholtzTable, h) -> T

Invert the strictly-increasing enthalpy grid: given sensible enthalpy `h`
[J/kg], return the temperature `T` [K] with `h_l(T) = h`. Binary search for
the bracketing interval on `h_l`, then linear interpolation. Out-of-range
`h` clamps to the table edges. This is the monotone, runaway-free `T(h)`
inverse that replaces the divide-by-cp update of the temperature equation.
"""
function htoT_helmholtz(table::HelmholtzTable, h::Real)
    hl = table.h_l
    n  = table.n_points
    hf = float(h)
    hf <= hl[1] && return table.T_grid[1]
    hf >= hl[n] && return table.T_grid[n]
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if hl[mid] <= hf
            lo = mid
        else
            hi = mid
        end
    end
    w = (hf - hl[lo]) / (hl[hi] - hl[lo])
    return (1.0 - w) * table.T_grid[lo] + w * table.T_grid[hi]
end

"""
    snapshot_phase_setups(phase_setups, table::HelmholtzTable)

Phase-1 helper: given a tuple of two `Phase` setups, return a new tuple
where each phase's `rho`, `mu`, `cp`, `k` are replaced by the values
looked up from `table` at `table.T_snapshot`. The first phase is
treated as the liquid branch, the second as the vapour branch — this
matches the two-phase indexing convention used elsewhere in
`build_phase` and the multiphase solver.

This function does not mutate; it constructs new `Phase` objects.
"""
function snapshot_phase_setups(phase_setups::NTuple{2, Phase}, table::HelmholtzTable)
    props = lookup_helmholtz(table, table.T_snapshot)
    return (
        Phase(rho = props.rho[1], mu = props.mu[1], cp = props.cp[1], k = props.k[1]),
        Phase(rho = props.rho[2], mu = props.mu[2], cp = props.cp[2], k = props.k[2]),
    )
end

# Hook into the build_multiphase pipeline (defined in 2_fluid_models.jl).
# When `fluid_properties = HelmholtzTable(...)` is present on
# `Fluid{Multiphase}`, the phase setups get rewritten with snapshotted
# constants before they hit `build_phase`.
_snapshot_from_fluid_properties(phase_setups, table::HelmholtzTable) =
    snapshot_phase_setups(phase_setups, table)

# Strip `fluid_properties` out of the per-property build pipeline — it's
# not a `Gravity`-like field, just a config object that the solver can
# read back later. The lookup arrays are moved to the mesh's backend so
# the per-cell live update can run as a GPU kernel.
build_property(table::HelmholtzTable, mesh) = adapt(_get_backend(mesh), table)

"""
    update_phase_properties_from_table!(phases, table::HelmholtzTable, T_field)

Phase-2 live update. For each cell, looks up `(rho, mu, cp, k)` from
`table` at temperature `T_field.values[i]` and writes the result into
the per-cell `ScalarField`s carried by the two phases. The first phase
receives the liquid branch, the second the vapour branch.

This routine assumes the phases were built with `build_phase_table_mode`
(i.e. their `rho/mu/cp/k` are `ScalarField`s rather than
`ConstantScalar`s). When that is not the case (no table set, or
constant-property mode), the function is a no-op.
"""
function update_phase_properties_from_table!(phases, table::HelmholtzTable, T_field, config)
    # Bail if the phases aren't field-backed (i.e. snapshot mode only).
    phases[1].rho isa ScalarField || return nothing

    (; backend, workgroup) = config.hardware
    Ts = T_field.values
    ndrange = length(Ts)

    kernel! = _interp_phase_properties!(_setup(backend, workgroup, ndrange)...)
    kernel!(
        phases[1].rho.values, phases[2].rho.values,
        phases[1].mu.values,  phases[2].mu.values,
        phases[1].cp.values,  phases[2].cp.values,
        phases[1].k.values,   phases[2].k.values,
        Ts,
        table.rho_l, table.rho_v, table.mu_l, table.mu_v,
        table.cp_l,  table.cp_v,  table.k_l,  table.k_v,
        table.T_range[1], table.T_range[2], table.n_points,
    )
    return nothing
end

@kernel inbounds=true function _interp_phase_properties!(
        rho1, rho2, mu1, mu2, cp1, cp2, k1, k2, Ts,
        rho_l, rho_v, mu_l, mu_v, cp_l, cp_v, k_l, k_v,
        T_lo, T_hi, npts)
    i = @index(Global)

    Δ = (T_hi - T_lo) / (npts - 1)
    T_clamped = clamp(Ts[i], T_lo, T_hi)
    pos  = (T_clamped - T_lo) / Δ
    i_lo = clamp(unsafe_trunc(Int, floor(pos)) + 1, 1, npts - 1)
    i_hi = i_lo + 1
    w    = pos - (i_lo - 1)
    w1   = 1.0 - w

    rho1[i] = w1*rho_l[i_lo] + w*rho_l[i_hi]
    rho2[i] = w1*rho_v[i_lo] + w*rho_v[i_hi]
    mu1[i]  = w1*mu_l[i_lo]  + w*mu_l[i_hi]
    mu2[i]  = w1*mu_v[i_lo]  + w*mu_v[i_hi]
    cp1[i]  = w1*cp_l[i_lo]  + w*cp_l[i_hi]
    cp2[i]  = w1*cp_v[i_lo]  + w*cp_v[i_hi]
    k1[i]   = w1*k_l[i_lo]   + w*k_l[i_hi]
    k2[i]   = w1*k_v[i_lo]   + w*k_v[i_hi]
end

# No-op for any other fluid_properties type
update_phase_properties_from_table!(phases, _other, T_field, config) = nothing

# =====================================================================
# Per-cell enthalpy <-> temperature conversion (for VariableSensibleEnthalpy)
# =====================================================================

"""
    interp_h_from_T!(h_field, T_field, table::HelmholtzTable, config)

Per-cell sensible enthalpy from temperature using the monotone `h_l(T)`
grid (liquid branch). Used to initialise `h` consistently from a `T` field.
"""
function interp_h_from_T!(h_field, T_field, table::HelmholtzTable, config)
    (; backend, workgroup) = config.hardware
    ndrange = length(T_field.values)
    kernel! = _interp_h_from_T!(_setup(backend, workgroup, ndrange)...)
    kernel!(h_field.values, T_field.values, table.h_l,
            table.T_range[1], table.T_range[2], table.n_points)
    return nothing
end

@kernel inbounds=true function _interp_h_from_T!(hs, Ts, h_l, T_lo, T_hi, npts)
    i = @index(Global)
    Δ    = (T_hi - T_lo) / (npts - 1)
    Tc   = clamp(Ts[i], T_lo, T_hi)
    pos  = (Tc - T_lo) / Δ
    i_lo = clamp(unsafe_trunc(Int, floor(pos)) + 1, 1, npts - 1)
    i_hi = i_lo + 1
    w    = pos - (i_lo - 1)
    hs[i] = (1.0 - w) * h_l[i_lo] + w * h_l[i_hi]
end

"""
    invert_T_from_h!(T_field, h_field, table::HelmholtzTable, config)

Per-cell temperature from sensible enthalpy by inverting the strictly
increasing `h_l` grid (binary search + linear interpolation). This is the
monotone, runaway-free `T(h)` update that replaces dividing by the spiky cp
in the temperature equation. Out-of-range `h` clamps to the table edges.
"""
function invert_T_from_h!(T_field, h_field, table::HelmholtzTable, config)
    (; backend, workgroup) = config.hardware
    ndrange = length(h_field.values)
    kernel! = _invert_T_from_h!(_setup(backend, workgroup, ndrange)...)
    kernel!(T_field.values, h_field.values, table.h_l, table.T_grid, table.n_points)
    return nothing
end

@kernel inbounds=true function _invert_T_from_h!(Ts, hs, h_l, T_grid, npts)
    i = @index(Global)
    h = hs[i]
    if h <= h_l[1]
        Ts[i] = T_grid[1]
    elseif h >= h_l[npts]
        Ts[i] = T_grid[npts]
    else
        lo = 1; hi = npts
        while hi - lo > 1
            mid = (lo + hi) >>> 1
            if h_l[mid] <= h
                lo = mid
            else
                hi = mid
            end
        end
        w = (h - h_l[lo]) / (h_l[hi] - h_l[lo])
        Ts[i] = (1.0 - w) * T_grid[lo] + w * T_grid[hi]
    end
end
