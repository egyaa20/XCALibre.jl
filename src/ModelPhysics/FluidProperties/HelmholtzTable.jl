export HelmholtzTable, lookup_helmholtz, snapshot_phase_setups
export update_phase_properties_from_table!, update_psi_drhodT_from_table!
export Ttoh_helmholtz, htoT_helmholtz
export interp_h_from_T!, invert_T_from_h!

"""
    HelmholtzTable(fluid; p_operating, T_range, p_range, n_T, n_p, ...)

Pre-tabulated thermophysical properties from the Helmholtz-energy EOS over a 2D
`(T, p)` grid. Replaces the earlier constant-pressure (T-only) table so the
solver can resolve the pressure dependence of density needed by the compressible
pressure equation.

Stored on the grid (each an `[n_T, n_p]` matrix):
- `rho`    : density [kg/m³]
- `mu`     : dynamic viscosity [Pa·s]
- `cp`     : isobaric specific heat [J/(kg·K)]
- `k`      : thermal conductivity [W/(m·K)]
- `h`      : sensible enthalpy ∫cp dT integrated **along T at each pressure** [J/kg]
- `psi`    : ψ = (∂ρ/∂p)_T  [kg/(m³·Pa)]   — implicit pressure-term coefficient
- `drhodT` : (∂ρ/∂T)_p      [kg/(m³·K)]    — explicit thermal-expansion source

The **temperature** grid is non-uniform: refined around the fluid critical
temperature `T_crit`, where cp peaks and density changes fastest. Spacing is the
base average everywhere except `|T − T_crit| ≤ mid_half` (×`mid_mult` denser) and
`|T − T_crit| ≤ fine_half` (×`fine_mult` denser). The **pressure** grid is uniform.

The table is indexed by **absolute** pressure. At runtime the cell pressure is
`p_operating + p_rgh` (the dynamic `p_rgh` is gauge), so `p_range` must bracket the
operating pressure plus the expected hydrodynamic swing.

# Single-phase note
For the intended supercritical operating range (every `p ≥ p_c`) there is no
liquid–vapour split, so a single property branch is stored. The live update writes
the same value into both `Phase` slots, preserving the two-phase solver interface.

# Inputs
- `fluid`        One of `H2()`, `H2_para()`, `N2()`.
- `p_operating`  Absolute operating pressure [Pa] (the table is indexed by `p_op + p_rgh`).
- `T_range`      `(T_min, T_max)` in K.
- `p_range`      `(p_min, p_max)` in Pa (default: ±3 % around `p_operating`).
- `n_T`          Base (average) number of T samples; refined zones add more.
- `n_p`          Number of pressure samples (uniform).
- `T_snapshot`   T at which `Phase` constants are snapshotted (default: T midpoint).
- `T_crit`       Critical temperature for the refinement centre (default: EOS value).
- `fine_half`, `mid_half`, `fine_mult`, `mid_mult`  Refinement window half-widths [K]
                 and density multipliers.
"""
struct HelmholtzTable{Fl<:HelmholtzEnergyFluid, V<:AbstractVector{Float64}, M<:AbstractMatrix{Float64}}
    fluid::Fl
    p_operating::Float64
    T_range::Tuple{Float64, Float64}
    p_range::Tuple{Float64, Float64}
    n_T::Int
    n_p::Int
    T_snapshot::Float64
    T_grid::V                 # non-uniform, refined near T_crit (length n_T)
    p_grid::V                 # uniform (length n_p)
    rho::M
    mu::M
    cp::M
    k::M
    h::M                      # ∫cp dT along T, per pressure column
    psi::M                    # (∂ρ/∂p)_T
    drhodT::M                 # (∂ρ/∂T)_p
end
Adapt.@adapt_structure HelmholtzTable

# Critical temperature used to centre the T-refinement (from the EOS constants).
_fluid_Tcrit(fluid::HelmholtzEnergyFluid) = helmholtz_constants(fluid).T_c

"""
    _refined_T_grid(T_lo, T_hi, T_crit; n_avg, fine_half, mid_half, fine_mult, mid_mult)

Build a strictly-increasing, non-uniform T grid. The base spacing `Δavg =
(T_hi-T_lo)/(n_avg-1)` applies away from `T_crit`; it is divided by `mid_mult`
inside `|T-T_crit| ≤ mid_half` and by `fine_mult` inside `|T-T_crit| ≤ fine_half`.
Built segment-by-segment then de-duplicated at the shared segment endpoints.
"""
function _refined_T_grid(T_lo::Float64, T_hi::Float64, T_crit::Float64;
                         n_avg::Int, fine_half::Float64, mid_half::Float64,
                         fine_mult::Real, mid_mult::Real)
    Δavg = (T_hi - T_lo) / (n_avg - 1)
    bnds = (T_lo,
            clamp(T_crit - mid_half,  T_lo, T_hi),
            clamp(T_crit - fine_half, T_lo, T_hi),
            clamp(T_crit + fine_half, T_lo, T_hi),
            clamp(T_crit + mid_half,  T_lo, T_hi),
            T_hi)
    spac = (Δavg, Δavg/mid_mult, Δavg/fine_mult, Δavg/mid_mult, Δavg)
    pts = Float64[]
    for s in 1:5
        a, b = bnds[s], bnds[s+1]
        b > a || continue
        n = max(2, round(Int, (b - a) / spac[s]) + 1)
        append!(pts, collect(range(a, b, length = n)))
    end
    sort!(pts)
    unique!(pts)               # drop the exact duplicates at shared segment ends
    return pts
end

function HelmholtzTable(fluid::HelmholtzEnergyFluid;
                        p_operating::Real,
                        T_range::Tuple{<:Real, <:Real},
                        p_range::Union{Tuple{<:Real, <:Real}, Nothing} = nothing,
                        n_T::Integer = 3000,
                        n_p::Integer = 10,
                        T_snapshot::Union{Real, Nothing} = nothing,
                        T_crit::Union{Real, Nothing} = nothing,
                        fine_half::Real = 1.0, mid_half::Real = 10.0,
                        fine_mult::Real = 5, mid_mult::Real = 3)

    T_lo, T_hi = float(T_range[1]), float(T_range[2])
    T_lo < T_hi || error("HelmholtzTable: T_range[1] must be < T_range[2], got $T_range")
    n_T >= 2 || error("HelmholtzTable: n_T must be >= 2, got $n_T")
    n_p >= 2 || error("HelmholtzTable: n_p must be >= 2, got $n_p")

    p_op = float(p_operating)
    p_lo, p_hi = isnothing(p_range) ? (0.97 * p_op, 1.03 * p_op) :
                                       (float(p_range[1]), float(p_range[2]))
    p_lo < p_hi || error("HelmholtzTable: p_range[1] must be < p_range[2], got $((p_lo, p_hi))")
    p_lo <= p_op <= p_hi || @warn "HelmholtzTable: p_operating=$p_op outside p_range=$((p_lo, p_hi))"

    T_snap = isnothing(T_snapshot) ? 0.5 * (T_lo + T_hi) : float(T_snapshot)
    T_lo <= T_snap <= T_hi || error("HelmholtzTable: T_snapshot=$T_snap outside T_range=$T_range")
    Tc = isnothing(T_crit) ? _fluid_Tcrit(fluid) : float(T_crit)

    T_grid = _refined_T_grid(T_lo, T_hi, Tc; n_avg = Int(n_T),
                             fine_half = float(fine_half), mid_half = float(mid_half),
                             fine_mult = fine_mult, mid_mult = mid_mult)
    nT     = length(T_grid)
    p_grid = collect(range(p_lo, p_hi, length = Int(n_p)))

    constants = helmholtz_constants(fluid)      # Float64 coefficients
    Mmol      = constants.M
    eos       = HelmholtzEnergy(name = fluid)

    rho    = zeros(Float64, nT, n_p); mu     = zeros(Float64, nT, n_p)
    cp     = zeros(Float64, nT, n_p); k      = zeros(Float64, nT, n_p)
    psi    = zeros(Float64, nT, n_p); drhodT = zeros(Float64, nT, n_p)

    @info "HelmholtzTable: building $(nT)×$(n_p) (T,p) table for $(typeof(fluid)), T ∈ $((T_lo, T_hi)) K (refined @ $(round(Tc, digits=2)) K), p ∈ $((p_lo, p_hi)) Pa"
    Threads.@threads for it in 1:nT
        T = T_grid[it]
        for jp in 1:n_p
            p   = p_grid[jp]
            res = eos(T, p)
            rho_mass    = res.rho[1]
            rho[it, jp] = rho_mass
            mu[it, jp]  = res.mu[1]
            cp[it, jp]  = res.cp[1]
            k[it, jp]   = res.k[1]
            rho_mol         = rho_mass / Mmol
            psi[it, jp]     = compute_psi(T, rho_mol, constants, fluid)
            drhodT[it, jp]  = compute_drho_dT_p(T, rho_mol, constants, fluid)
        end
    end

    # Sensible enthalpy h(T;p) = ∫_{T_lo}^{T} cp(T',p) dT' — integrated along T
    # independently for each pressure column (trapezoidal, strictly increasing in T).
    h = zeros(Float64, nT, n_p)
    for jp in 1:n_p
        @inbounds for it in 2:nT
            h[it, jp] = h[it-1, jp] +
                0.5 * (cp[it-1, jp] + cp[it, jp]) * (T_grid[it] - T_grid[it-1])
        end
    end

    return HelmholtzTable(fluid, p_op, (T_lo, T_hi), (p_lo, p_hi), nT, Int(n_p), T_snap,
                          T_grid, p_grid, rho, mu, cp, k, h, psi, drhodT)
end

# ---------------------------------------------------------------------------
# Index helpers (CPU). T grid is non-uniform → binary search; p grid uniform.
# ---------------------------------------------------------------------------
@inline function _locate_T(T_grid, nT, T)
    Tc = isfinite(T) ? clamp(T, T_grid[1], T_grid[nT]) : T_grid[1]
    lo, hi = 1, nT
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        @inbounds if T_grid[mid] <= Tc
            lo = mid
        else
            hi = mid
        end
    end
    @inbounds w = (Tc - T_grid[lo]) / (T_grid[hi] - T_grid[lo])
    return lo, hi, w
end

@inline function _locate_p(p_lo, p_hi, np, p)
    pc  = isfinite(p) ? clamp(p, p_lo, p_hi) : p_lo
    Δp  = (p_hi - p_lo) / (np - 1)
    pos = (pc - p_lo) / Δp
    jlo = clamp(unsafe_trunc(Int, floor(pos)) + 1, 1, np - 1)
    return jlo, jlo + 1, pos - (jlo - 1)
end

@inline function _bilinear(A, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    @inbounds begin
        lo = (1.0 - wp) * A[iT_lo, jp_lo] + wp * A[iT_lo, jp_hi]
        hi = (1.0 - wp) * A[iT_hi, jp_lo] + wp * A[iT_hi, jp_hi]
    end
    return (1.0 - wT) * lo + wT * hi
end

"""
    lookup_helmholtz(table, T, p=table.p_operating)

Bilinearly interpolate every tabulated property at `(T, p)` (absolute pressure).
Out-of-range `T`/`p` clamp to the table edges. Returns a `NamedTuple` of scalars.
"""
function lookup_helmholtz(table::HelmholtzTable, T::Real, p::Real = table.p_operating)
    iT_lo, iT_hi, wT = _locate_T(table.T_grid, table.n_T, float(T))
    jp_lo, jp_hi, wp = _locate_p(table.p_range[1], table.p_range[2], table.n_p, float(p))
    bil(A) = _bilinear(A, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    return (rho = bil(table.rho), mu = bil(table.mu), cp = bil(table.cp),
            k = bil(table.k), h = bil(table.h), psi = bil(table.psi),
            drhodT = bil(table.drhodT))
end

"""
    Ttoh_helmholtz(table, T, p=table.p_operating) -> h

Sensible enthalpy [J/kg] at `(T, p)` by bilinear interpolation. Used to set the
inlet enthalpy for the h-formulation.
"""
function Ttoh_helmholtz(table::HelmholtzTable, T::Real, p::Real = table.p_operating)
    iT_lo, iT_hi, wT = _locate_T(table.T_grid, table.n_T, float(T))
    jp_lo, jp_hi, wp = _locate_p(table.p_range[1], table.p_range[2], table.n_p, float(p))
    return _bilinear(table.h, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
end

"""
    htoT_helmholtz(table, h, p=table.p_operating) -> T

Invert `h(T;p)` at fixed pressure `p`: interpolate the h-column to the pressure,
then binary-search the strictly-increasing column for `h`. Out-of-range `h` clamps.
"""
function htoT_helmholtz(table::HelmholtzTable, h::Real, p::Real = table.p_operating)
    nT = table.n_T
    jp_lo, jp_hi, wp = _locate_p(table.p_range[1], table.p_range[2], table.n_p, float(p))
    hcol(it) = @inbounds (1.0 - wp) * table.h[it, jp_lo] + wp * table.h[it, jp_hi]
    hf = float(h)
    hf <= hcol(1)  && return table.T_grid[1]
    hf >= hcol(nT) && return table.T_grid[nT]
    lo, hi = 1, nT
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if hcol(mid) <= hf
            lo = mid
        else
            hi = mid
        end
    end
    w = (hf - hcol(lo)) / (hcol(hi) - hcol(lo))
    return (1.0 - w) * table.T_grid[lo] + w * table.T_grid[hi]
end

"""
    snapshot_phase_setups(phase_setups, table::HelmholtzTable)

Phase-1 helper: replace each phase's `rho`/`mu`/`cp`/`k` with the values looked up
at `(T_snapshot, p_operating)`. Supercritical → both phases get identical values.
"""
function snapshot_phase_setups(phase_setups::NTuple{2, Phase}, table::HelmholtzTable)
    p = lookup_helmholtz(table, table.T_snapshot, table.p_operating)
    ph = Phase(rho = p.rho, mu = p.mu, cp = p.cp, k = p.k)
    return (ph, ph)
end

_snapshot_from_fluid_properties(phase_setups, table::HelmholtzTable) =
    snapshot_phase_setups(phase_setups, table)

# Move the lookup arrays (and grids) onto the mesh backend so the per-cell live
# update runs as a GPU kernel.
build_property(table::HelmholtzTable, mesh) = adapt(_get_backend(mesh), table)

# ===========================================================================
# Live per-cell property update — bilinear in (T, p_abs = p_operating + p_rgh)
# ===========================================================================
"""
    update_phase_properties_from_table!(phases, table, T_field, p_rgh_field, config)

Phase-2 live update. For each cell, bilinearly looks up `(rho, mu, cp, k)` at
`(T, p_operating + p_rgh)` and writes the result into both phases' per-cell
`ScalarField`s (single-phase supercritical → identical values). No-op when the
phases are not field-backed (snapshot mode).
"""
function update_phase_properties_from_table!(phases, table::HelmholtzTable, T_field, p_rgh_field, config)
    phases[1].rho isa ScalarField || return nothing

    (; backend, workgroup) = config.hardware
    ndrange = length(T_field.values)

    kernel! = _interp_phase_properties!(_setup(backend, workgroup, ndrange)...)
    kernel!(
        phases[1].rho.values, phases[2].rho.values,
        phases[1].mu.values,  phases[2].mu.values,
        phases[1].cp.values,  phases[2].cp.values,
        phases[1].k.values,   phases[2].k.values,
        T_field.values, p_rgh_field.values,
        table.T_grid, table.rho, table.mu, table.cp, table.k,
        table.p_range[1], table.p_range[2], table.n_T, table.n_p, table.p_operating,
    )
    return nothing
end

@kernel inbounds=true function _interp_phase_properties!(
        rho1, rho2, mu1, mu2, cp1, cp2, k1, k2, Ts, p_rgh,
        T_grid, rho, mu, cp, k, p_lo, p_hi, nT, np, p_op)
    i = @index(Global)

    iT_lo, iT_hi, wT = _kT_locate(T_grid, nT, Ts[i])
    jp_lo, jp_hi, wp = _kp_locate(p_lo, p_hi, np, p_op + p_rgh[i])

    r = _kbilinear(rho, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    m = _kbilinear(mu,  iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    c = _kbilinear(cp,  iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    kk = _kbilinear(k,  iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)

    rho1[i] = r; rho2[i] = r
    mu1[i]  = m; mu2[i]  = m
    cp1[i]  = c; cp2[i]  = c
    k1[i]   = kk; k2[i]  = kk
end

# No-op for any other fluid_properties type (e.g. NISTTable handles its own).
update_phase_properties_from_table!(phases, _other, T_field, p_rgh_field, config) = nothing

"""
    update_psi_drhodT_from_table!(psi_field, drhodT_field, table, T_field, p_rgh_field, config)

Interpolate the compressibility ψ=(∂ρ/∂p)_T and thermal-expansion rate (∂ρ/∂T)_p
into per-cell fields at `(T, p_operating + p_rgh)`, for the compressible pressure
equation (implicit `Time(ψ, p_rgh)` and explicit `(∂ρ/∂T)·dT/dt`).
"""
function update_psi_drhodT_from_table!(psi_field, drhodT_field, table::HelmholtzTable,
                                       T_field, p_rgh_field, config)
    (; backend, workgroup) = config.hardware
    ndrange = length(T_field.values)
    kernel! = _interp_psi_drhodT!(_setup(backend, workgroup, ndrange)...)
    kernel!(psi_field.values, drhodT_field.values,
            T_field.values, p_rgh_field.values,
            table.T_grid, table.psi, table.drhodT,
            table.p_range[1], table.p_range[2], table.n_T, table.n_p, table.p_operating)
    return nothing
end

@kernel inbounds=true function _interp_psi_drhodT!(
        psi_out, drhodT_out, Ts, p_rgh, T_grid, psi, drhodT, p_lo, p_hi, nT, np, p_op)
    i = @index(Global)
    iT_lo, iT_hi, wT = _kT_locate(T_grid, nT, Ts[i])
    jp_lo, jp_hi, wp = _kp_locate(p_lo, p_hi, np, p_op + p_rgh[i])
    psi_out[i]    = _kbilinear(psi,    iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    drhodT_out[i] = _kbilinear(drhodT, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
end

# =====================================================================
# Device-side index helpers (used inside @kernel bodies).
# =====================================================================
@inline function _kT_locate(T_grid, nT, T)
    # NaN/Inf-safe: a diverged field must clamp to a table edge, never feed a
    # garbage index to the (inbounds) matrix access — that is an illegal GPU
    # address, not a recoverable error.
    Tc = isfinite(T) ? clamp(T, T_grid[1], T_grid[nT]) : T_grid[1]
    lo = 1; hi = nT
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if T_grid[mid] <= Tc
            lo = mid
        else
            hi = mid
        end
    end
    w = (Tc - T_grid[lo]) / (T_grid[hi] - T_grid[lo])
    return lo, hi, w
end

@inline function _kp_locate(p_lo, p_hi, np, p)
    # NaN/Inf-safe (see _kT_locate): unsafe_trunc on a non-finite `pos` yields a
    # garbage column index → illegal GPU address. Clamp non-finite p to the edge.
    pc  = isfinite(p) ? clamp(p, p_lo, p_hi) : p_lo
    Δp  = (p_hi - p_lo) / (np - 1)
    pos = (pc - p_lo) / Δp
    jlo = clamp(unsafe_trunc(Int, floor(pos)) + 1, 1, np - 1)
    return jlo, jlo + 1, pos - (jlo - 1)
end

@inline function _kbilinear(A, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
    lo = (1.0 - wp) * A[iT_lo, jp_lo] + wp * A[iT_lo, jp_hi]
    hi = (1.0 - wp) * A[iT_hi, jp_lo] + wp * A[iT_hi, jp_hi]
    return (1.0 - wT) * lo + wT * hi
end

# =====================================================================
# Per-cell enthalpy <-> temperature conversion (for VariableSensibleEnthalpy).
# Both directions are pressure-aware: h = h(T, p_op + p_rgh).
# =====================================================================
"""
    interp_h_from_T!(h_field, T_field, p_rgh_field, table, config)

Per-cell sensible enthalpy from `(T, p_operating + p_rgh)` by bilinear interpolation.
"""
function interp_h_from_T!(h_field, T_field, p_rgh_field, table::HelmholtzTable, config)
    (; backend, workgroup) = config.hardware
    ndrange = length(T_field.values)
    kernel! = _interp_h_from_T!(_setup(backend, workgroup, ndrange)...)
    kernel!(h_field.values, T_field.values, p_rgh_field.values,
            table.T_grid, table.h,
            table.p_range[1], table.p_range[2], table.n_T, table.n_p, table.p_operating)
    return nothing
end

@kernel inbounds=true function _interp_h_from_T!(hs, Ts, p_rgh, T_grid, h, p_lo, p_hi, nT, np, p_op)
    i = @index(Global)
    iT_lo, iT_hi, wT = _kT_locate(T_grid, nT, Ts[i])
    jp_lo, jp_hi, wp = _kp_locate(p_lo, p_hi, np, p_op + p_rgh[i])
    hs[i] = _kbilinear(h, iT_lo, iT_hi, wT, jp_lo, jp_hi, wp)
end

"""
    invert_T_from_h!(T_field, h_field, p_rgh_field, table, config)

Per-cell temperature from `(h, p_operating + p_rgh)`: interpolate the enthalpy
column to the cell pressure, then binary-search the strictly-increasing column for
`h` (the monotone, runaway-free `T(h,p)` update). Out-of-range `h` clamps.
"""
function invert_T_from_h!(T_field, h_field, p_rgh_field, table::HelmholtzTable, config)
    (; backend, workgroup) = config.hardware
    ndrange = length(h_field.values)
    kernel! = _invert_T_from_h!(_setup(backend, workgroup, ndrange)...)
    kernel!(T_field.values, h_field.values, p_rgh_field.values,
            table.T_grid, table.h,
            table.p_range[1], table.p_range[2], table.n_T, table.n_p, table.p_operating)
    return nothing
end

@kernel inbounds=true function _invert_T_from_h!(Ts, hs, p_rgh, T_grid, h, p_lo, p_hi, nT, np, p_op)
    i = @index(Global)
    jp_lo, jp_hi, wp = _kp_locate(p_lo, p_hi, np, p_op + p_rgh[i])
    # interpolated h column at this cell's pressure (monotone increasing in T)
    h1  = (1.0 - wp) * h[1,  jp_lo] + wp * h[1,  jp_hi]
    hN  = (1.0 - wp) * h[nT, jp_lo] + wp * h[nT, jp_hi]
    hi_target = hs[i]
    if hi_target <= h1
        Ts[i] = T_grid[1]
    elseif hi_target >= hN
        Ts[i] = T_grid[nT]
    else
        lo = 1; hi = nT
        while hi - lo > 1
            mid  = (lo + hi) >>> 1
            hmid = (1.0 - wp) * h[mid, jp_lo] + wp * h[mid, jp_hi]
            if hmid <= hi_target
                lo = mid
            else
                hi = mid
            end
        end
        hlo = (1.0 - wp) * h[lo, jp_lo] + wp * h[lo, jp_hi]
        hhi = (1.0 - wp) * h[hi, jp_lo] + wp * h[hi, jp_hi]
        w = (hi_target - hlo) / (hhi - hlo)
        Ts[i] = (1.0 - w) * T_grid[lo] + w * T_grid[hi]
    end
end
