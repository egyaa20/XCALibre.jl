export NISTTable, NISTTable_N2_3p5MPa

using KernelAbstractions

# =============================================================================
# NISTTable — hardcoded isobaric NIST data, drop-in replacement for
# HelmholtzTable when you want to bypass the Helmholtz EOS computation and
# work directly with tabulated reference data.
#
# Unlike HelmholtzTable (which holds separate liquid/vapour branches), the
# tabulated data here is for a SINGLE supercritical branch — both phases
# in the Multiphase pipeline get the same property values, so the case
# behaves as effectively single-phase.
#
# Out-of-range temperatures are clamped (same as HelmholtzTable), but a
# WARNING is emitted on first detection of clamping per session — so you
# can tell when the table is no longer valid.
# =============================================================================

"""
    NISTTable

Hardcoded isobaric NIST table.

Fields:
- `pressure`  : table pressure [Pa]
- `T_range`   : (T_min, T_max) [K]
- `n_points`  : number of grid points
- `T_snapshot`: representative T for `Phase` snapshot construction [K]
- `T_grid`    : T grid array
- `rho`, `mu`, `cp`, `k` : property arrays at each grid point (single branch)

Methods follow the same `snapshot_phase_setups` / `build_property` /
`update_phase_properties_from_table!` interface as `HelmholtzTable`.
"""
struct NISTTable{V}
    pressure::Float64
    T_range::Tuple{Float64, Float64}
    n_points::Int
    T_snapshot::Float64
    T_grid::V
    rho::V
    mu::V
    cp::V
    k::V
end
Adapt.@adapt_structure NISTTable

# Warn-once-per-session flag for out-of-range detection.
const _nist_oor_warned = Ref(false)

# -----------------------------------------------------------------------------
# Hardcoded data: N₂ at P = 3.5 MPa, T ∈ [75, 100] K (NIST REFPROP isobar)
#
# Subsampled from the user-supplied 0.1-K table at 1-K spacing (26 points).
# Properties are in SI:
#   rho [kg/m³], mu [Pa·s], cp [J/(kg·K)], k [W/(m·K)]
# -----------------------------------------------------------------------------

const NIST_N2_3p5MPa_T   = [75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0,
                            83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0,
                            91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0,
                            99.0, 100.0]

const NIST_N2_3p5MPa_RHO = [824.60, 820.34, 816.05, 811.73, 807.38, 802.99,
                            798.57, 794.11, 789.61, 785.07, 780.49, 775.87,
                            771.20, 766.48, 761.71, 756.89, 752.02, 747.08,
                            742.09, 737.03, 731.90, 726.70, 721.42, 716.06,
                            710.62, 705.08]

const NIST_N2_3p5MPa_CP  = [2000.5, 2003.1, 2006.0, 2009.1, 2012.5, 2016.1,
                            2020.1, 2024.3, 2028.9, 2033.9, 2039.2, 2044.9,
                            2051.1, 2057.7, 2064.8, 2072.4, 2080.6, 2089.3,
                            2098.7, 2108.9, 2119.7, 2131.4, 2144.0, 2157.5,
                            2172.1, 2187.9]

const NIST_N2_3p5MPa_MU  = [1.8624e-4, 1.7887e-4, 1.7194e-4, 1.6541e-4,
                            1.5926e-4, 1.5346e-4, 1.4797e-4, 1.4279e-4,
                            1.3787e-4, 1.3321e-4, 1.2878e-4, 1.2457e-4,
                            1.2056e-4, 1.1674e-4, 1.1309e-4, 1.0960e-4,
                            1.0627e-4, 1.0307e-4, 1.0001e-4, 9.7066e-5,
                            9.4238e-5, 9.1516e-5, 8.8892e-5, 8.6360e-5,
                            8.3914e-5, 8.1548e-5]

const NIST_N2_3p5MPa_K   = [0.15313, 0.15120, 0.14928, 0.14736, 0.14545,
                            0.14353, 0.14162, 0.13971, 0.13780, 0.13589,
                            0.13398, 0.13207, 0.13016, 0.12826, 0.12636,
                            0.12446, 0.12257, 0.12067, 0.11878, 0.11688,
                            0.11498, 0.11309, 0.11119, 0.10929, 0.10739,
                            0.10548]

"""
    NISTTable_N2_3p5MPa(; T_snapshot = 87.5) -> NISTTable

Pre-built NIST table for nitrogen at 3.5 MPa, T ∈ [75, 100] K (26 pts).
"""
function NISTTable_N2_3p5MPa(; T_snapshot::Real = 87.5)
    T_lo, T_hi = NIST_N2_3p5MPa_T[1], NIST_N2_3p5MPa_T[end]
    T_lo <= T_snapshot <= T_hi ||
        error("NISTTable_N2_3p5MPa: T_snapshot=$T_snapshot outside ($T_lo, $T_hi) K")
    _nist_oor_warned[] = false   # reset flag for new session
    @info "NISTTable: hardcoded N₂ @ 3.5 MPa, $(length(NIST_N2_3p5MPa_T)) pts, T ∈ ($T_lo, $T_hi) K"
    return NISTTable(
        3.5e6,
        (T_lo, T_hi),
        length(NIST_N2_3p5MPa_T),
        float(T_snapshot),
        NIST_N2_3p5MPa_T,
        NIST_N2_3p5MPa_RHO,
        NIST_N2_3p5MPa_MU,
        NIST_N2_3p5MPa_CP,
        NIST_N2_3p5MPa_K,
    )
end

"""
    lookup_nist(table::NISTTable, T)

Linear interpolation of `(rho, mu, cp, k)` at temperature `T`.
Out-of-range temperatures clamp.
"""
function lookup_nist(table::NISTTable, T::Real)
    T_lo, T_hi = table.T_range
    Tf = float(T)
    T_clamped = clamp(Tf, T_lo, T_hi)
    Δ = (T_hi - T_lo) / (table.n_points - 1)
    pos = (T_clamped - T_lo) / Δ
    i_lo = clamp(Int(floor(pos)) + 1, 1, table.n_points - 1)
    i_hi = i_lo + 1
    w = pos - (i_lo - 1)
    @inline interp(a) = (1.0 - w) * a[i_lo] + w * a[i_hi]
    return (rho = interp(table.rho), mu = interp(table.mu),
            cp  = interp(table.cp),  k  = interp(table.k))
end

# -----------------------------------------------------------------------------
# Multiphase pipeline integration (mirror of HelmholtzTable hooks)
# -----------------------------------------------------------------------------

"""
    snapshot_phase_setups(phase_setups::NTuple{2,Phase}, table::NISTTable)

Both phases get the SAME property snapshot at `table.T_snapshot` —
the case behaves as effectively single-phase supercritical fluid.
"""
function snapshot_phase_setups(phase_setups::NTuple{2, Phase}, table::NISTTable)
    p = lookup_nist(table, table.T_snapshot)
    phase_set = Phase(rho = p.rho, mu = p.mu, cp = p.cp, k = p.k)
    return (phase_set, phase_set)
end

_snapshot_from_fluid_properties(phase_setups, table::NISTTable) =
    snapshot_phase_setups(phase_setups, table)

build_property(table::NISTTable, mesh) = adapt(_get_backend(mesh), table)

"""
    update_phase_properties_from_table!(phases, table::NISTTable, T_field, config)

Per-cell property update. Both phases receive identical values from the
single-branch NIST table. Emits a warning once per session if any cell T
leaves `T_range`.
"""
function update_phase_properties_from_table!(phases, table::NISTTable, T_field, config)
    phases[1].rho isa ScalarField || return nothing
    # Out-of-range monitor (cheap host-side reduction over GPU array)
    if !_nist_oor_warned[]
        T_max = maximum(T_field.values)
        T_min = minimum(T_field.values)
        if T_max > table.T_range[2] || T_min < table.T_range[1]
            @warn "NISTTable: T outside table range — clamping! " *
                  "T ∈ ($T_min, $T_max) K vs table range $(table.T_range) K. " *
                  "Subsequent clampings will be silent."
            _nist_oor_warned[] = true
        end
    end
    (; backend, workgroup) = config.hardware
    Ts = T_field.values
    ndrange = length(Ts)
    kernel! = _interp_nist_properties!(_setup(backend, workgroup, ndrange)...)
    kernel!(
        phases[1].rho.values, phases[2].rho.values,
        phases[1].mu.values,  phases[2].mu.values,
        phases[1].cp.values,  phases[2].cp.values,
        phases[1].k.values,   phases[2].k.values,
        Ts,
        table.rho, table.mu, table.cp, table.k,
        table.T_range[1], table.T_range[2], table.n_points,
    )
    return nothing
end

@kernel inbounds=true function _interp_nist_properties!(
        rho1, rho2, mu1, mu2, cp1, cp2, k1, k2, Ts,
        rho_tab, mu_tab, cp_tab, k_tab, T_lo, T_hi, npts)
    i = @index(Global)
    Δ = (T_hi - T_lo) / (npts - 1)
    T_clamped = clamp(Ts[i], T_lo, T_hi)
    pos  = (T_clamped - T_lo) / Δ
    i_lo = clamp(unsafe_trunc(Int, floor(pos)) + 1, 1, npts - 1)
    i_hi = i_lo + 1
    w    = pos - (i_lo - 1)
    w1   = 1.0 - w
    # Identical assignment to both phases (single-branch supercritical fluid)
    r = w1*rho_tab[i_lo] + w*rho_tab[i_hi]
    m = w1*mu_tab[i_lo]  + w*mu_tab[i_hi]
    c = w1*cp_tab[i_lo]  + w*cp_tab[i_hi]
    κ = w1*k_tab[i_lo]   + w*k_tab[i_hi]
    rho1[i] = r;  rho2[i] = r
    mu1[i]  = m;  mu2[i]  = m
    cp1[i]  = c;  cp2[i]  = c
    k1[i]   = κ;  k2[i]   = κ
end
