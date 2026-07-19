# =============================================================================
# V&V platform — shared helpers
# =============================================================================

using Statistics

vv_grid(name) = joinpath(pkgdir(XCALibre, "examples/0_GRIDS"), name)

l2(err)   = sqrt(sum(abs2, err) / length(err))
linf(err) = maximum(abs, err)

# Observed order of accuracy from errors on two meshes/dts refined by `r` (default 2).
observed_order(e_coarse, e_fine; r=2) = log(e_coarse / e_fine) / log(r)

"""
    vv_check(name, value; below=nothing, between=nothing) -> Bool

Quantitative check with a diagnosable log line: prints the measured value next to
its bound so a red test is debuggable from the log alone. Returns the Bool for
use inside `@test`.
"""
function vv_check(name, value; below=nothing, between=nothing)
    ok = true
    if below !== nothing
        ok &= value < below
        @info "V&V check" name value bound=below pass=(value < below)
    end
    if between !== nothing
        lo, hi = between
        ok &= lo < value < hi
        @info "V&V check" name value range=(lo, hi) pass=(lo < value < hi)
    end
    return ok
end

# Max velocity magnitude over the field (CPU copy — test-side only).
function max_umag(U)
    ux = Array(U.x.values); uy = Array(U.y.values); uz = Array(U.z.values)
    return maximum(sqrt.(ux .^ 2 .+ uy .^ 2 .+ uz .^ 2))
end

# ---------------------------------------------------------------------------
# Structured-quad column statistics (front tracking for advection/settling).
# ---------------------------------------------------------------------------

"""
    column_means(mesh, vals; axis=1) -> (coords, means)

Group cells by their centre coordinate along `axis` (1=x, 2=y) and return the
sorted unique coordinates with the mean of `vals` in each column/row.
"""
function column_means(mesh, vals; axis=1)
    coords = [c.centre[axis] for c in mesh.cells]
    ks = round.(coords; digits=9)
    uq = sort(unique(ks))
    means = [mean(vals[ks .== k]) for k in uq]
    return uq, means
end

# ---------------------------------------------------------------------------
# Boundary-patch audits (CPU mesh) — power the conservation/balance contracts.
# ---------------------------------------------------------------------------

function patch_range(mesh, patch::Symbol)
    for b in mesh.boundaries
        b.name == patch && return b.IDs_range
    end
    error("V&V: boundary patch :$patch not found in mesh")
end

patch_area(mesh, patch::Symbol) = sum(mesh.faces[i].area for i in patch_range(mesh, patch))

"Area-weighted mean of owner-cell values over a boundary patch."
function patch_mean(mesh, patch::Symbol, vals)
    s = 0.0; A = 0.0
    for i in patch_range(mesh, patch)
        f = mesh.faces[i]
        s += f.area * vals[f.ownerCells[1]]
        A += f.area
    end
    return s / A
end

"""
    front_position(coords, means, level) -> x

Position where the column-mean profile crosses `level`, scanning from the
high-value side: last coordinate with mean ≥ level, linearly interpolated
toward the next column. NaN if the profile never reaches `level`.
"""
function front_position(coords, means, level)
    i = findlast(>=(level), means)
    i === nothing && return NaN
    i == length(means) && return coords[end]
    x1, x2 = coords[i], coords[i+1]
    m1, m2 = means[i], means[i+1]
    m2 == m1 && return x1
    return x1 + (level - m1) / (m2 - m1) * (x2 - x1)
end
