export save_checkpoint, load_checkpoint!

"""
    save_checkpoint(filename; time=0.0, fields...)

Write a JLD2 snapshot of the cell-centre values of the given fields.
Each keyword pair must name an `AbstractScalarField` or `AbstractVectorField`
from the live model; vector fields are stored as three components
(`name__x`, `name__y`, `name__z`). GPU arrays are copied to host before writing.

# Example

```julia
save_checkpoint("snap.jld2"; time = 4.0,
    U     = model.momentum.U,
    p_rgh = model.fluid.p_rgh,
    alpha = model.fluid.alpha,
    T     = model.energy.T,
    k     = model.turbulence.k,
    omega = model.turbulence.omega,
    nut   = model.turbulence.nut,
)
```
"""
function save_checkpoint(filename::AbstractString; time::Real=0.0, fields...)
    JLD2.jldopen(filename, "w") do io
        io["time"] = Float64(time)
        for (name, fld) in pairs(fields)
            _save_field(io, String(name), fld)
        end
    end
    return filename
end

_save_field(io, name, f::AbstractScalarField) = (io[name] = Array(f.values))

function _save_field(io, name, f::AbstractVectorField)
    io[name * "__x"] = Array(f.x.values)
    io[name * "__y"] = Array(f.y.values)
    io[name * "__z"] = Array(f.z.values)
end

"""
    load_checkpoint!(filename; fields...) -> time::Float64

Reload field values from a `save_checkpoint` snapshot into the supplied
live fields. Each keyword field must match the type and size of what was
saved. Returns the simulation time recorded in the snapshot (or `0.0` if
absent). Works for both CPU and GPU fields — `copyto!` performs the
host→device transfer when the live field lives on a GPU.

# Example

```julia
t0 = load_checkpoint!("snap.jld2";
    U     = model.momentum.U,
    p_rgh = model.fluid.p_rgh,
    alpha = model.fluid.alpha,
    T     = model.energy.T,
    k     = model.turbulence.k,
    omega = model.turbulence.omega,
    nut   = model.turbulence.nut,
)
```
"""
function load_checkpoint!(filename::AbstractString; fields...)
    return JLD2.jldopen(filename, "r") do io
        for (name, fld) in pairs(fields)
            _load_field!(io, String(name), fld)
        end
        haskey(io, "time") ? Float64(io["time"]) : 0.0
    end
end

function _check_size(name, n_disk, n_live)
    n_disk == n_live || throw(ArgumentError(
        "Checkpoint size mismatch on field `$name`: file has $n_disk cells, " *
        "model has $n_live cells. The checkpoint was saved on a different " *
        "mesh — regenerate it (rerun prep_heatflux.jl on the current mesh)."))
end

function _load_field!(io, name, f::AbstractScalarField)
    arr = io[name]
    _check_size(name, length(arr), length(f.values))
    copyto!(f.values, arr)
    return nothing
end

function _load_field!(io, name, f::AbstractVectorField)
    for (suffix, comp) in (("__x", f.x), ("__y", f.y), ("__z", f.z))
        arr = io[name * suffix]
        _check_size(name * suffix, length(arr), length(comp.values))
        copyto!(comp.values, arr)
    end
    return nothing
end
