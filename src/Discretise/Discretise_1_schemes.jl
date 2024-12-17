export scheme!, scheme_source!

#= NOTE:
In source scheme the following indices are used and should be used with care:
cID - Index of the cell outer loop. Use to index "b" 
cIndex - Index of the cell based on sparse matrix. Use to index "nzval_array"
=#


function calculate_distance(point1::SVector{N, T}, point2::SVector{N, T}, method::String = "euclidean") where {N, T <: Real}
    if method == "euclidean"
        return norm(point1 - point2)
    elseif method == "manhattan"
        return sum(abs.(point1 - point2)) 
    elseif method == "chebysev"
        return maximum(abs.(point1 - point2))
    elseif method == "euclidean_squared"
        return dot(point1 - point2, point1 - point2)

        # return sum((point1 - point2).^2)*0.5
    else
        error("Unsupported method: $method. Use \"euclidean\" or \"manhattan\".")
    end
end



# TIME 

                # SteadyState
@inline function scheme!(
    term::Operator{F,P,I,Time{SteadyState}}, 
    nzval_array, cell, face,  cellN, ns, cIndex, nIndex, fID, prev, runtime)  where {F,P,I}
    # nothing
    0.0, 0.0 # add types if this approach works
end
@inline scheme_source!(
    term::Operator{F,P,I,Time{SteadyState}}, cell, cID, cIndex, prev, runtime)  where {F,P,I} = begin
    0.0, 0.0
end

                    ## Euler
@inline function scheme!(
    term::Operator{F,P,I,Time{Euler}}, 
    nzval_array, cell, face,  cellN, ns, cIndex, nIndex, fID, prev, runtime)  where {F,P,I}
    # nothing
    0.0, 0.0 # add types if this approach works
end
@inline scheme_source!(
    term::Operator{F,P,I,Time{Euler}}, cell, cID, cIndex, prev, runtime)  where {F,P,I} = begin
        volume = cell.volume
        vol_rdt = volume/runtime.dt
        
        # Increment sparse and b arrays 
        ac = vol_rdt
        b = prev[cID]*vol_rdt
        return ac, b
end

# LAPLACIAN

@inline function scheme!(
    term::Operator{F,P,I,Laplacian{Linear}}, 
    nzval_array, cell, face,  cellN, ns, cIndex, nIndex, fID, prev, runtime
    )  where {F,P,I}

    (; area, normal, delta, e) = face # maybe define delta differently e.g. not perpendicular across the face but somehow else
    Sf = ns*area*normal
    e = ns*e
    Ef = ((Sf⋅Sf)/(Sf⋅e))*e
    Ef_mag = norm(Ef)
    ap = term.sign*(term.flux[fID] * Ef_mag)/delta

    # ap = term.sign*(term.flux[fID] * area)/delta
    # Increment sparse array
    ac = -ap
    an = ap
    return ac, an
end
@inline scheme_source!(
    term::Operator{F,P,I,Laplacian{Linear}}, cell, cID, cIndex, prev, runtime)  where {F,P,I} = begin
    0.0, 0.0
end

# DIVERGENCE

# Linear
@inline function scheme!(
    term::Operator{F,P,I,Divergence{Linear}}, 
    nzval_array, cell, face, cellN, ns, cIndex, nIndex, fID, prev, runtime
    )  where {F,P,I}

    # Retrieve mesh centre values
    xf = face.centre
    xC = cell.centre
    xN = cellN.centre

    # [weighted scheme]
    # Calculate weights using normal functions
    # weight = norm(xf - xC)/norm(xN - xC) 
    # one_minus_weight = one(eltype(weight)) - weight 

    # # Calculate required increment
    # ap = term.sign*(term.flux[fID]*ns)
    # ac = ap*one_minus_weight
    # an = ap*weight
    # return ac, an


    # Calculate weights using normal functions
    # weight = calculate_distance(xf, xC, "euclidean")/calculate_distance(xN, xC, "euclidean")

    weight = face.weight
    one_minus_weight = one(eltype(weight)) - weight 

    # Calculate required increment
    ap = term.sign*(term.flux[fID]*ns)
    ac = ap*one_minus_weight
    an = ap*weight
    return ac, an
end
@inline scheme_source!(
    term::Operator{F,P,I,Divergence{Linear}}, cell, cID, cIndex, prev, runtime) where {F,P,I} = begin
    0.0, 0.0
end

# Upwind
@inline function scheme!(
    term::Operator{F,P,I,Divergence{Upwind}}, 
    nzval_array, cell, face, cellN, ns, cIndex, nIndex, fID, prev, runtime
    )  where {F,P,I}
    # Calculate required increment
    println("[SCHEME NAME] > Upwind")
    ap = term.sign*(term.flux[fID]*ns)
    ac = max(ap, 0.0) 
    an = -max(-ap, 0.0)
    return ac, an
end
@inline scheme_source!(
    term::Operator{F,P,I,Divergence{Upwind}}, cell, cID, cIndex, prev, runtime) where {F,P,I} = begin
    0.0, 0.0
end

# LUST
@inline function scheme!(
    term::Operator{F,P,I,Divergence{LUST}}, 
    nzval_array, cell, face, cellN, ns, cIndex, nIndex, fID, prev, runtime
    )  where {F,P,I}
    # Retrieve mesh centre values
    println("[SCHEME NAME] > LUST")
    xf = face.centre
    xC = cell.centre
    xN = cellN.centre
    
    # Calculate weights using normal functions
    weight = norm(xf - xC)/norm(xN - xC)
    one_minus_weight = one(eltype(weight)) - weight

    # Calculate coefficients
    ap = term.sign*(term.flux[fID]*ns)
    acLinear = ap*one_minus_weight
    anLinear = ap*weight
    acUpwind = max(ap, 0.0) 
    anUpwind = -max(-ap, 0.0)
    ac = 0.75*acLinear + 0.25*acUpwind
    an = 0.75*anLinear + 0.25*anUpwind
    return ac, an
end
@inline scheme_source!(
    term::Operator{F,P,I,Divergence{LUST}}, cell, cID, cIndex, prev, runtime) where {F,P,I} = begin
    0.0, 0.0
end

# BoundedUpwind
@inline function scheme!(
    
    term::Operator{F,P,I,Divergence{BoundedUpwind}}, 
    nzval_array, cell, face, cellN, ns, cIndex, nIndex, fID, prev, runtime
    )  where {F,P,I}
    # Calculate required increment
    println("[SCHEME NAME] > BoundedUpwind")

    volume = cell.volume
    ap = term.sign*(term.flux[fID]*ns)
    ac = max(ap, 0.0) - term.flux[fID]#*volume 
    an = -max(-ap, 0.0)
    return ac, an
end
@inline scheme_source!(
    term::Operator{F,P,I,Divergence{BoundedUpwind}}, cell, cID, cIndex, prev, runtime) where {F,P,I} = begin
    0.0, 0.0
end


# IMPLICIT SOURCE
@inline function scheme!(
    term::Operator{F,P,I,Si}, 
    nzval_array, cell, face,  cellN, ns, cIndex, nIndex, fID, prev, runtime
    )  where {F,P,I}
    0.0, 0.0
end
@inline scheme_source!(
    term::Operator{F,P,I,Si}, cell, cID, cIndex, prev, runtime)  where {F,P,I} = begin
    println("[SCHEME NAME] > IMPLICIT SOURCE")

    # Retrieve and calculate flux for cell 
    flux = term.sign*term.flux[cID]*cell.volume # indexed with cID
    ac = flux # indexed with cIndex
    ac, 0.0
end
