module MeshCompare
using Test, StaticArrays, Statistics
export compare_meshes2, face_key

"""
    face_key(nodes::AbstractVector{<:Integer})

Return a canonical NTuple of the node ids that define a face, sorted
once so it is order‑independent and hash‑able.
"""
face_key(nodes) = ntuple(i->nodes[sortperm(nodes)[i]], length(nodes))

"""
    compare_meshes(unv::UNV3D_mesh, med::MED3D_mesh;
                   atol = 1e-9, rtol = 1e-6)

Run a battery of equivalence checks between a UNV and MED mesh.
Returns a tuple `(pass::Bool, report::Vector{String})`.
"""
function compare_meshes2(unv, med; atol=1e-9, rtol=1e-6)
    failures = String[]

    ## A. Top‑level counts
    if length(unv.nodes) != length(med.nodes)
        push!(failures,
              "Node count: $(length(unv.nodes)) vs $(length(med.nodes))")
    end
    if length(unv.cells) != length(med.cells)
        push!(failures,
              "Cell count: $(length(unv.cells)) vs $(length(med.cells))")
    end
    if length(unv.faces) != length(med.faces)
        push!(failures,
              "Face count: $(length(unv.faces)) vs $(length(med.faces))")
    end

    ## B. Node coordinates
    @inbounds for (i,(p₁,p₂)) in enumerate(zip(unv.nodes, med.nodes))
        isapprox(p₁.coords, p₂.coords; atol, rtol) && continue
        push!(failures, "Node $i coords differ → $(p₁.coords) vs $(p₂.coords)")
    end

    ## C. Cell‑node connectivity (fast hash)
    cell_hash(v) = hash(v.nodes_range) ⊻ hash(copy(v.node_ids))
    h_unv = map(cell_hash, unv.cells)
    h_med = map(cell_hash, med.cells)
    h_unv == h_med ||
        push!(failures, "Cell‑node connectivity differs (hash mismatch)")

    ## D+E. Face dictionary
    unv_face_idx = Dict(face_key(nodes(unv, f)) => i for (i, f) in enumerate(unv.faces))
    med_to_unv = Dict{Int,Int}()

    for (i,f) in enumerate(med.faces)
        key = face_key(nodes(med, f))
        if !haskey(unv_face_idx, key)
            push!(failures, "MED face $i unmatched")
        else
            med_to_unv[i] = unv_face_idx[key]
        end
    end
    if length(med_to_unv) != length(unv.faces)
        push!(failures, "Not all UNV faces were mapped")
    end

    ## E. Per‑face numeric props
    for (mid,uid) in pairs(med_to_unv)
        fm, fu = med.faces[mid], unv.faces[uid]
        sort(fm.ownerCells) == sort(fu.ownerCells) ||
            push!(failures, "Face $mid owners differ")
        isapprox(fm.centre, fu.centre; atol, rtol) ||
            push!(failures, "Face $mid centre differs")
        isapprox(fm.area, fu.area; rtol)                 ||
            push!(failures, "Face $mid area differs")
        n_ok = isapprox(fm.normal, fu.normal; atol, rtol) ||
               isapprox(fm.normal, -fu.normal; atol, rtol)
        n_ok || push!(failures, "Face $mid normal differs (sign‑insensitive)")
    end

    ## F. Cell‑face connectivity
    for (cid,(c1,c2)) in enumerate(zip(unv.cells, med.cells))
        faces1 = map(med_to_unv, c2.face_ids)
        Set(faces1) == Set(c1.face_ids) ||
            push!(failures, "Cell $cid face set differs")
        c1.nsign == c2.nsign ||
            push!(failures, "Cell $cid nsign differs")
    end

    ## G. Cell volumes & centroids
    for (i,(c₁,c₂)) in enumerate(zip(unv.cells, med.cells))
        isapprox(c₁.centre, c₂.centre; atol, rtol) ||
            push!(failures, "Cell $i centre differs")
        isapprox(c₁.volume, c₂.volume; rtol) ||
            push!(failures, "Cell $i volume differs")
    end

    ## H. Boundaries
    bnd_unv = filter(f->f.ownerCells[2]==0, unv.faces)
    bnd_med = filter(f->f.ownerCells[2]==0, med.faces)
    length(bnd_unv) == length(bnd_med) ||
        push!(failures, "Boundary face count differs")
    # Further boundary patch checks...

    return isempty(failures), failures
end

end # module MeshCompare
