using Revise
Revise.track("src/MED/MED.jl")
include("../src/MED/MED.jl")

include("meshChecker.jl")
include("meshQuality.jl")
include("mesh_comparison.jl")
# includet("../src/MED/MED.jl")

# Revise.track("src/MED/MED_2_builder.jl")
# Revise.track("src/MED/MED_1_reader.jl")

using XCALibre
using LinearAlgebra
using StaticArrays
using Statistics #to compute mean values
# using XCALibre.MED


# using HDF5


using .meshCheck
using .meshQuality
using .MED
using Test, .MeshCompare


grids_dir = pkgdir(XCALibre, "examples/testing_grids")

grid_med = "med2.med"
# grid_med = "multiblock.med"
# grid_med = "no_auto_create.med"
# grid_med = "med.med"
# grid_unv = "multiblock.unv"
grid_unv = "med2.unv"
# grid_unv = "mesh_for_testing2.unv"
# grid_unv = "block_unv.unv"

mesh_med = joinpath(grids_dir, grid_med)
mesh_unv = joinpath(grids_dir, grid_unv)

# h5open(mesh_med,"r") do f
#     dset = f["/ENS_MAA/Mesh_1/-0000000000000000001-0000000000000000001/MAI/TE4/NOD"]
#     raw = read(dset)                     # Vector{Int64} of length 88
#     n_per_elem = 4
#     n_elems    = length(raw) ÷ n_per_elem  # should be 22
#     mat = reshape(raw, n_per_elem, n_elems)
#     @show size(mat)      # (4, 22)
#     @show mat[:,1]       # first tet’s 4 node IDs
#     @show mat[:,2]       # second tet’s 4 node IDs
#   end
  
med = MED3D_mesh(mesh_med, scale=0.001)
# unv = UNV3D_mesh(mesh_unv, scale=0.001)





# function match_faces_by_centroid(med, unv; tol = 1e-8)
#     matched = []
#     unmatched_med = []
#     unmatched_unv = collect(1:length(unv.faces))

#     for (i, f_med) in enumerate(med.faces)
#         found = false
#         for (j, f_unv) in enumerate(unv.faces)
#             if norm(f_med.centre .- f_unv.centre) < tol
#                 push!(matched, (i, j))
#                 deleteat!(unmatched_unv, findfirst(==(j), unmatched_unv))
#                 found = true
#                 break
#             end
#         end
#         if !found
#             push!(unmatched_med, i)
#         end
#     end

#     println("\nStarting matching by face centre with tolerance = $(tol)\n")
#     println("✅ Matched Faces: $(length(matched))")
#     println("❌ Unmatched MED Faces: $(length(unmatched_med))")
#     println("❌ Unmatched UNV Faces: $(length(unmatched_unv))")

#     for (m, u) in matched
#         fn_med = med.face_nodes[med.faces[m].nodes_range]
#         fn_unv = unv.face_nodes[unv.faces[u].nodes_range]
#         println("MED[$m] ↔ UNV[$u]")
#         println("  MED Nodes: $fn_med")
#         println("  UNV Nodes: $fn_unv")
#         println("  Match: ", fn_med == fn_unv, "\n")
#     end

#     return matched, unmatched_med, unmatched_unv
# end


# function match_cells_by_centroid(med, unv; tol=1e-8)

#     # Build KD-tree or brute-force matching by comparing coordinates
#     matched_pairs = Dict{Int, Int}()
#     unmatched_med = Set(1:length(med.cells))
#     unmatched_unv = Set(1:length(unv.cells))

#     # Build lookup from unv cell centres → cell index
#     unv_centres = [unv.cells[i].centre for i in 1:length(unv.cells)]
#     unv_lookup = [(i, c) for (i, c) in enumerate(unv_centres)]

#     println("Starting matching by centroid with tolerance = $tol")

#     for (med_idx, med_cell) in enumerate(med.cells)
#         med_centre = med_cell.centre
#         # Try to find a UNV cell with close enough centre
#         found = false
#         for (unv_idx, unv_centre) in unv_lookup
#             if norm(med_centre - unv_centre) < tol
#                 matched_pairs[med_idx] = unv_idx
#                 delete!(unmatched_med, med_idx)
#                 delete!(unmatched_unv, unv_idx)
#                 found = true
#                 break
#             end
#         end
#         if !found
#             println("⚠️  No match found for MED cell $med_idx with centre $med_centre")
#         end
#     end

#     println("\n✅ Matched Cells: $(length(matched_pairs))")
#     println("❌ Unmatched MED Cells: $(length(unmatched_med))")
#     println("❌ Unmatched UNV Cells: $(length(unmatched_unv))")

#     return matched_pairs, unmatched_med, unmatched_unv
# end



# matches, missing_med, missing_unv = match_faces_by_delta(med, unv, tol=1e-8)  #match_cells_by_centroid

# # Check nodeID match for a few pairs:
# for (med_id, unv_id) in Iterators.take(matches, 5)
#     med_nodes = med.cell_nodes[med.cells[med_id].nodes_range]
#     unv_nodes = unv.cell_nodes[unv.cells[unv_id].nodes_range]
#     println("MED[$med_id] ↔ UNV[$unv_id]")
#     println("  MED Nodes: $med_nodes")
#     println("  UNV Nodes: $unv_nodes")
#     println("  Match: ", med_nodes == unv_nodes)
# end


# p_unv = unv.cells[5].nodes_range
# p_med = med.cells[5].nodes_range

# nodes_unv = unv.nodes[p_unv]

# # 1) Inspect the two ranges
# fr_med =   med.cells[5].faces_range   # e.g.  A:B
# fr_unv =   unv.cells[2].faces_range   # e.g.  C:D
# println("MED faces_range for cell 5: ", fr_med)
# println("UNV faces_range for cell 2: ", fr_unv)

# # 2) Pull out the actual face IDs
# med_fids = med.cell_faces[fr_med]     # these are indices into med.faces
# unv_fids = unv.cell_faces[fr_unv]     # indices into unv.faces
# println("MED face IDs: ", med_fids)
# println("UNV face IDs: ", unv_fids)

# # 3) Compare as sets
# println("Only in MED: ", setdiff(med_fids, unv_fids))
# println("Only in UNV: ", setdiff(unv_fids, med_fids))

# # 4) For each mismatched face, print a bit more:
# for fid in setdiff(med_fids, unv_fids)
#     println(" MED-only fid=", fid,
#             " centre=", med.faces[fid].centre,
#             " nodes=", collect(med.face_nodes[ med.faces[fid].nodes_range ]))
# end
# for fid in setdiff(unv_fids, med_fids)
#     println(" UNV-only fid=", fid,
#             " centre=", unv.faces[fid].centre,
#             " nodes=", collect(unv.face_nodes[ unv.faces[fid].nodes_range ]))
# end



# function verify_meshs(med::Mesh3, unv::Mesh3)
#     ## 1) build the cell→cell map by centroid
#     med_centres = [c.centre for c in med.cells]
#     unv_centres = [c.centre for c in unv.cells]
#     mapping = Dict{Int,Int}()
#     for (i, cm) in enumerate(med_centres)
#         # build a 1-D Vector of distances
#         dists = [norm(cm - cu) for cu in unv_centres]
#         mapping[i] = argmin(dists)
#     end

#     ## 2) cell‐level checks inside that same loop
#     for (i,j) in mapping
#         # a) connectivity
#         conn_med = sort(med.cell_nodes[ med.cells[i].nodes_range ])
#         conn_unv = sort(unv.cell_nodes[ unv.cells[j].nodes_range ])
#         @assert conn_med == conn_unv "Cell $i↔$j connectivity mismatch"

#         # b) faces_range
#         fr_med = med.cells[i].faces_range
#         fr_unv = unv.cells[j].faces_range
#         @assert fr_med == fr_unv "Cell $i↔$j face‐range mismatch"
#     end
#     println("✅ All cell checks passed")

#     ## 3) build the face→face map by centroid
#     med_fc = [mean(med.nodes[n].coords for n in med.face_nodes[f.nodes_range]) for f in med.faces]
#     unv_fc = [mean(unv.nodes[n].coords for n in unv.face_nodes[f.nodes_range]) for f in unv.faces]
#     face_map = Dict(i => argmin(norm.(med_fc[i] .- unv_fc)) for i in eachindex(med_fc))

#     ## 4) face‐level checks in its own loop
#     for (k,l) in face_map
#         # a) normals
#         dn = dot(med.faces[k].normal, unv.faces[l].normal)
#         @assert abs(abs(dn) - 1.0) < 1e-8 "Face $k↔$l normal mismatch"

#         # b) ownerCells
#         @assert med.faces[k].ownerCells == unv.faces[l].ownerCells "Face $k↔$l ownerCells mismatch"
#     end
#     println("✅ All face checks passed")

#     ## 5) boundary patches
#     for n in eachindex(med.boundaries)
#         bm = med.boundaries[n]
#         bu = unv.boundaries[n]
#         med_ids = med.cell_faces[bm.IDs_range]
#         unv_ids = unv.cell_faces[bu.IDs_range]
#         mapped_ids = Set(face_map[f] for f in med_ids)
#         @assert mapped_ids == Set(unv_ids) "Boundary $(bm.name) mismatch"
#         println("✅ Boundary ", bm.name, " matches")
#     end

#     println("🎉 All verifications passed!")
# end


# verify_meshs(med, unv)





# med_centres = [c.centre for c in med.cells]
# med_vols    = [c.volume for c in med.cells]

# unv_centres = [c.centre for c in unv.cells]
# unv_vols    = [c.volume  for c in unv.cells]

# sort(med_vols) == sort(unv_vols)         # volumes agree?
# maximum(abs.(sort(med_vols) .- sort(unv_vols)))  # largest volume error

# mapping = Dict{Int,Int}()
# for (i,cm) in enumerate(med_centres)
#     dists = [norm(cm - cu) for cu in unv_centres]
#     j = argmin(dists)
#     mapping[i] = j
# end
# # now mapping[i] tells you which unv.cells[j] is closest to med.cells[i]


# for (i,j) in mapping
#     println("MED $i → UNV $j, Δcentre = ", norm(med_centres[i] - unv_centres[j]))
# end





# for (i,j) in mapping
#     med_conn = sort(med.cells[i].nodes)
#     unv_conn = sort(unv.cells[j].nodes)
#     println("Cell $i vs $j conn equal? ", med_conn == unv_conn)
# end


report_mesh_quality(unv)
compute_mesh_quality_stats(unv)
general_mesh_check(1*0.2*0.2, unv)

println("<<<-------------------------------->>>")
println("<<<-------------------------------->>>")


report_mesh_quality(med)
compute_mesh_quality_stats(med)
general_mesh_check(1*0.2*0.2, med)


# print_element_counts(mesh)
# compute_mesh_quality_stats(mesh)


# @testset "UNV vs MED" begin
#     ok, report = compare_meshes2(unv, med)
#     @test ok |> (@info "Comparison report" failures = report)
# end






# mesh_dev = mesh
# # mesh_dev = adapt(CUDABackend(), mesh) # uncomment to run on GPU

# # Inlet conditions
# velocity = [5, 0.0, 0.0]
# noSlip = [0.0, 0.0, 0.0]
# nu = 1e-3
# Re = (0.2*velocity[1])/nu
# display(Re)

# model = Physics(
#     time = Steady(),
#     fluid = Fluid{Incompressible}(nu = nu),
#     turbulence = RANS{Laminar}(),
#     energy = Energy{Isothermal}(), 
#     domain = mesh_dev
#     )

# @assign! model momentum U ( 
#     Dirichlet(:inlet, velocity),
#     Neumann(:outlet, 0.0),
#     Wall(:cylinder, noSlip),
#     Neumann(:bottom, 0.0),
#     Neumann(:top, 0.0)
#     # Neumann(:front, 0.0),
#     # Neumann(:back, 0.0)
# )

# @assign! model momentum p (
#     Neumann(:inlet, 0.0),
#     Dirichlet(:outlet, 0.0),
#     Neumann(:cylinder, 0.0),
#     Neumann(:bottom, 0.0),
#     Neumann(:top, 0.0)
#     # Neumann(:front, 0.0),
#     # Neumann(:back, 0.0)
    
# )

# solvers = (
#     U = set_solver(
#         model.momentum.U;
#         solver      = BicgstabSolver, # BicgstabSolver, GmresSolver
#         preconditioner = Jacobi(),
#         convergence = 1e-7,
#         relax       = 1.0,
#         rtol = 1e-4,
#         atol = 1e-5
#     ),
#     p = set_solver(
#         model.momentum.p;
#         solver      = CgSolver, # BicgstabSolver, GmresSolver
#         preconditioner = Jacobi(), #NormDiagonal(),
#         convergence = 1e-7,
#         relax       = 0.8,
#         rtol = 1e-4,
#         atol = 1e-5
#     )
# )

# schemes = (
#     U = set_schemes(divergence=Upwind, gradient=Midpoint),
#     p = set_schemes(gradient=Midpoint)
    
# )


# runtime = set_runtime(iterations=500, write_interval=20, time_step=1) 


# hardware = set_hardware(backend=CPU(), workgroup=1024)

# config = Configuration(
#     solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware)

# GC.gc(true)

# initialise!(model.momentum.U, velocity)
# initialise!(model.momentum.p, 0.0)

# residuals = run!(model, config)