
using XCALibre
using LinearAlgebra

function tetra_checker()
    cell = mesh.cells[1]

    faces_pointer = cell.faces_range
    nodes_array = mesh.cell_nodes[faces_pointer]
    
    vertex_ids = nodes_array
    cell_coords = [mesh.nodes[ID].coords for ID in vertex_ids]

    result = is_tetrahedron(cell_coords) # if true then it is tetrahedra element
    return result

    # would not work if it is a mixed mesh (e.g. tetra with prisms)
    # improvement: loop over all cells and sort tetra elements and other types of elements?
end

function volume_check_simple(expected_volume) #sums individual volumes by taking the volume value from mesh data
    sum_of_volumes = 0.0

    for cell ∈ mesh.cells
        sum_of_volumes += cell.volume
    end

    error = 100*(expected_volume-sum_of_volumes)/expected_volume
    return sum_of_volumes, error
end




function volume_check_complex(expected_volume) 
#decompises non-tetrahedra mesh into tetrahedra if required => computes cells' volume using faces area data from mesh => sums it up

# Loop over all cells
    # Extract centroid coordinates
    # Loop over all faces
        # Split quadrilateral face into triangles
        # Connect newly created triangular surface's vertices to the centroid
        # Compute the volume for newly created tetrahedra and save it (using compute_tetrahedron_volume function)
    # Compare the volume of tetrahedras to the volume of polyhedra element


    total_volume = 0.0

    for cell ∈ mesh.cells

        centroid = cell.centre
        cell_volume = 0.0
        faces_pointer = cell.faces_range
        faces_array = mesh.cell_faces[faces_pointer]

        is_mesh_tetrahedral = tetra_checker()

        #if mesh is non tetra, then decompose it into tetra
        #if it is tetra, compute volumes straight away

        if (is_mesh_tetrahedral)
            for face ∈ faces_array
                nodes_pointer = mesh.faces[face].nodes_range
                nodes_array = mesh.face_nodes[nodes_pointer]
    
                
                nodes = [mesh.nodes[ID] for ID in nodes_array]

                tetra = [nodes[1], nodes[2], nodes[3], nodes[4]]
                vol = compute_tetrahedron_volume(tetra)

                total_volume += vol
            end
        else
            for face ∈ faces_array
                nodes_pointer = mesh.faces[face].nodes_range
                nodes_array = mesh.face_nodes[nodes_pointer]
    

                nodes = [mesh.nodes[ID] for ID in nodes_array]
                n = length(nodes)
    
                if n < 3
                    error("Face has less than 3 nodes")
                end
    
                v1 = nodes[1].coords #define vertex
                for i in 1:(n-2)
                    v2 = nodes[i+1].coords
                    v3 = nodes[i+2].coords
                    
                    tetra = [centroid, v1, v2, v3]
                    
                    vol = compute_tetrahedron_volume(tetra)
                    cell_volume += vol
                end
            end
    
            total_volume += cell_volume # add tetra element volume to the total calculation
        end
    end

    error = 100*(expected_volume-total_volume)/expected_volume
    return total_volume, error
end



function nsign_check()

    for (face_ID, face) in enumerate(mesh.faces)

        cell1_ID = face.ownerCells[1]
        cell2_ID = face.ownerCells[2]


        boundary_face = false

        if cell1_ID == cell2_ID
            boundary_face = true
            cell = mesh.cells[cell1_ID]

            # the difference is that for boundary faces, mesh.cell_faces will only return internal faces (so findfirst() is useless)
            # do we even need to check those?
            continue
        else
            cell1_faces_range = mesh.cells[cell1_ID].faces_range
            cell2_faces_range = mesh.cells[cell2_ID].faces_range #pointer
            
            face_IDs_1 = mesh.cell_faces[cell1_faces_range]
            face_IDs_2 = mesh.cell_faces[cell2_faces_range]

            face_index_1 = findfirst(==(face_ID), face_IDs_1) #looks for index in array
            face_index_2 = findfirst(==(face_ID), face_IDs_2) #looks for index in array

            actual_nsign_1 = mesh.cell_nsign[cell1_faces_range][face_index_1]
            actual_nsign_2 = mesh.cell_nsign[cell2_faces_range][face_index_2]

            if cell1_ID > cell2_ID
                if actual_nsign_1 != -1
                    error("Expected nsign to be -1 since cell1_ID ($cell1_ID) > cell2_ID ($cell2_ID), but got $actual_nsign_1")
                end
            elseif cell2_ID > cell1_ID
                if actual_nsign_2 != -1
                    error("Expected nsign to be -1 since cell2_ID ($cell2_ID) > cell1_ID ($cell1_ID), but got $actual_nsign_2")
                end
            end
        end
    end

    return true
end


function e_value_check() 
#takes coordinates of two centroids => computes unit vector between them => compares with e value from the mesh

    for (face_id, face) in enumerate(mesh.faces)
        face_id = 40
        face = mesh.faces[face_id]
        cell1_ID = face.ownerCells[1]
        cell2_ID = face.ownerCells[2]

        boundary_face = false

        if cell1_ID == cell2_ID
            boundary_face = true
            cell = mesh.cells[cell1_ID]
            
            # e is expected to point outwards of the domain e.g. be the same as normal
            # check if dot product equals to one?

            result = dot(face.e, face.normal)

            if result != 1.0
                error("Unit vector e at boundary face is not in the direction of face's normal")
            end
            
        else
            cell1 = mesh.cells[cell1_ID]
            cell2 = mesh.cells[cell2_ID]

            centroid1 = cell1.centre
            centroid2 = cell2.centre

            if cell1_ID < cell2_ID
                direction = centroid2 .- centroid1
            else
                direction = centroid1 .- centroid2
            end

            magnitude = norm(direction)
            unit_vector = direction / magnitude
            if (unit_vector != face.e)
                error("Unit vector e is not properly calculated")
            end
        end
    end

    return true
end


function compute_tetrahedron_volume(cell_vertices) #self-explanitory
    A, B, C, D = cell_vertices

    AB = B .- A
    AC = C .- A
    AD = D .- A

    cross_product = cross(AC, AD)
    tripple_product = dot(AB, cross_product)
    volume = abs(tripple_product) / 6.0

    return volume
end

function is_tetrahedron(cell_vertices) #logic: check if there are 4 vertices and that the volume can be computed (found this method online)
    if length(cell_vertices) != 4
        return false
    else
        volume = compute_tetrahedron_volume(cell_vertices)

        epsilon = 1e-10 #floating-point error threshold
        return volume > epsilon #make sure we avoid floating-point errors
    end
    
end



function general_mesh_check(expected_volume)

    #Stage 1: tetra check

    println("[Checking if mesh is tetrahedral....]")
    is_mesh_tetrahedral = tetra_checker()
    if is_mesh_tetrahedral
        println("Mesh is tetrahedral")
    else
        println("Mesh is not tetrahedral... Decomposition might be required for some checks")
    end

    #Stage 2: simple volume check based on cell volume data recorded in mesh

    println("\n[Running simple volume inspection...]")
    simple_volume_result, simple_volume_error = volume_check_simple(expected_volume)
    println("Simple volume inspection completed.\nExpected volume: $expected_volume\nComputed volume: $simple_volume_result\nError: $simple_volume_error %")
    
    #Stage 3: complex volume check (split mesh into tetra if required and calculate individual volumes based on geometry)
    println("\n[Running complex volume inspection...]")
    complex_volume_result, complex_volume_error = volume_check_complex(expected_volume)
    println("Complex volume inspection completed.\nExpected volume: $expected_volume\nComputed volume: $complex_volume_result\nError: $complex_volume_error %")
   
    #Stage 4: n_signs check
    println("\n[Running n signs inspection...]")
    nsign_inspection_state = nsign_check()
    if nsign_inspection_state
        println("n sign inspection completed. Everything is correct.")
    end

    #Stage 5: e vectors check
    println("\n[Running e vectors inspection...]")
    e_inspection_state = e_value_check()
    if e_inspection_state
        println("e vectors inspection completed. Everything is correct.")
    end
end





grids_dir = pkgdir(XCALibre, "examples/testing_grids")

grid = "mesh_test.unv"

mesh_file = joinpath(grids_dir, grid)

mesh = UNV3D_mesh(mesh_file, scale=0.001)


general_mesh_check(1*0.2*0.2)









mesh_dev = mesh
# mesh_dev = adapt(CUDABackend(), mesh) # uncomment to run on GPU

# Inlet conditions
velocity = [5, 0.0, 0.0]
noSlip = [0.0, 0.0, 0.0]
nu = 1e-3
Re = (0.2*velocity[1])/nu
display(Re)

model = Physics(
    time = Steady(),
    fluid = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(), 
    domain = mesh_dev
    )

@assign! model momentum U ( 
    Dirichlet(:inlet, velocity),
    Neumann(:outlet, 0.0),
    Wall(:cylinder, noSlip),
    Neumann(:bottom, 0.0),
    Neumann(:top, 0.0)
    # Neumann(:front, 0.0),
    # Neumann(:back, 0.0)
)

@assign! model momentum p (
    Neumann(:inlet, 0.0),
    Dirichlet(:outlet, 0.0),
    Neumann(:cylinder, 0.0),
    Neumann(:bottom, 0.0),
    Neumann(:top, 0.0)
    # Neumann(:front, 0.0),
    # Neumann(:back, 0.0)
    
)

solvers = (
    U = set_solver(
        model.momentum.U;
        solver      = BicgstabSolver, # BicgstabSolver, GmresSolver
        preconditioner = Jacobi(),
        convergence = 1e-7,
        relax       = 1.0,
        rtol = 1e-4,
        atol = 1e-5
    ),
    p = set_solver(
        model.momentum.p;
        solver      = CgSolver, # BicgstabSolver, GmresSolver
        preconditioner = Jacobi(), #NormDiagonal(),
        convergence = 1e-7,
        relax       = 0.8,
        rtol = 1e-4,
        atol = 1e-5
    )
)

schemes = (
    U = set_schemes(divergence=Upwind, gradient=Midpoint),
    p = set_schemes(gradient=Midpoint)
    
)


runtime = set_runtime(iterations=500, write_interval=20, time_step=1) 


hardware = set_hardware(backend=CPU(), workgroup=1024)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware)

GC.gc(true)

initialise!(model.momentum.U, velocity)
initialise!(model.momentum.p, 0.0)

residuals = run!(model, config)