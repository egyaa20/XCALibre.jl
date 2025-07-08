export coupled!



# Note: THIS IS A VERY ROUGH DRAFT!


# fluid_face_lookup = Dict{Vector{Float64}, Int}()

# for pointer in fluid_mesh_faces_pointer
#     fluid_id = fluid_mesh.cell_faces[pointer]
#     coords = fluid_mesh.faces[fluid_id].centre
    
#     fluid_face_lookup[coords] = fluid_id
# end

# face_mapping = Vector{Vector{Int}}()

# for pointer in solid_mesh_faces_pointer
#     solid_id = solid_mesh.cell_faces[pointer]
#     solid_coords = solid_mesh.faces[solid_id].centre
    
#     if haskey(fluid_face_lookup, solid_coords)
#         matching_fluid_id = fluid_face_lookup[solid_coords]
        
#         push!(face_mapping, [solid_id, matching_fluid_id])
#     else
#         println("Warning: No matching fluid face found for solid face ID: $solid_id")
#     end
# end

# println(face_mapping)





function coupled!(mp, configs; output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0) 

    residuals = setup_coupled_solver(
        COUPLED, mp, configs; 
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops
        )
    return residuals
end

# Setup for all compressible algorithms
function setup_coupled_solver(
    solver_variant, mp, configs; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    ) 

    # println(configs[2]) #fluid
    # println(configs[1]) #solid
    # println(mp.couplings.Coupling1.physics1) #model for solid

    solid_mesh = mp.couplings.Coupling1.physics1.domain
    fluid_mesh = mp.couplings.Coupling1.physics2.domain

    solid_model = mp.couplings.Coupling1.physics1
    fluid_model = mp.couplings.Coupling1.physics2

    solid_mesh_faces_pointer = mp.couplings.Coupling1.interface1.IDs_range
    fluid_mesh_faces_pointer = mp.couplings.Coupling1.interface2.IDs_range

    println(solid_mesh_faces_pointer)
    println(fluid_mesh_faces_pointer)

    face_ID_example = solid_mesh.cell_faces[solid_mesh_faces_pointer[1]]
    face_centre_coords_example = solid_mesh.faces[face_ID_example].centre
    println(face_centre_coords_example)



    solid_config = configs[1]
    fluid_config = configs[2]

    @info "WE ARE HERE!!!"

    num_couplings = length(mp.couplings) #for solvers loop

    solver_variants = [PISO, LAPLACET] #please replace this line ASAP..... (solver_variant list)
    #the one stepped here needs to be the last one
    

    # sequence = 1 #1 = forward, 2 = backward

    # Initialise LAPLACE solver
    T_eqn, outputWriter, R_T, time = setup_transient_laplace_solver(
        LAPLACET, solid_model, solid_config;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        coupling=true
        )

    


    # LAPLACE.step()
    LAPLACET(
        solid_model, solid_model.energy, T_eqn, solid_config; 
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        outputWriter, R_T, time, isCoupled=true)






    # DONE > 0) Call solvers' initialisation
    # 1) Map Interface
    # 2) Define Order (e.g. region 1-2 then region 2-3 etc)
    # DONE >  3) Set Sequence = Forward
    # DONE >  4) Do SOLVER1.step()


    # Basically that's it....






    # (; solvers, schemes, runtime, hardware, boundaries) = config

    # @info "THIS IS COUPLED ALGO!!!!"
    # @info "Extracting configuration and input fields..."

    # # model = adapt(hardware.backend, model_in)
    # (; U, p, Uf, pf) = model.momentum
    # (; rho) = model.fluid
    # mesh = model.domain

    # @info "Pre-allocating fields..."
    
    # ∇p = Grad{schemes.p.gradient}(p)
    # mdotf = FaceScalarField(mesh)
    # rhorDf = FaceScalarField(mesh)
    # initialise!(rhorDf, 1.0)
    # mueff = FaceScalarField(mesh)
    # mueffgradUt = VectorField(mesh)
    # divHv = ScalarField(mesh)

    # @info "Defining models..."

    # U_eqn = (
    #     Time{schemes.U.time}(rho, U)
    #     + Divergence{schemes.U.divergence}(mdotf, U) 
    #     - Laplacian{schemes.U.laplacian}(mueff, U) 
    #     == 
    #     - Source(∇p.result)
    #     + Source(mueffgradUt)
    # ) → VectorEquation(U, boundaries.U)

    # if typeof(model.fluid) <: WeaklyCompressible

    #     p_eqn = (
    #         - Laplacian{schemes.p.laplacian}(rhorDf, p) == - Source(divHv)
    #     ) → ScalarEquation(p, boundaries.p)

    # elseif typeof(model.fluid) <: Compressible

    #     pconv = FaceScalarField(mesh)
    #     p_eqn = (
    #         Laplacian{schemes.p.laplacian}(rhorDf, p) 
    #         - Divergence{schemes.p.divergence}(pconv, p) == Source(divHv)
    #     ) → ScalarEquation(p, boundaries.p)

    # end

    # @info "Initialising preconditioners..."

    # @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    # @reset p_eqn.preconditioner = set_preconditioner(solvers.p.preconditioner, p_eqn)

    # @info "Pre-allocating solvers..."
     
    # @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    # @reset p_eqn.solver = _workspace(solvers.p.solver, _b(p_eqn))
  
    # @info "Initialising energy model..."
    # energyModel = initialise(model.energy, model, mdotf, rho, p_eqn, config)

    # @info "Initialising turbulence model..."
    # turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)

    residuals  = solver_variant(
        mp, configs;
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        time, num_couplings, solver_variants)

    return residuals    
end # end function

function COUPLED(
    mp, configs; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0,
    time, num_couplings, solver_variants
    )


    # (; iterations, write_interval, dt) = runtime

    iterations = 10
    dt = 0.1

    sequence = 1 #1 = forward, 2 = backward
    solver_index = 1
    current_solver = solver_variants[solver_index]

    # println(solver_variants)
    # println(length(solver_variants))

    @info "Starting COUPLED loops..."

    for iteration ∈ 1:iterations
        time = iteration * dt

        for coupling ∈ 1:num_couplings
        # for coupling ∈ 1:2
            # Solver selection logic; NOT ADAPTED FOR KERNEL, or is it????
            if sequence == 1
                solver_index += 1
            else
                solver_index -= 1
            end

            if solver_index > length(solver_variants)
                solver_index = 1
            elseif solver_index == 0
                solver_index = length(solver_variants)
            end


            interface_update()
            # solver_step()
        end
        sequence = 3 - sequence #returns 1 if it is 2 now; if it is 2 then returns 1
    end


    return nothing


    
    # return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p, e=R_e)

    # ITERATIONS LOOP
        # SOLVERS LOOP (per number of couplings e.g. 2 solvers coupled = 1 solver loop)
            # Update_Interface() e.g. T1=T2

            # SolverX.step() e.g. Solver2
        # Reverse Sequence
        # Update residuals?


    # That's Basically it...






    
    # # Extract model variables and configuration
    # (; U, p, Uf, pf) = model.momentum
    # (; nu, nuf, rho, rhof) = model.fluid

    # mesh = model.domain
    # p_model = p_eqn.model
    # (; solvers, schemes, runtime, hardware, boundaries) = config
    # (; iterations, write_interval) = runtime
    # (; backend) = hardware
    
    # # rho = get_flux(U_eqn, 1)
    # mdotf = get_flux(U_eqn, 2)
    # mueff = get_flux(U_eqn, 3)
    # mueffgradUt = get_source(U_eqn, 2)
    # rhorDf = get_flux(p_eqn, 1)
    # divHv = get_source(p_eqn, 1)
    # if typeof(model.fluid) <: Compressible
    #     pconv = get_flux(p_eqn, 2)
    # end

    # outputWriter = initialise_writer(output, model.domain)
    
    # @info "Allocating working memory..."

    # # Define aux fields 
    # gradU = Grad{schemes.U.gradient}(U)
    # gradUT = T(gradU)
    # # Uf = FaceVectorField(mesh)
    # S = StrainRate(gradU, gradUT, U, Uf)

    # n_cells = length(mesh.cells)
    # # pf = FaceScalarField(mesh)
    # nueff = FaceScalarField(mesh)
    # prevpf = FaceScalarField(mesh)
    # gradpf = FaceVectorField(mesh)
    # Hv = VectorField(mesh)
    # rD = ScalarField(mesh)
    # Psi = ScalarField(mesh)
    # Psif = FaceScalarField(mesh)

    # mugradUTx = FaceScalarField(mesh)
    # mugradUTy = FaceScalarField(mesh)
    # mugradUTz = FaceScalarField(mesh)

    # divmugradUTx = ScalarField(mesh)
    # divmugradUTy = ScalarField(mesh)
    # divmugradUTz = ScalarField(mesh)

    # # Pre-allocate auxiliary variables
    # TF = _get_float(mesh)
    # # prev = zeros(TF, n_cells)
    # # prev = _convert_array!(prev, backend) 
    # prev = KernelAbstractions.zeros(backend, TF, n_cells) 

    # # Pre-allocate vectors to hold residuals 
    # R_ux = ones(TF, iterations)
    # R_uy = ones(TF, iterations)
    # R_uz = ones(TF, iterations)
    # R_p = ones(TF, iterations)
    # R_e = ones(TF, iterations)
    
    # # Initial calculations
    # time = zero(TF) # assuming time=0
    # interpolate!(Uf, U, config)   
    # correct_boundaries!(Uf, U, boundaries.U, time, config) 
    # grad!(∇p, pf, p, boundaries.p, time, config)
    # thermo_Psi!(model, Psi); thermo_Psi!(model, Psif, config);
    # @. rho.values = Psi.values * p.values
    # @. rhof.values = Psif.values * pf.values
    # flux!(mdotf, Uf, rhof, config)

    # update_nueff!(nueff, nu, model.turbulence, config)
    # @. mueff.values = nueff.values * rhof.values

    # @info "Starting CSIMPLE loops..."

    # progress = Progress(iterations; dt=1.0, showspeed=true)

    # xdir, ydir, zdir = XDir(), YDir(), ZDir()

    # for iteration ∈ 1:iterations
    #     time = iteration

    #     ## CHECK GRADU AND EXPLICIT STRESSES
    #     # grad!(gradU, Uf, U, boundaries.U, time, config) # calculated in `turbulence!``

    #     explicit_shear_stress!(mugradUTx, mugradUTy, mugradUTz, mueff, gradU, config)
    #     div!(divmugradUTx, mugradUTx, config)
    #     div!(divmugradUTy, mugradUTy, config)
    #     div!(divmugradUTz, mugradUTz, config)
        
    #     @. mueffgradUt.x.values = divmugradUTx.values
    #     @. mueffgradUt.y.values = divmugradUTy.values
    #     @. mueffgradUt.z.values = divmugradUTz.values

    #     # Set up and solve momentum equations
        
    #     rx, ry, rz = solve_equation!(U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config)
    #     energy!(energyModel, model, prev, mdotf, rho, mueff, time, config)
    #     thermo_Psi!(model, Psi); thermo_Psi!(model, Psif, config);

    #     # Pressure correction
    #     inverse_diagonal!(rD, U_eqn, config)
    #     interpolate!(rhorDf, rD, config)
    #     @. rhorDf.values *= rhof.values

    #     remove_pressure_source!(U_eqn, ∇p, config)
    #     H!(Hv, U, U_eqn, config)
        
    #     # Interpolate faces
    #     interpolate!(Uf, Hv, config) # Careful: reusing Uf for interpolation
    #     correct_boundaries!(Uf, Hv, boundaries.U, time, config)

    #     if typeof(model.fluid) <: Compressible
    #         flux!(pconv, Uf, config)
    #         @. pconv.values *= Psif.values
    #         flux!(mdotf, Uf, config)
    #         @. mdotf.values *= rhof.values
    #         interpolate!(pf, p, config)
    #         correct_boundaries!(pf, p, boundaries.p, time, config)
    #         @. mdotf.values -= mdotf.values*Psif.values*pf.values/rhof.values
    #         div!(divHv, mdotf, config)

    #     elseif typeof(model.fluid) <: WeaklyCompressible
    #         flux!(mdotf, Uf, config)
    #         @. mdotf.values *= rhof.values
    #         div!(divHv, mdotf, config)
    #     end

    #     # Pressure calculations
    #     rp = 0.0
    #     @. prev = p.values
    #     @. prevpf.values = pf.values
    #     if typeof(model.fluid) <: Compressible
    #         # Ensure diagonal dominance for hyperbolic equations
    #         rp = solve_equation!(p_eqn, p, boundaries.p, solvers.p, config; ref=nothing, irelax=solvers.U.relax)
    #     elseif typeof(model.fluid) <: WeaklyCompressible
    #         rp = solve_equation!(p_eqn, p, boundaries.p, solvers.p, config; ref=nothing)
    #     end

    #     if !isnothing(solvers.p.limit)
    #         pmin = solvers.p.limit[1]; pmax = solvers.p.limit[2]
    #         clamp!(p.values, pmin, pmax)
    #     end

    #     explicit_relaxation!(p, prev, solvers.p.relax, config)

    #     grad!(∇p, pf, p, boundaries.p, time, config) 
    #     limit_gradient!(schemes.p.limiter, ∇p, p, config)

    #     # non-orthogonal correction
    #     for i ∈ 1:ncorrectors
    #         discretise!(p_eqn, p, config)       
    #         apply_boundary_conditions!(p_eqn, boundaries.p, nothing, time, config)
    #         setReference!(p_eqn, pref, 1, config)
    #         nonorthogonal_face_correction(p_eqn, ∇p, rhorDf, config)
    #         update_preconditioner!(p_eqn.preconditioner, p.mesh, config)
    #         rp = solve_system!(p_eqn, solvers.p, p, nothing, config)
    #         explicit_relaxation!(p, prev, solvers.p.relax, config)
            
    #         grad!(∇p, pf, p, boundaries.p, time, config) 
    #         limit_gradient!(schemes.p.limiter, ∇p, p, config)
    #     end

    #     if typeof(model.fluid) <: Compressible
    #         rhorelax = 1.0 #0.01
    #         @. rho.values = rho.values * (1-rhorelax) + Psi.values * p.values * rhorelax
    #         @. rhof.values = rhof.values * (1-rhorelax) + Psif.values * pf.values * rhorelax
    #     else
    #         @. rho.values = Psi.values * p.values
    #         @. rhof.values = Psif.values * pf.values
    #     end

    #     # Velocity and boundaries correction
    #     # correct_face_interpolation!(pf, p, Uf) # not needed added upwind interpolation
    #     # correct_boundaries!(pf, p, boundaries.p, time, config)
    #     # pgrad = face_normal_gradient(p, pf)

    #     if typeof(model.fluid) <: Compressible
    #         # @. mdotf.values += (pconv.values*(pf.values) - pgrad.values*rhorDf.values)  
    #         correct_mass_flux(mdotf, p, rhorDf, config)
    #         @. mdotf.values += pconv.values*(pf.values)
    #     elseif typeof(model.fluid) <: WeaklyCompressible
    #         # @. mdotf.values -= pgrad.values*rhorDf.values
    #         correct_mass_flux(mdotf, p, rhorDf, config)
    #     end

    #     correct_velocity!(U, Hv, ∇p, rD, config)
    #     # interpolate!(Uf, U, config)
    #     # correct_boundaries!(Uf, U, boundaries.U, time, config)
        
    #     turbulence!(turbulenceModel, model, S, prev, time, config) 
    #     update_nueff!(nueff, nu, model.turbulence, config)

    #     @. mueff.values = rhof.values*nueff.values

    #     R_ux[iteration] = rx
    #     R_uy[iteration] = ry
    #     R_uz[iteration] = rz
    #     R_p[iteration] = rp

    #     Uz_convergence = true
    #     if typeof(mesh) <: Mesh3
    #         Uz_convergence = rz <= solvers.U.convergence
    #     end

    #     if (R_ux[iteration] <= solvers.U.convergence && 
    #         R_uy[iteration] <= solvers.U.convergence && 
    #         Uz_convergence &&
    #         R_p[iteration] <= solvers.p.convergence &&
    #         turbulenceModel.state.converged)

    #         progress.n = iteration
    #         finish!(progress)
    #         @info "Simulation converged in $iteration iterations!"
    #         if !signbit(write_interval)
    #             save_output(model, outputWriter, time, config)
    #         end
    #         break
    #     end

    #     ProgressMeter.next!(
    #         progress, showvalues = [
    #             (:iter,iteration),
    #             (:Ux, R_ux[iteration]),
    #             (:Uy, R_uy[iteration]),
    #             (:Uz, R_uz[iteration]),
    #             (:p, R_p[iteration]),
    #             turbulenceModel.state.residuals...,
    #             energyModel.state.residuals
    #             ]
    #         )

    #     if iteration%write_interval + signbit(write_interval) == 0      
    #         save_output(model, outputWriter, time, config)
    #     end

    # end # end for loop

    # return (Ux=R_ux, Uy=R_uy, Uz=R_uz, p=R_p, e=R_e)
end







function solver_step()
    
end



function interface_update(Field1, Field2, face_mapping)
    for face ∈ 1:face_mapping

    end
end
















# function explicit_shear_stress!(mugradUTx::FaceScalarField, mugradUTy::FaceScalarField, mugradUTz::FaceScalarField, mueff, gradU, config)
#     (; hardware) = config
#     (; backend, workgroup) = hardware

#     (; faces, boundary_cellsID) = mugradUTx.mesh

#     n_faces = length(faces)
#     n_bfaces = length(boundary_cellsID)
#     n_ifaces = n_faces - n_bfaces

#     ndrange = n_ifaces
#     kernel! = _explicit_shear_stress_internal!(_setup(backend, workgroup, ndrange)...)
#     kernel!(mugradUTx, mugradUTy, mugradUTz, mueff, gradU, faces, n_bfaces)
#     # KernelAbstractions.synchronize(backend)

#     ndrange=n_bfaces
#     kernel! = _explicit_shear_stress_boundaries!(_setup(backend, workgroup, ndrange)...)
#     kernel!(mugradUTx, mugradUTy, mugradUTz, mueff, gradU, faces)
#     # KernelAbstractions.synchronize(backend)
# end

# @kernel function _explicit_shear_stress_internal!(
#     mugradUTx, mugradUTy, mugradUTz, mueff, gradU, faces,n_bfaces)
#     i = @index(Global)

#     fID = i + n_bfaces
#     face = faces[fID]
#     (; area, normal, ownerCells) = face 
#     cID1 = ownerCells[1]
#     cID2 = ownerCells[2]
#     gradUf = 0.5*(gradU[cID1] + gradU[cID2]) # should this be the transpose of gradU?
#     gradUf_projection = gradUf*normal
#     trace = 2/3*sum(diag(gradUf))
#     mueffi = mueff[fID]
#     mugradUTx[fID] = mueffi*(gradUf_projection[1] - trace)*area
#     mugradUTy[fID] = mueffi*(gradUf_projection[2] - trace)*area
#     mugradUTz[fID] = mueffi*(gradUf_projection[3] - trace)*area
# end

# @kernel function _explicit_shear_stress_boundaries!(
#     mugradUTx, mugradUTy, mugradUTz, mueff, gradU, faces)
#     fID = @index(Global)

#     face = faces[fID]
#     (; area, normal, ownerCells) = face 
#     cID1 = ownerCells[1]
#     gradUi = gradU[cID1]
#     trace = 2/3*sum(diag(gradUi))
#     gradUi_projection = gradUi*normal
#     mueffi = mueff[fID]
#     mugradUTx[fID] = mueffi*(gradUi_projection[1] - trace)*area
#     mugradUTy[fID] = mueffi*(gradUi_projection[2] - trace)*area
#     mugradUTz[fID] = mueffi*(gradUi_projection[3] - trace)*area
# end
