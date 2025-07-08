export simple!

mutable struct SolverAux{E}
    mdotf::Any
    nueff::Any
    rDf::Any
    divHv::Any
    S::Any
    Hv::Any
    rD::Any
    prev::Any
    R_ux::Any
    R_uy::Any
    R_uz::Any
    R_p::Any
    time::Any
    extras::E # A container for any number of extra, optional variables
end

"""
    simple!(model_in, config; 
        output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0)

Incompressible variant of the SIMPLE algorithm to solving coupled momentum and mass conservation equations.

# Input arguments

- `model` reference to a `Physics` model defined by the user.
- `config` Configuration structure defined by the user with solvers, schemes, runtime and hardware structures configuration details.
- `output` select the format used for simulation results from `VTK()` or `OpenFOAM` (default = `VTK()`)
- `pref` Reference pressure value for cases that do not have a pressure defining BC. Incompressible solvers only (default = `nothing`)
- `ncorrectors` number of non-orthogonality correction loops (default = `0`)
- `inner_loops` number to inner loops used in transient solver based on PISO algorithm (default = `0`)

# Output

This function returns a `NamedTuple` for accessing the residuals (e.g. `residuals.Ux`) with the following entries:

- `Ux` Vector of x-velocity residuals for each iteration.
- `Uy` Vector of y-velocity residuals for each iteration.
- `Uz` Vector of y-velocity residuals for each iteration.
- `p` Vector of pressure residuals for each iteration.

"""
function simple!(
    model, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0
    )

    residuals = setup_incompressible_solvers(
        SIMPLE, model, config; 
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        coupling=false
        )

    return residuals
end

# Setup for all incompressible algorithms
function setup_incompressible_solvers(
    solver_variant, model, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0, coupling=false
    ) 


    (; solvers, schemes, runtime, hardware, boundaries) = config


    (; nu) = model.medium
    (; iterations, write_interval) = runtime
    (; backend) = hardware




    @info "Extracting configuration and input fields..."

    (; U, p, Uf, pf) = model.momentum
    mesh = model.domain

    @info "Pre-allocating fields..."
    
    ∇p = Grad{schemes.p.gradient}(p)
    mdotf = FaceScalarField(mesh)
    rDf = FaceScalarField(mesh)
    initialise!(rDf, 1.0)
    nueff = FaceScalarField(mesh)
    divHv = ScalarField(mesh)

    @info "Defining models..."

    U_eqn = (
        Time{schemes.U.time}(U)
        + Divergence{schemes.U.divergence}(mdotf, U) 
        - Laplacian{schemes.U.laplacian}(nueff, U) 
        == 
        - Source(∇p.result)
    ) → VectorEquation(U, boundaries.U)

    p_eqn = (
        - Laplacian{schemes.p.laplacian}(rDf, p) == - Source(divHv)
    ) → ScalarEquation(p, boundaries.p)

    @info "Initialising preconditioners..."

    # @reset U_eqn.preconditioner = set_preconditioner(
    #                 solvers.U.preconditioner, U_eqn, boundaries.U, config)
    # @reset p_eqn.preconditioner = set_preconditioner(
    #                 solvers.p.preconditioner, p_eqn, boundaries.p, config)

    @reset U_eqn.preconditioner = set_preconditioner(solvers.U.preconditioner, U_eqn)
    @reset p_eqn.preconditioner = set_preconditioner(solvers.p.preconditioner, p_eqn)

    @info "Pre-allocating solvers..."

    @reset U_eqn.solver = _workspace(solvers.U.solver, _b(U_eqn, XDir()))
    @reset p_eqn.solver = _workspace(solvers.p.solver, _b(p_eqn))

    @info "Initialising turbulence model..."
    turbulenceModel = initialise(model.turbulence, model, mdotf, p_eqn, config)


    # NEW INITIALISATION LOGIC



    
    mdotf = get_flux(U_eqn, 2)
    nueff = get_flux(U_eqn, 3)
    rDf = get_flux(p_eqn, 1)
    divHv = get_source(p_eqn, 1)

    outputWriter = initialise_writer(output, model.domain)
    
    @info "Allocating working memory..."

    # Define aux fields 
    gradU = Grad{schemes.U.gradient}(U)
    gradUT = T(gradU)
    S = StrainRate(gradU, gradUT, U, Uf)

    n_cells = length(mesh.cells)
    Hv = VectorField(mesh)
    rD = ScalarField(mesh)

    # Pre-allocate auxiliary variables
    TF = _get_float(mesh)
    # prev = zeros(TF, n_cells)
    # prev = _convert_array!(prev, backend) 
    prev = KernelAbstractions.zeros(backend, TF, n_cells) 

    # Pre-allocate vectors to hold residuals 
    R_ux = ones(TF, iterations)
    R_uy = ones(TF, iterations)
    R_uz = ones(TF, iterations)
    R_p = ones(TF, iterations)
    
    # Initial calculations
    time = zero(TF) # assuming time=0
    interpolate!(Uf, U, config)   
    correct_boundaries!(Uf, U, boundaries.U, time, config)
    flux!(mdotf, Uf, config)
    grad!(∇p, pf, p, boundaries.p, time, config)
    limit_gradient!(schemes.p.limiter, ∇p, p, config)


    update_nueff!(nueff, nu, model.turbulence, config)







    xdir, ydir, zdir = XDir(), YDir(), ZDir()

    extra_vars = NamedTuple()

    

    # IF PISO
    if typeof(model.time) <: Transient
        Uf = FaceVectorField(mesh)
        pf = FaceScalarField(mesh)
        cellsCourant =adapt(backend, zeros(TF, length(mesh.cells)))

        extra_vars = (
            cellsCourant = cellsCourant,
        )
    end



    # Instantiate the SolverAux struct with all the required variables.
    aux = SolverAux(
        mdotf,
        nueff,
        rDf,
        divHv,
        S,
        Hv,
        rD,
        prev,
        R_ux,
        R_uy,
        R_uz,
        R_p,
        time,
        extra_vars
    )


    # .


    residuals  = solver_variant(
        model, turbulenceModel, ∇p, U_eqn, p_eqn, config; 
        output=output,
        pref=pref, 
        ncorrectors=ncorrectors, 
        inner_loops=inner_loops,
        aux=aux,
        xdir, ydir, zdir, outputWriter, isCoupled=false)

    return residuals
end # end function

function SIMPLE(
    model, turbulenceModel, ∇p, U_eqn, p_eqn, config; 
    output=VTK(), pref=nothing, ncorrectors=0, inner_loops=0,
    aux,
    xdir, ydir, zdir, outputWriter, isCoupled
    )

    # Extract model variables and configuration
    (; U, p, Uf, pf) = model.momentum # cell center and face values
    (; nu) = model.medium
    mesh = model.domain
    (; solvers, schemes, runtime, hardware, boundaries) = config
    (; iterations, write_interval) = runtime
    (; backend) = hardware
    
    # mdotf = get_flux(U_eqn, 2)
    # nueff = get_flux(U_eqn, 3)
    # rDf = get_flux(p_eqn, 1)
    # divHv = get_source(p_eqn, 1)

    # outputWriter = initialise_writer(output, model.domain)
    
    # @info "Allocating working memory..."

    # # Define aux fields 
    # gradU = Grad{schemes.U.gradient}(U)
    # gradUT = T(gradU)
    # S = StrainRate(gradU, gradUT, U, Uf)

    # n_cells = length(mesh.cells)
    # Hv = VectorField(mesh)
    # rD = ScalarField(mesh)

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
    
    # # Initial calculations
    # time = zero(TF) # assuming time=0
    # interpolate!(Uf, U, config)   
    # correct_boundaries!(Uf, U, boundaries.U, time, config)
    # flux!(mdotf, Uf, config)
    # grad!(∇p, pf, p, boundaries.p, time, config)
    # limit_gradient!(schemes.p.limiter, ∇p, p, config)


    # update_nueff!(nueff, nu, model.turbulence, config)


    # @info "Starting LAPLACE loops..."

    # if (isCoupled)
    #     println("COUPLED!!!!")

    #     rt = solve_equation!(T_eqn, T, boundaries.T, solvers, config; time=time)

    #     energy!(model.energy, model, T, rDf, rhocp, k, kf, cp, rho, material, config) # does nothing for diffusion energy model

    #     # R_T[iteration] = rt 
        
    # else
    #     progress = Progress(iterations; dt=1.0, showspeed=true)



    @info "Starting SIMPLE loops..."
    
    if (isCoupled)
        println("COUPLED!!!!")
    else
        progress = Progress(iterations; dt=1.0, showspeed=true)

        # xdir, ydir, zdir = XDir(), YDir(), ZDir()

        for iteration ∈ 1:iterations
            aux.time = iteration

            rx, ry, rz = solve_equation!(U_eqn, U, boundaries.U, solvers.U, xdir, ydir, zdir, config)
            
            # Pressure correction
            inverse_diagonal!(aux.rD, U_eqn, config)
            interpolate!(aux.rDf, aux.rD, config)
            # correct_boundaries!(rDf, rD, rD.BCs, time, config) # ADDED FOR PERIODIC BCS
            remove_pressure_source!(U_eqn, ∇p, config)
            H!(aux.Hv, U, U_eqn, config)
            
            # Interpolate faces
            interpolate!(Uf, aux.Hv, config) # Careful: reusing Uf for interpolation
            correct_boundaries!(Uf, aux.Hv, boundaries.U, aux.time, config)

            # old approach
            # div!(divHv, Uf, config) 

            # new approach
            flux!(aux.mdotf, Uf, config)
            div!(aux.divHv, aux.mdotf, config)
            
            # Pressure calculations
            @. aux.prev = p.values
            rp = solve_equation!(p_eqn, p, boundaries.p, solvers.p, config; ref=pref)
            explicit_relaxation!(p, aux.prev, solvers.p.relax, config)
            
            grad!(∇p, pf, p, boundaries.p, aux.time, config)  #Leave time for test
            limit_gradient!(schemes.p.limiter, ∇p, p, config)

            # non-orthogonal correction
            for i ∈ 1:ncorrectors
                # @. prev = p.values
                discretise!(p_eqn, p, config)       
                apply_boundary_conditions!(p_eqn, boundaries.p, nothing, aux.time, config)
                # setReference!(p_eqn, pref, 1, config)
                nonorthogonal_face_correction(p_eqn, ∇p, aux.rDf, config)
                # update_preconditioner!(p_eqn.preconditioner, p.mesh, config)
                rp = solve_system!(p_eqn, solvers.p, p, nothing, config)
                explicit_relaxation!(p, aux.prev, solvers.p.relax, config)
                grad!(∇p, pf, p, boundaries.p, aux.time, config) 
                limit_gradient!(schemes.p.limiter, ∇p, p, config)
            end

            # explicit_relaxation!(p, prev, solvers.p.relax, config)

            # Velocity and boundaries correction

            # old approach
            # correct_velocity!(U, Hv, ∇p, rD, config)
            # interpolate!(Uf, U, config)
            # correct_boundaries!(Uf, U, boundaries.U, time, config)
            # flux!(mdotf, Uf, config) 

            # new approach

            # 1. using velocity from momentum equation
            # interpolate!(Uf, U, config)
            # correct_boundaries!(Uf, U, boundaries.U, time, config)
            # flux!(mdotf, Uf, config)
            correct_mass_flux(aux.mdotf, p, aux.rDf, config)
            correct_velocity!(U, aux.Hv, ∇p, aux.rD, config)

            turbulence!(turbulenceModel, model, aux.S, aux.prev, aux.time, config) 
            update_nueff!(aux.nueff, nu, model.turbulence, config)

            aux.R_ux[iteration] = rx
            aux.R_uy[iteration] = ry
            aux.R_uz[iteration] = rz
            aux.R_p[iteration] = rp

            Uz_convergence = true
            if typeof(mesh) <: Mesh3
                Uz_convergence = rz <= solvers.U.convergence
            end

            if (aux.R_ux[iteration] <= solvers.U.convergence && 
                aux.R_uy[iteration] <= solvers.U.convergence && 
                Uz_convergence &&
                aux.R_p[iteration] <= solvers.p.convergence &&
                turbulenceModel.state.converged)

                progress.n = iteration
                finish!(progress)
                @info "Simulation converged in $iteration iterations!"
                if !signbit(write_interval)
                    save_output(model, outputWriter, aux.time, config)
                end
                break
            end

            ProgressMeter.next!(
                progress, showvalues = [
                    (:iter,iteration),
                    (:Ux, aux.R_ux[iteration]),
                    (:Uy, aux.R_uy[iteration]),
                    (:Uz, aux.R_uz[iteration]),
                    (:p, aux.R_p[iteration]),
                    turbulenceModel.state.residuals...
                    ]
                )

            if iteration%write_interval + signbit(write_interval) == 0      
                save_output(model, outputWriter, aux.time, config)
            end

        end # end for loop
    end
    
    return (Ux=aux.R_ux, Uy=aux.R_uy, Uz=aux.R_uz, p=aux.R_p)
end

### TEMP LOCATION FOR PROTOTYPING - NONORTHOGONAL CORRECTION 

function nonorthogonal_face_correction(eqn, grad, flux, config)
    mesh = grad.mesh
    (; faces, cells, boundary_cellsID) = mesh

    (; hardware) = config
    (; backend, workgroup) = hardware

    (; b) = eqn.equation
    
    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_ifaces
    kernel! = _nonorthogonal_face_correction(_setup(backend, workgroup, ndrange)...)
    kernel!(b, grad, flux, faces, cells, n_bfaces)
    # KernelAbstractions.synchronize(backend)
end

@kernel function _nonorthogonal_face_correction(b, grad, flux, faces, cells, n_bfaces)
    i = @index(Global)
    fID = i + n_bfaces
    face = faces[fID]
    (; ownerCells, area, normal, e, delta) = face
    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    cell1 = cells[cID1]
    cell2 = cells[cID2]

    xf = face.centre
    xC = cell1.centre
    xN = cell2.centre
    
    # Calculate weights using normal functions
    # weight = norm(xf - xC)/norm(xN - xC)
    # weight = norm(xf - xN)/norm(xN - xC)

    dPN = cell2.centre - cell1.centre

    (; values) = grad.field
    weight, df = correction_weight(cells, faces, fID)
    # weight = face.weight
    gradi = weight*grad[cID1] + (1.0 - weight)*grad[cID2]
    gradf = gradi + ((values[cID2] - values[cID1])/delta - (gradi⋅e))*e
    # gradf = gradi

    Sf = area*normal
    # Ef = ((Sf⋅Sf)/(Sf⋅e))*e # original
    Ef = dPN*(norm(normal)^2/(dPN⋅normal))*area
    T_hat = Sf - Ef # original
    faceCorrection = flux[fID]*gradf⋅T_hat

    Atomix.@atomic b[cID1] += faceCorrection #*cell1.volume
    Atomix.@atomic b[cID2] -= faceCorrection #*cell2.volume # should this be -ve?

    # Atomix.@atomic b[cID1] -= faceCorrection #*cell1.volume
    # Atomix.@atomic b[cID2] += faceCorrection #*cell2.volume # should this be -ve?
        
end

# +- => good match
# -+ => looks worse at edges for gradient
# -- => looks bad on top-right corner for gradient
# ++ => looks bad on left grad and oscillations on the right

function correction_weight(cells, faces, fi)
    (; ownerCells, centre) = faces[fi]
    cID1 = ownerCells[1]
    cID2 = ownerCells[2]
    c1 = cells[cID1].centre
    c2 = cells[cID2].centre
    c1_f = centre - c1
    c1_c2 = c2 - c1
    q = (c1_f⋅c1_c2)/(c1_c2⋅c1_c2)
    f_prime = c1 - q*(c1 - c2)
    w = norm(c2 - f_prime)/norm(c2 - c1)
    df = centre - f_prime
    return w, df
end

### TEMP LOCATION FOR PROTOTYPING

function correct_mass_flux(mdotf, p, rDf, config)
    # sngrad = FaceScalarField(mesh)
    (; faces, cells, boundary_cellsID) = mdotf.mesh
    (; hardware) = config
    (; backend, workgroup) = hardware

    n_faces = length(faces)
    n_bfaces = length(boundary_cellsID)
    n_ifaces = n_faces - n_bfaces

    ndrange = n_ifaces # length(n_ifaces) was a BUG! should be n_ifaces only!!!!
    kernel! = _correct_mass_flux(_setup(backend, workgroup, ndrange)...)
    kernel!(mdotf, p, rDf, faces, cells, n_bfaces)
    # KernelAbstractions.synchronize(backend)
end

@kernel function _correct_mass_flux(mdotf, p, rDf, faces, cells, n_bfaces)
    i = @index(Global)
    fID = i + n_bfaces

    @inbounds begin 
        face = faces[fID]
        (; area, normal, ownerCells, delta) = face 
        cID1 = ownerCells[1]
        cID2 = ownerCells[2]
        p1 = p[cID1]
        p2 = p[cID2]
        face_grad = area*(p2 - p1)/delta # best option so far!
        mdotf[fID] -= face_grad*rDf[fID]
    end
end