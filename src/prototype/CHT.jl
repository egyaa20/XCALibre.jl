# using Revise
# using Adapt
using XCALibre

using CUDA

grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
# grid = "backwardFacingStep_10mm.unv"
grid = "summer_2d_5x10.unv"
# grid = "summer_3d_extruded_pipe.unv"
mesh_file = joinpath(grids_dir, grid)

mesh = UNV2D_mesh(mesh_file, scale=0.001)
# mesh = UNV3D_mesh(mesh_file, scale=0.001)

backend = CPU(); workgroup = 1024; activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)

mesh_dev_1 = adapt(backend, mesh)
mesh_dev_2 = adapt(backend, mesh)






solid_model = Physics(
    time = Transient(),
    medium = Solid{Uniform}(k = 16.2),
    energy = Energy{CryogenicConduction}(material = :Aluminium, rho = 8000.0),
    domain = mesh_dev_1
    )

fluid_model = Physics(
    time = Transient(),
    fluid = Fluid{Incompressible}(nu = nu),
    turbulence = RANS{KOmega}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev_2
    )

BCs = assign(
    region = mesh_dev_1,
    (
        T = [     
            Dirichlet(:inlet, 250),
            Zerogradient(:outlet),    
            Zerogradient(:walls)
            # Dirichlet(:walls, 50)      
        ],
    )
)





# TBD: interfacing BC
# TBD: allow LAPLACET and SIMPLE/PISO for coupling
# TBD: scalable multiphysics solver (calls LAPLACET, PISO, exchanges data over the interface)




interface12 = 1 #Dummy interface 

mp = MultiPhysics(
    Coupling1 = Coupling(solid_model, fluid_model, interface12)
#   Coupling2 = Coupling(fluid_model, solid_model_2, interface23),
)







solvers = (
    T = SolverSetup(
        solver      = Cg(), # Bicgstab(), Gmres()
        preconditioner = Jacobi(), # Jacobi(), #NormDiagonal(), # DILU()
        convergence = 1e-8,
        relax       = 0.8,
        rtol = 1e-4,
        atol = 1e-5
    )
)

schemes = (
    T = Schemes(laplacian = Linear)
)


runtime = Runtime(
    iterations=10, write_interval=1, time_step=0.1) #0.1 * 100 = 10 sec

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)

GC.gc(true)

initialise!(model.energy.T, 100.0)

# residuals = run!(model, config)