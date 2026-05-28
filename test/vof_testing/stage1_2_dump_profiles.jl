# Helper for plot_verification.py — runs the three α-transport configurations
# on quad40.unv (matching stages 1 & 2) and dumps the centre-row α profiles
# to JSON for plotting.
#
# Output: alpha_profiles_quad40.json
#         {
#           "upwind":     {"x": [...], "alpha": [...]},
#           "vanleer":    {"x": [...], "alpha": [...]},
#           "compressed": {"x": [...], "alpha": [...]}
#         }

using XCALibre
using JSON

GRID         = "quad40.unv"
Δx           = 25.0
U_FREESTREAM = 1.0
X_INTERFACE  = 500.0
SAMPLING_Y   = 500.0
t_END        = 100.0
dt           = 1.0
ITERATIONS   = round(Int, t_END / dt)

CONFIGS = [
    ("upwind",     true,  0.0),
    ("vanleer",    false, 0.0),
    ("compressed", false, 1.0),
]

function run_and_extract(pure_upwind::Bool, c_alpha::Float64)
    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh_file = joinpath(grids_dir, GRID)
    mesh = UNV2D_mesh(mesh_file, scale=1.0)

    backend  = CPU(); workgroup = AutoTune(); activate_multithread(backend)
    hardware = Hardware(backend=backend, workgroup=workgroup)
    mesh_dev = adapt(backend, mesh)

    gravity = Gravity([0.0, 0.0, 0.0])
    model = Physics(
        time = Transient(),
        fluid = Fluid{Multiphase}(
            model = VOF(sigma=0.0, cAlpha=c_alpha, cycles=1, smooth=0,
                        pure_upwind=pure_upwind),
            phases = (Phase(rho=1.0, mu=1e-6), Phase(rho=1.0, mu=1e-6)),
            gravity = gravity,
        ),
        turbulence = RANS{Laminar}(),
        energy     = Energy{Isothermal}(),
        domain     = mesh_dev,
    )

    inlet_velocity = [U_FREESTREAM, 0.0, 0.0]
    BCs = assign(region = mesh_dev, (
        U = [
            Dirichlet(:inlet, inlet_velocity), Zerogradient(:outlet),
            Zerogradient(:top), Zerogradient(:bottom),
        ],
        p_rgh = [
            Zerogradient(:inlet), Dirichlet(:outlet, 0.0),
            Zerogradient(:top), Zerogradient(:bottom),
        ],
        alpha = [
            Dirichlet(:inlet, 1.0), Zerogradient(:outlet),
            Zerogradient(:top), Zerogradient(:bottom),
        ],
    ))

    schemes = (
        U     = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
        p     = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
        p_rgh = Schemes(time=Euler, gradient=Gauss,    laplacian=Linear),
        alpha = Schemes(time=Euler, divergence=Upwind, laplacian=Linear),
    )
    solvers = (
        U = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                        convergence=1e-7, relax=1.0, rtol=0.0, atol=1e-7),
        p_rgh = SolverSetup(solver=Cg(), preconditioner=Jacobi(),
                            convergence=1e-7, relax=1.0, rtol=0.0, atol=1e-9, itmax=5000),
        alpha = SolverSetup(solver=Bicgstab(), preconditioner=Jacobi(),
                            convergence=1e-7, relax=1.0, rtol=0.0, atol=1e-7),
    )

    runtime = Runtime(iterations=ITERATIONS, time_step=dt, write_interval=ITERATIONS)
    config = Configuration(
        solvers=solvers, schemes=schemes, runtime=runtime,
        hardware=hardware, boundaries=BCs)

    GC.gc()
    initialise!(model.fluid.p_rgh, 0.0)
    initialise!(model.momentum.U,  inlet_velocity)
    initialise!(model.fluid.alpha, 0.0)
    setField_Expression!(
        mesh        = mesh,
        field       = model.fluid.alpha,
        condition   = (x, y, z) -> x < X_INTERFACE,
        value_true  = 1.0, value_false = 0.0,
    )

    run!(model, config)

    cells_cpu = Array(mesh.cells)
    alpha_cpu = Array(model.fluid.alpha.values)
    cy_all   = [cell.centre[2] for cell in cells_cpu]
    y_target = cy_all[argmin(abs.(cy_all .- SAMPLING_Y))]

    sample = Tuple{Float64,Float64}[]
    for (i, cell) in enumerate(cells_cpu)
        cx, cy, _ = cell.centre
        if abs(cy - y_target) < Δx/4
            push!(sample, (cx, alpha_cpu[i]))
        end
    end
    sort!(sample, by = first)
    return [s[1] for s in sample], [s[2] for s in sample]
end

profiles = Dict{String, Dict{String, Vector{Float64}}}()
for (label, pure_upwind, c_alpha) in CONFIGS
    @info "Running $label  (pure_upwind=$pure_upwind, cAlpha=$c_alpha)"
    xs, αs = run_and_extract(pure_upwind, c_alpha)
    profiles[label] = Dict("x" => xs, "alpha" => αs)
end

out = "alpha_profiles_quad40.json"
open(out, "w") do io
    JSON.print(io, profiles)
end
@info "Wrote $out"
