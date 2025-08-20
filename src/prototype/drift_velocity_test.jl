using XCALibre
# using Accessors
using KernelAbstractions
using StaticArrays
using LinearAlgebra



grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")

# grid = "finer_mesh_laplace.unv"
grid = "laplace_unit_3by3.unv"

mesh_file = joinpath(grids_dir, grid)

mesh = UNV2D_mesh(mesh_file)


backend = CPU(); workgroup = 1024; activate_multithread(backend)

hardware = Hardware(backend=backend, workgroup=workgroup)

mesh_dev = adapt(backend, mesh)



model = Physics(
    time = Transient(),
    fluid = Fluid{Multiphase}(
        phases = ( #first phase is liquid, second if vapour - common assumption
            Phase(eos=ConstEos(1.0), mu=ConstMu(1.8e-5)),       #air
            Phase(eos=ConstEos(1000.0), mu=ConstMu(1.0e-3))     #water
        ),
    ),
    turbulence = RANS{Laminar}(),
    energy = Energy{Isothermal}(),
    domain = mesh_dev
    )

BCs = assign(
    region = mesh_dev,
    (
        U = [     
            Dirichlet(:left_wall, 50.0),
            Zerogradient(:right_wall),
            Dirichlet(:bottom_wall, 10.0),
            Zerogradient(:upper_wall)
        ],
        p = [     
            Dirichlet(:left_wall, 50.0),
            Zerogradient(:right_wall),
            Dirichlet(:bottom_wall, 10.0),
            Zerogradient(:upper_wall)
        ],
        alpha = [     
            Dirichlet(:left_wall, 50.0),
            Zerogradient(:right_wall),
            Dirichlet(:bottom_wall, 10.0),
            Zerogradient(:upper_wall)
        ],
    )
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
    U = Schemes(time=Euler, divergence = Upwind),
    p = Schemes(time=Euler),
    alpha = Schemes(time=Euler, divergence = Upwind)
)

iterations=10

runtime = Runtime(
    iterations=iterations, write_interval=1, time_step=1)

config = Configuration(
    solvers=solvers, schemes=schemes, runtime=runtime, hardware=hardware, boundaries=BCs)






# function construct_acceleration_field(U_m, ∇U_m, U_m_prev, dt) #, workgroup
#     grad_U = ∇U_m.result
#     a = VectorField(mesh)

#     # g = props.DriftVelocity.gravity
#     g = [0.0, -9.81, 0.0]

#     x0, y0, z0 = g[1], g[2], g[3]

#     x = ScalarField(mesh)
#     y = ScalarField(mesh)
#     z = ScalarField(mesh)
#     initialise!(x, x0)
#     initialise!(y, y0)
#     initialise!(z, z0)
#     G = VectorField(x, y, z, mesh)

#     dUdt = VectorField(mesh)
    
#     for i in eachindex(U_m)
#         dUdt[i] = (U_m[i] - U_m_prev[i]) / dt
#         a[i] = G[i] - (grad_U[i] * U_m[i]) - dUdt[i]
#     end
    
#     return a
# end


# function get_E!(E, rho_p, d_p, mu_p, rho, a_field)
#     rho_p = rho_p.values
#     d_p = d_p.values
#     mu_p = mu_p.values
#     rho = rho.values

#     for i in eachindex(E)
#         E[i] = (rho_p[i]*((d_p[i])^2))/(18*mu_p[i]) * ((rho_p[i]-rho[i])/rho_p[i]) * a_field[i]
#     end

#     return nothing
# end

# function compute_RE!(RE, rho_q, v_pq, d_p, mu_q)
#     rho_q = rho_q.values
#     d_p = d_p.values
#     mu_q = mu_q.values
#     RE = RE.values
    
#     for i in eachindex(RE)
#         RE[i] = (rho_q[i] * d_p[i] * norm(v_pq[i]))/mu_q[i]
#     end

#     return nothing
# end

# function find_Vpq!(rho_q, rho_p, rho, mu_q, mu_p, d_p, a_field, v_pq_prev, isInitialisation, v_low, v_high, v_pq, mesh)
#     @. mu_q.values = mu_q.values * rho_q.values #conversion again..... NU -> MU
#     @. mu_p.values = mu_p.values * rho_p.values #conversion again.....
    
#     E = VectorField(mesh)
#     expr = ScalarField(mesh)
#     v_pq_trial = VectorField(mesh)
#     RE_trial = ScalarField(mesh)

#     get_E!(E, rho_p, d_p, mu_p, rho, a_field)

#     @. expr.values = mu_q.values/(0.0183 * rho_q.values * d_p.values)
    
#     for i in eachindex(v_pq_trial)
#         v_pq_trial[i] = sqrt.(E[i]*expr.values[i])
#     end

#     compute_RE!(RE_trial, rho_q, v_pq_trial, d_p, mu_q)

#     K = ScalarField(mesh)
#     @. K.values = 0.15 * ((rho_q.values * d_p.values)/mu_q.values)^0.687

#     for i in eachindex(RE_trial)
#         if RE_trial[i] > 1000
#             v_pq[i] = v_pq_trial[i]
#         else
#             if isInitialisation
#                 v_pq[i] = Vpq_bisect(K[i], E[i], v_low, v_high)
#             else
#                 v_pq[i] = Vpq_newton_raphson(E[i], K[i], v_pq_prev[i])
#             end
#         end
#     end
# end



# function Vpq_newton_raphson(E, K, v; max_iter=20, tol=1.0e-7)
#     for it in 1:max_iter # E IS A VECTOR FIELD
#         v_pow_0_687 = v .^ 0.687
#         f = v .+ K .* v_pow_0_687 .* v .- E
#         # f = v + K*(v^1.687) - E

#         f_prime = 1.0 .+ 1.687 .* K .* v_pow_0_687
#         # f_prime = 1 + 1.687 * K * (v^1.687)

#         v_new = v .- f ./ f_prime # division by zero in one direction is possible
#         # v_new = v - (f/f_prime)

#         err = abs.((v_new .- v) ./ (v_new .+ eps())) # Add eps to prevent division by zero
#         # err = abs((v_new - v) / v_new)

#         if maximum(err) < tol
#             v_pq = v_new
#             println(v_pq)

#             return v_pq
#         end

#         v = v_new
#     end
# end


# function bisect_f(v, K, E)
#     return v .+ K .* (v .^ 1.687) .- E
# end

# function Vpq_bisect(K, E, v_low, v_high; max_iter=50, tol=1.0e-7) #v_mid=SVector(0.0, 0.0, 0.0)
#     for it in 1:max_iter
#         v_mid = (v_low .+ v_high) ./ 2.0
        
#         f_mid = bisect_f(v_mid, K, E)
#         f_low = bisect_f(v_low, K, E)

#         vectorised_check = sign.(f_mid) .== sign.(f_low)
#         v_low = ifelse.(vectorised_check, v_mid, v_low)
#         v_high = ifelse.(vectorised_check, v_high, v_mid)

#         # if sign.(f_mid) == sign.(f_low)
#         #     v_low = v_mid
#         # else
#         #     v_high = v_mid
#         # end

#         err = abs.(v_high .- v_low)

#         if maximum(err) < tol
#             v_pq = v_mid

#             return v_pq
#         end
#     end
# end

# function update_velocities!(U_m, v_pq, v_dr_p, v_dr_q, v_p, v_q, rho, rho_q, rho_p, alpha, mesh)
#     alpha = alpha.values
#     rho = rho.values
#     rho_q = rho_q.values
#     rho_p = rho_p.values

#     C_p = ScalarField(mesh)
#     C_q = ScalarField(mesh)

#     @. C_p.values = (alpha * rho_p) / rho          #liquid
#     @. C_q.values = ((1.0-alpha) * rho_q) / rho    #vapour

    
#     for i in eachindex(U_m)
#         v_dr_p[i] = (1.0 + C_p.values[i]) * v_pq[i]
#         v_dr_q[i] = -(1.0 + C_q.values[i]) * v_pq[i]

#         v_p[i] = v_dr_p[i] - U_m[i]
#         v_q[i] = v_dr_q[i] - U_m[i]
#     end
# end




dt = 0.5

U_m = VectorField(mesh)
U_mf = FaceVectorField(mesh)
U_m_prev = VectorField(mesh)
v_pq_prev = VectorField(mesh)
v_pq = VectorField(mesh)


initialise!(U_m, [5.0, 0.0, 0.0])
initialise!(U_m_prev, [4.0, 0.0, 0.0])
initialise!(v_pq_prev, [0.01, 0.0, 0.0])

∇U_m = Grad{schemes.U.gradient}(U_m)

a_field = construct_acceleration_field(U_m, ∇U_m, U_m_prev, dt)

rho_q = ScalarField(mesh)
rho_p = ScalarField(mesh)
mu_q = ScalarField(mesh)
mu_p = ScalarField(mesh)

alpha = ScalarField(mesh)
initialise!(alpha, 0.5)

isInitialisation = true
d_p = ScalarField(mesh)
v_pq = VectorField(mesh)
# mu_p = ScalarField(mesh)

initialise!(rho_q, 0.1)
initialise!(rho_p, 0.2)
initialise!(mu_q, 1.0e-5)
initialise!(mu_p, 5.0e-5)
initialise!(d_p, 1.0e-5)

rho = ScalarField(mesh)
E = VectorField(mesh)
initialise!(rho, 0.25)
initialise!(E, [1.0, 1.0, 1.0])

d_p[2] = 1.0


find_Vpq!(rho_q, rho_p, rho, mu_q, mu_p, d_p, a_field, v_pq_prev, isInitialisation, SVector(1.0e-9, 1.0e-9, 1.0e-9), SVector(10.0, 10.0, 10.0), v_pq, mesh)


v_dr_p = VectorField(mesh)
v_dr_q = VectorField(mesh)
v_p = VectorField(mesh)
v_q = VectorField(mesh)

update_velocities!(U_m, v_pq, v_dr_p, v_dr_q, v_p, v_q, rho, rho_q, rho_p, alpha, mesh)
