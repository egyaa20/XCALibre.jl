# ============================================================
# Stage 1 — Terminal-velocity / drag-closure unit tests
# ============================================================
#
# Purpose
# -------
# Validate the building blocks of the mixture slip-velocity framework in
# `Solvers_2_MULTIPHASE.jl` WITHOUT running a full CFD simulation:
#   * Reference Schiller-Naumann closure (pure math, hand-verifiable)
#   * `compute_Ur!` kernel iterated on a small mesh to a fixed point
#   * `compute_DUmDt!` kernel with known U fields
#   * `blend_properties!` linearity
#   * `turbulent_dispersion!` dispatch (Laminar → no-op)
#
# Physical reference case: 200 µm glass beads in water
#   rho_c = 1000  (water, continuous)
#   rho_d = 2500  (glass beads, dispersed)
#   mu_c  = 1.0e-3
#   d     = 200e-6
#   g     = 9.81 m/s²
#
# Hand calculations (dilute α_d → 0):
#   Stokes:        u_t = g*(ρ_d-ρ_c)*d² / (18*μ_c) ≈ 3.27e-2 m/s  (≈33 mm/s)
#   Schiller-Naumann fixed point (Re ≈ 6.5, f≈1.59): u_t ≈ 2.06e-2 m/s (≈21 mm/s)
#
# IMPORTANT — hardcoded constants to change in Solvers_2_MULTIPHASE.jl
# -------------------------------------------------------------------
# These tests do NOT touch the solver file; they invoke the kernels directly
# with the test-chosen values. The CFD stages (2 & 3) DO require the solver
# constants to be edited. See stage2 / stage3 headers.
#
# Warning about the current `compute_Ur!` formulation
# ---------------------------------------------------
# The kernel uses
#     buoyancy = (ρ_d - ρ_m) / (ρ_m + C_vm*ρ_c)
# rather than the standard Manninen
#     buoyancy_Manninen = (ρ_d - ρ_m) / ρ_d
# At this density ratio and C_vm=0.5 in the dilute limit:
#     solver = (2500-1000)/(1000 + 0.5*1000) = 1.0
#     Manninen = (2500-1000)/2500            = 0.6
# So the solver returns u_t ≈ 1.67× larger than the Manninen form.
# The tests below compare to the *solver's own* formulation, so they
# ought to pass — but they also document the expected Stokes/Manninen
# values side-by-side so the discrepancy is visible.
# ============================================================

using XCALibre
using Test
using LinearAlgebra
using StaticArrays
using KernelAbstractions

"Average of non-zero entries (ignore ghost / untouched cells if any)."
function mean_nonzero(v)
    s, n = 0.0, 0
    @inbounds for x in v
        if abs(x) > 1e-14
            s += x; n += 1
        end
    end
    n == 0 ? 0.0 : s / n
end

# ------------------------------------------------------------
# Pure-Julia reference closure (hand-verifiable)
# ------------------------------------------------------------

"Schiller-Naumann drag factor (Re<1000) clipped at 1."
function sn_factor(Re)
    f = Re < 1000 ? 1 + 0.15*Re^0.687 : 0.0183*Re
    return max(f, 1.0)
end

"""
    stokes_terminal(ρd, ρc, d, μc, g)

Closed-form Stokes terminal velocity (Re→0, α_d→0).
"""
stokes_terminal(ρd, ρc, d, μc, g) = abs(g)*(ρd - ρc)*d^2 / (18*μc)

"""
    sn_terminal_manninen(ρd, ρc, d, μc, g; tol=1e-12, itmax=500)

Fixed-point Schiller-Naumann terminal velocity for *standard Manninen*
buoyancy = (ρd - ρm)/ρd, evaluated at α_d → 0 so ρm = ρc.
Useful as an independent check when the solver's buoyancy is switched to
the Manninen form.
"""
function sn_terminal_manninen(ρd, ρc, d, μc, g; tol=1e-12, itmax=500)
    ρm = ρc                                    # dilute limit
    τd = ρd * d^2 / (18*μc)
    buoy = (ρd - ρm) / ρd
    u = stokes_terminal(ρd, ρc, d, μc, g)
    for _ in 1:itmax
        Re = ρc*abs(u)*d/μc
        f  = sn_factor(Re)
        u_new = (τd / f) * buoy * abs(g)
        abs(u_new - u) < tol && return u_new
        u = u_new
    end
    return u
end

"""
    sn_terminal_solver(ρd, ρc, d, μc, g, Cvm; tol=1e-12, itmax=500)

Fixed-point Schiller-Naumann terminal velocity for the *solver's* buoyancy
form (ρd - ρm)/(ρm + Cvm·ρc), evaluated at α_d → 0 (ρm = ρc).
"""
function sn_terminal_solver(ρd, ρc, d, μc, g, Cvm; tol=1e-12, itmax=500)
    ρm = ρc
    τd = ρd * d^2 / (18*μc)
    buoy = (ρd - ρm) / (ρm + Cvm*ρc)
    u = stokes_terminal(ρd, ρc, d, μc, g)
    for _ in 1:itmax
        Re = ρc*abs(u)*d/μc
        f  = sn_factor(Re)
        u_new = (τd / f) * buoy * abs(g)
        abs(u_new - u) < tol && return u_new
        u = u_new
    end
    return u
end


@testset "Stage 1 — terminal velocity & drag closure" begin

    # Canonical inputs (200 µm glass beads in water)
    ρd, ρc, μc, d, gmag, Cvm = 2500.0, 1000.0, 1.0e-3, 2.0e-4, 9.81, 0.5

    @testset "Reference closures (hand calc)" begin
        u_st = stokes_terminal(ρd, ρc, d, μc, gmag)
        @test isapprox(u_st, 3.27e-2; rtol=2e-2)   # ≈33 mm/s

        u_sn_m = sn_terminal_manninen(ρd, ρc, d, μc, gmag)
        @test 1.5e-2 ≤ u_sn_m ≤ 3.0e-2             # SN < Stokes, ≈20–25 mm/s

        u_sn_s = sn_terminal_solver(ρd, ρc, d, μc, gmag, Cvm)
        @test u_sn_s > u_sn_m                      # solver's form is larger
        @test u_sn_s < stokes_terminal(ρd, ρc, d, μc, gmag) * 2

        # Self-consistency: recompute residual of fixed point
        τd   = ρd * d^2 / (18*μc)
        Re   = ρc*u_sn_s*d/μc
        buoy = (ρd - ρc) / (ρc + Cvm*ρc)
        @test isapprox(u_sn_s, τd/sn_factor(Re) * buoy * gmag; rtol=1e-6)
    end

    @testset "sn_factor monotone & clipped" begin
        @test sn_factor(0.0) == 1.0
        @test sn_factor(1e-6) ≈ 1.0 atol=1e-4   # floor approached, not exact
        @test sn_factor(1.0) > 1.0
        @test sn_factor(100.0) > sn_factor(1.0)
    end

    # ------------------------------------------------------------
    # Kernel-level tests — invoke actual solver kernels on a mesh
    # ------------------------------------------------------------
    grids_dir = pkgdir(XCALibre, "examples/0_GRIDS")
    mesh_file = joinpath(grids_dir, "unit_test_stage2.unv")
    mesh  = UNV2D_mesh(mesh_file, scale=1.0)
    backend = CPU()
    workgroup = AutoTune()
    activate_multithread(backend)
    mesh_dev = adapt(backend, mesh)
    hardware = Hardware(backend=backend, workgroup=workgroup)

    # Minimal config stub carrying hardware
    config = (hardware=hardware,)

    @testset "blend_properties! linearity" begin
        alpha = ScalarField(mesh_dev)
        prop  = ScalarField(mesh_dev)

        # α = 0 → property == prop_1 (secondary)
        initialise!(alpha, 0.0)
        XCALibre.Solvers.blend_properties!(prop, alpha, 1000.0, 1.225)
        @test all(prop.values .≈ 1.225)

        # α = 1 → property == prop_0 (main)
        initialise!(alpha, 1.0)
        XCALibre.Solvers.blend_properties!(prop, alpha, 1000.0, 1.225)
        @test all(prop.values .≈ 1000.0)

        # α = 0.5 → midpoint
        initialise!(alpha, 0.5)
        XCALibre.Solvers.blend_properties!(prop, alpha, 1000.0, 1.225)
        @test all(prop.values .≈ 0.5*(1000.0 + 1.225))
    end

    @testset "compute_DUmDt! with constant U and zero gradU → zero" begin
        U      = VectorField(mesh_dev)
        U_prev = VectorField(mesh_dev)
        gradU  = Grad{Gauss}(U)
        DUmDt  = VectorField(mesh_dev)

        initialise!(U,      [0.0, 0.1, 0.0])
        initialise!(U_prev, [0.0, 0.1, 0.0])     # ∂U/∂t = 0
        # gradU defaults to zero after construction

        XCALibre.Solvers.compute_DUmDt!(DUmDt, U, U_prev, gradU, 1.0e-3, config)
        @test maximum(abs, DUmDt.x.values) < 1e-12
        @test maximum(abs, DUmDt.y.values) < 1e-12
        @test maximum(abs, DUmDt.z.values) < 1e-12
    end

    @testset "compute_DUmDt! with a dU/dt jump" begin
        U      = VectorField(mesh_dev)
        U_prev = VectorField(mesh_dev)
        gradU  = Grad{Gauss}(U)
        DUmDt  = VectorField(mesh_dev)

        initialise!(U,      [0.0, 0.2, 0.0])
        initialise!(U_prev, [0.0, 0.1, 0.0])
        dt = 0.01
        XCALibre.Solvers.compute_DUmDt!(DUmDt, U, U_prev, gradU, dt, config)
        expected_dudt = (0.2 - 0.1) / dt  # = 10.0
        @test all(isapprox.(DUmDt.y.values, expected_dudt; rtol=1e-10))
    end

    @testset "compute_Ur! fixed-point → solver terminal velocity" begin
        # Fields
        alpha = ScalarField(mesh_dev)
        rho   = ScalarField(mesh_dev)
        Ur    = VectorField(mesh_dev)
        DUmDt = VectorField(mesh_dev)

        initialise!(alpha, 1.0 - 1e-3)   # dilute dispersed phase
        # After blend_properties!(rho, alpha, rho_c, rho_d) with α ≈ 1 → ρ_m ≈ ρ_c
        XCALibre.Solvers.blend_properties!(rho, alpha, ρc, ρd)
        initialise!(Ur,    [0.0, 0.0, 0.0])
        initialise!(DUmDt, [0.0, 0.0, 0.0])

        g_vec = SVector{3,Float64}(0.0, -gmag, 0.0)
        τd    = ρd*d^2 / (18*μc)

        # Iterate kernel to fixed point (Ur depends on |Ur| via Re → Picard)
        for _ in 1:200
            XCALibre.Solvers.compute_Ur!(Ur, alpha, rho, g_vec, DUmDt,
                                         ρc, ρd, μc, d, τd, Cvm, config)
        end

        u_expected = sn_terminal_solver(ρd, ρc, d, μc, gmag, Cvm)
        uy = Ur.y.values
        # All interior cells should have the same downward terminal velocity
        # (negative because g is -y and u_d,rel is downward). Note that Ur is
        # the dispersed-phase relative velocity, scaled by 1/α_c ≈ 1.
        @test maximum(abs, Ur.x.values) < 1e-10
        @test maximum(abs, Ur.z.values) < 1e-10
        @test isapprox(-u_expected, mean_nonzero(uy); rtol=5e-2)

        # Stokes upper bound (solver's form happens to exceed Stokes)
        @test abs(mean_nonzero(uy)) < 10*stokes_terminal(ρd, ρc, d, μc, gmag)

        @info "compute_Ur! fixed point" u_computed=abs(mean_nonzero(uy)) u_solver_expected=u_expected u_stokes=stokes_terminal(ρd, ρc, d, μc, gmag)
    end # compute_Ur

    @testset "turbulent_dispersion! dispatch" begin
        # Laminar branch must be a no-op
        alpha = ScalarField(mesh_dev)
        Ur    = VectorField(mesh_dev)
        ∇alpha = Grad{Gauss}(alpha)
        initialise!(alpha, 0.5)
        initialise!(Ur, [1.0, 2.0, 3.0])
        # ∇alpha defaults to zero

        turb_laminar = Laminar()
        XCALibre.Solvers.turbulent_dispersion!(Ur, alpha, ∇alpha, turb_laminar, 0.7, config)
        @test all(Ur.x.values .≈ 1.0)
        @test all(Ur.y.values .≈ 2.0)
        @test all(Ur.z.values .≈ 3.0)
    end
end
