using Test, StaticArrays, LinearAlgebra

# Convention:
#   phase[1] = water = continuous,  rho_1 ≈ 1000, mu_1 ≈ 1e-3
#   phase[2] = air   = dispersed,   rho_2 ≈ 1.2
#   alpha = 1.0 → pure water,  alpha = 0.0 → pure air
#   rho_m = rho_1*alpha + rho_2*(1 - alpha)

function _scalar_blend(alpha, rho_1, rho_2)
    return rho_1 * alpha + rho_2 * (1.0 - alpha)
end

function _scalar_DUmDt(U, U_prev, gradU, dt)
    dUdt = (U - U_prev) / dt

    ux, uy, uz = U
    dudx, dudy, dudz = gradU[1,:]
    dvdx, dvdy, dvdz = gradU[2,:]
    dwdx, dwdy, dwdz = gradU[3,:]

    conv_x = ux*dudx + uy*dudy + uz*dudz
    conv_y = ux*dvdx + uy*dvdy + uz*dvdz
    conv_z = ux*dwdx + uy*dwdy + uz*dwdz

    return dUdt + SVector(conv_x, conv_y, conv_z)
end

function _scalar_Ur(rho_m, rho_1, rho_2, mu_1, d, g, DUmDt, Ur_prev)
    # rho_1 = water (continuous), rho_2 = air (dispersed)
    rho_c = rho_1
    rho_d = rho_2
    mu_c  = mu_1

    tau_d = (rho_d * d^2) / (18 * mu_c)

    Ur_mag = norm(Ur_prev)
    Re_p   = rho_c * Ur_mag * d / mu_c

    f_drag = if Re_p < 1000
        1.0 + 0.15 * Re_p^0.687
    else
        0.0183 * Re_p
    end
    f_drag = max(f_drag, 1.0)

    a_eff    = g - DUmDt
    buoyancy = (rho_d - rho_m) / rho_d

    return (tau_d / f_drag) * buoyancy * a_eff
end

# ─────────────────────────────────────────────────────────────
# Physical constants — matches your framework
# ─────────────────────────────────────────────────────────────
const RHO_1  = 1000.0          # water density
const RHO_2  = 1.2             # air density
const MU_1   = 1e-3            # water dynamic viscosity
const D      = 1e-3            # bubble diameter
const G      = SVector(0.0, -9.81, 0.0)

@testset "blend_properties convention" begin

    @testset "alpha=1 → pure water" begin
        rho_m = _scalar_blend(1.0, RHO_1, RHO_2)
        @test rho_m ≈ RHO_1
    end

    @testset "alpha=0 → pure air" begin
        rho_m = _scalar_blend(0.0, RHO_1, RHO_2)
        @test rho_m ≈ RHO_2
    end

    @testset "alpha=0.5 → arithmetic mean" begin
        rho_m = _scalar_blend(0.5, RHO_1, RHO_2)
        @test rho_m ≈ 0.5 * RHO_1 + 0.5 * RHO_2
    end
end

@testset "DUmDt" begin

    @testset "Zero acceleration — steady uniform flow" begin
        U      = SVector(1.0, 0.0, 0.0)
        U_prev = SVector(1.0, 0.0, 0.0)
        gradU  = zeros(3, 3)
        @test _scalar_DUmDt(U, U_prev, gradU, 0.01) ≈ SVector(0.0, 0.0, 0.0) atol=1e-14
    end

    @testset "Pure temporal acceleration, no convection" begin
        # dU/dt = (2-1)/0.5 = 2 in x
        result = _scalar_DUmDt(SVector(2.0,0.0,0.0), SVector(1.0,0.0,0.0), zeros(3,3), 0.5)
        @test result ≈ SVector(2.0, 0.0, 0.0) atol=1e-14
    end

    @testset "Pure convective acceleration, steady flow" begin
        u     = 3.0
        gradU = [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        result = _scalar_DUmDt(SVector(u,0.0,0.0), SVector(u,0.0,0.0), gradU, 0.01)
        @test result ≈ SVector(u, 0.0, 0.0) atol=1e-13
    end

    @testset "Superposition of temporal + convective" begin
        gradU = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 0.0]
        # dU/dt=(1,0,0), (U·∇)U=(2,2,0)
        result = _scalar_DUmDt(SVector(2.0,1.0,0.0), SVector(1.0,1.0,0.0), gradU, 1.0)
        @test result ≈ SVector(3.0, 2.0, 0.0) atol=1e-13
    end
end

@testset "compute_Ur" begin

    @testset "Direction — air bubbles must rise in water (alpha≈1 region)" begin
        # alpha=1 → pure water, single air bubble rising
        rho_m   = _scalar_blend(1.0, RHO_1, RHO_2)  # = 1000 kg/m³
        Ur      = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0,0.0,0.0), SVector(0.0,0.0,0.0))

        # buoyancy = (rho_2 - rho_m)/rho_2 = (1.2-1000)/1.2 < 0
        # a_eff_y  = -9.81
        # Ur_y     = tau_d * negative * negative = positive → bubble moves UP ✓
        @test Ur[2] > 0.0
        println("alpha=1.0 (pure water): Ur_y = $(Ur[2]) m/s  (>0 = rising ✓)")
    end

    @testset "Direction — water drops must sink in air (alpha≈0 region)" begin
        # alpha=0 → pure air, water droplet sinking
        rho_m   = _scalar_blend(0.0, RHO_1, RHO_2)  # = 1.2 kg/m³
        Ur      = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0,0.0,0.0), SVector(0.0,0.0,0.0))

        # buoyancy = (1.2 - 1.2)/1.2 = 0 exactly in pure air
        # (dispersed phase IS the mixture → no slip)
        @test abs(Ur[2]) < 1e-10
        println("alpha=0.0 (pure air): Ur_y = $(Ur[2]) m/s  (≈0 = no slip in pure dispersed ✓)")
    end

    @testset "Stokes limit — Re=0, analytical terminal velocity in pure water" begin
        # alpha=1.0: single air bubble in water, Re→0 → f_drag=1
        rho_m    = _scalar_blend(1.0, RHO_1, RHO_2)
        Ur       = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0,0.0,0.0), SVector(0.0,0.0,0.0))

        # Stokes: v_t = d²(ρ_2 - ρ_1)·g / (18·μ_1)  [negative = downward for drops, positive for bubbles]
        Ur_stokes = D^2 * (RHO_2 - RHO_1) * (-9.81) / (18 * MU_1)
        println("Stokes terminal velocity: Ur_y=$(Ur[2])  expected=$(Ur_stokes)")

        @test Ur[2] ≈ Ur_stokes rtol=1e-10
        @test Ur[2] > 0.0   # must be upward
    end

    @testset "DUmDt — downward mixture acceleration reduces bubble slip" begin
        rho_m   = _scalar_blend(1.0, RHO_1, RHO_2)
        Ur_prev = SVector(0.0, 0.0, 0.0)

        # a_eff = g - DUmDt
        # DUmDt downward (-5) → a_eff_y = -9.81+5 = -4.81  → smaller |Ur|
        # DUmDt upward  (+5) → a_eff_y = -9.81-5 = -14.81 → larger  |Ur|
        Ur_base = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0,  0.0,0.0), Ur_prev)
        Ur_down = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0, -5.0,0.0), Ur_prev)
        Ur_up   = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0, +5.0,0.0), Ur_prev)

        @test abs(Ur_down[2]) < abs(Ur_base[2])
        @test abs(Ur_up[2])   > abs(Ur_base[2])

        tau_d = (RHO_2 * D^2) / (18 * MU_1)
        b     = (RHO_2 - RHO_1) / RHO_2
        @test Ur_down[2] ≈ tau_d * b * (-9.81 + 5.0) rtol=1e-10
        @test Ur_up[2]   ≈ tau_d * b * (-9.81 - 5.0) rtol=1e-10

        println("Ur_down=$(Ur_down[2])  Ur_base=$(Ur_base[2])  Ur_up=$(Ur_up[2])")
    end

    @testset "Wrong denominator regression — rho_m vs rho_d gives 833x error" begin
        rho_m  = _scalar_blend(1.0, RHO_1, RHO_2)  # alpha=1, pure water
        tau_d  = (RHO_2 * D^2) / (18 * MU_1)

        Ur_correct = tau_d * (RHO_2 - rho_m) / RHO_2  * (-9.81)  # denominator = rho_2 (air)
        Ur_wrong   = tau_d * (RHO_2 - rho_m) / rho_m  * (-9.81)  # denominator = rho_m (water)

        @test Ur_correct / Ur_wrong ≈ rho_m / RHO_2 rtol=1e-10   # ≈ 833
        @test abs(Ur_correct) > 100 * abs(Ur_wrong)
        println("Denominator bug ratio = $(Ur_correct/Ur_wrong)x  (≈ $(rho_m/RHO_2))")
    end

    @testset "Schiller-Naumann — high Re transitions to linear drag" begin
        rho_m   = _scalar_blend(1.0, RHO_1, RHO_2)
        Ur_prev = SVector(10.0, 0.0, 0.0)   # Re = 1000*10*1e-3/1e-3 = 10000 → linear regime

        Re_p = RHO_1 * 10.0 * D / MU_1
        @test Re_p > 1000

        tau_d         = (RHO_2 * D^2) / (18 * MU_1)
        f_drag_linear = 0.0183 * Re_p
        buoyancy      = (RHO_2 - rho_m) / RHO_2
        Ur_expected_y = (tau_d / f_drag_linear) * buoyancy * (-9.81)

        Ur = _scalar_Ur(rho_m, RHO_1, RHO_2, MU_1, D, G, SVector(0.0,0.0,0.0), Ur_prev)
        @test Ur[2] ≈ Ur_expected_y rtol=1e-10
    end
end