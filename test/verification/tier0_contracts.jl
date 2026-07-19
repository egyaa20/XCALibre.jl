# =============================================================================
# Tier 0 — Property DATA CONTRACTS (Class A: formulation-independent)
#
# Tests the property source through its public lookup interface. Nothing here
# depends on which PDE terms consume which column — if the pressure equation
# changes shape (ψ-split, lumped ∂ρ/∂t, implicit, …), every test below still
# runs unchanged and still means the same thing.
# =============================================================================

T_CRIT_N2 = 126.192

vv_table = HelmholtzTable(
    N2();
    p_operating = 3.5e6,
    T_range     = (75.0, 375.0),
    p_range     = (3.4e6, 3.6e6),
    n_T         = 300,
    n_p         = 6,
    T_snapshot  = 78.0,
)
vv_eos = HelmholtzEnergy(name = N2())

@testset "Tier 0 — property data contracts" begin

    # -------------------------------------------------------------------------
    # 0.1 Interpolation fidelity: table lookup vs direct EOS at OFF-GRID (T,p).
    # Tolerances are interpolation-error budgets for the 300-pt base grid:
    # away from the cp spike the properties are smooth (linear interp ~1e-4..1e-3);
    # inside T_c ± 2 K curvature dominates and a linear table is honest to ~15 %.
    # -------------------------------------------------------------------------
    @testset "0.1 lookup vs direct EOS" begin
        far_T  = [83.3, 101.7, 118.9, 139.3, 187.4, 291.6]
        near_T = [125.3, 126.8]
        ps     = [3.43e6, 3.512e6, 3.578e6]

        for p in ps, (Ts, rtol) in ((far_T, 2e-3), (near_T, 0.15))
            for T in Ts
                lut = lookup_helmholtz(vv_table, T, p)
                ref = vv_eos(T, p)
                @test isapprox(lut.rho, ref.rho[1]; rtol=rtol)
                @test isapprox(lut.mu,  ref.mu[1];  rtol=max(rtol, 5e-3))
                @test isapprox(lut.cp,  ref.cp[1];  rtol=max(rtol, 1e-2))
                @test isapprox(lut.k,   ref.k[1];   rtol=max(rtol, 5e-3))
            end
        end
    end

    # -------------------------------------------------------------------------
    # 0.2 Thermodynamic self-consistency BETWEEN exposed quantities, by finite
    # differences ON THE TABLE ITSELF. Catches unit swaps / return-order swaps
    # anywhere in the chain (e.g. the historical β↔s functor swap) without
    # caring how the table is built.
    # -------------------------------------------------------------------------
    @testset "0.2 cross-derivative identities" begin
        # cp ≈ ∂h/∂T at fixed p (away from the spike)
        for p in (3.45e6, 3.5e6), T in (85.0, 105.0, 145.0, 200.0, 300.0)
            δ = 0.25
            dhdT = (Ttoh_helmholtz(vv_table, T + δ, p) -
                    Ttoh_helmholtz(vv_table, T - δ, p)) / (2δ)
            cp = lookup_helmholtz(vv_table, T, p).cp
            @test isapprox(dhdT, cp; rtol=1e-2)
        end

        # h strictly monotone in T; T ↔ h round trip exact to interpolation
        for p in (3.42e6, 3.5e6, 3.58e6)
            hs = [Ttoh_helmholtz(vv_table, T, p) for T in 76.0:3.0:374.0]
            @test all(diff(hs) .> 0)
            for T in 76.0:24.5:374.0
                @test abs(htoT_helmholtz(vv_table, Ttoh_helmholtz(vv_table, T, p), p) - T) < 1e-8
            end
        end

        # ψ ≈ ∂ρ/∂p|T and (∂ρ/∂T)|p — auto-skipped if the source stops
        # exposing them (they are OPTIONAL columns; the continuity contract
        # tests in tier 5 are the formulation-level check).
        probe = lookup_helmholtz(vv_table, 150.0, 3.5e6)
        if hasproperty(probe, :psi)
            for T in (90.0, 150.0, 250.0), p in (3.46e6, 3.54e6)
                δp = 5.0e3
                fd = (lookup_helmholtz(vv_table, T, p + δp).rho -
                      lookup_helmholtz(vv_table, T, p - δp).rho) / (2δp)
                ψ = lookup_helmholtz(vv_table, T, p).psi
                @test ψ > 0
                @test isapprox(fd, ψ; rtol=2e-2)
            end
        end
        if hasproperty(probe, :drhodT)
            for T in (90.0, 150.0, 250.0)
                δT = 0.3
                fd = (lookup_helmholtz(vv_table, T + δT, 3.5e6).rho -
                      lookup_helmholtz(vv_table, T - δT, 3.5e6).rho) / (2δT)
                dρdT = lookup_helmholtz(vv_table, T, 3.5e6).drhodT
                @test dρdT < 0
                @test isapprox(fd, dρdT; rtol=5e-2)
            end
        end
    end

    # -------------------------------------------------------------------------
    # 0.3 Poisoned-input safety. Contract: no lookup may ever throw or index
    # out of bounds (memory safety — the illegal-GPU-address class); lookups
    # keyed on T/p must return FINITE edge-clamped values. h-inversion with a
    # NaN enthalpy may propagate NaN (garbage in) but must not crash.
    # -------------------------------------------------------------------------
    @testset "0.3 poisoned inputs" begin
        bad = (NaN, Inf, -Inf, -1.0e300)
        for x in bad
            r1 = lookup_helmholtz(vv_table, x, 3.5e6)   # poisoned T
            @test isfinite(r1.rho) && isfinite(r1.cp)
            r2 = lookup_helmholtz(vv_table, 100.0, x)    # poisoned p
            @test isfinite(r2.rho) && isfinite(r2.cp)
            @test isfinite(Ttoh_helmholtz(vv_table, x, 3.5e6))
            @test (htoT_helmholtz(vv_table, x, 3.5e6); true)   # must not throw
            @test (htoT_helmholtz(vv_table, 1.0e5, x); true)   # must not throw
        end
    end

    # -------------------------------------------------------------------------
    # 0.4 Grid contract: strictly increasing; refinement expressed as RATIOS of
    # local spacing vs the far-field base (constructor promise: ×5 within
    # T_c ± 1 K, ×3 within T_c ± 10 K), tolerant of segment rounding.
    # -------------------------------------------------------------------------
    @testset "0.4 refined-grid contract" begin
        g = vv_table.T_grid
        d = diff(g)
        @test all(d .> 0)
        mids = 0.5 .* (g[1:end-1] .+ g[2:end])
        base = median(d[mids .> 160.0])
        fine = median(d[abs.(mids .- T_CRIT_N2) .< 0.9])
        mid  = median(d[2.0 .< abs.(mids .- T_CRIT_N2) .< 9.0])
        @test vv_check("fine-band refinement ratio", base / fine; between=(3.5, 7.0))
        @test vv_check("mid-band refinement ratio",  base / mid;  between=(2.2, 4.5))
    end

    # -------------------------------------------------------------------------
    # 0.5 Blend identities (value conventions: phase-1 is the α=1 phase).
    # Guarded — the blend helpers are solver-internal and may be renamed.
    # -------------------------------------------------------------------------
    @testset "0.5 blend identities" begin
        if isdefined(XCALibre.Solvers, :blend_properties!)
            mesh_b  = UNV2D_mesh(vv_grid("quad20.unv"), scale=0.001)
            fld     = ScalarField(mesh_b)
            alpha_b = ScalarField(mesh_b)
            for (a, expected) in ((1.0, 1000.0), (0.0, 1.2), (0.5, 500.6))
                fill!(alpha_b.values, a)
                XCALibre.Solvers.blend_properties!(fld, alpha_b, 1000.0, 1.2)
                @test all(isapprox.(fld.values, expected; rtol=1e-14))
            end
        else
            @test_skip "blend_properties! not defined in current solver"
        end
    end
end
