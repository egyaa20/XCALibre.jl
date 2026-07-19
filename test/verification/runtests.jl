# =============================================================================
# V&V platform entry point — see DESIGN.md
#
# Standalone:
#     julia --project=. test/verification/runtests.jl            # default tiers
#     julia --project=. test/verification/runtests.jl 0,3        # explicit tiers
#     XCAL_VV_TIERS=0,3,5 julia --project=. test/verification/runtests.jl
#
# Tiers: 0 property contracts | 1 operators | 2 equations (MMS) | 3 coupling |
#        4 mixture | 5 energy+density | 6 boiling | 7 validation/robustness
# Unimplemented tiers are silently absent; tier 5 self-skips until the
# multiphase solver regains energy support.
# =============================================================================

using XCALibre
using Test
using LinearAlgebra
using Statistics

vv_tiers = let
    raw = !isempty(ARGS) ? ARGS[1] : get(ENV, "XCAL_VV_TIERS", "0,2,3,4,5")
    Set(strip.(split(raw, ",")))
end

include("common.jl")

@testset verbose = true "V&V platform (tiers: $(join(sort(collect(vv_tiers)), ",")))" begin
    tier_files = Dict(
        "0" => "tier0_contracts.jl",
        "2" => "tier2_equations.jl",
        "3" => "tier3_coupling.jl",
        "4" => "tier4_mixture.jl",
        "5" => "tier5_energy_density.jl",
    )
    for tier in sort(collect(vv_tiers))
        file = get(tier_files, tier, nothing)
        if file === nothing
            @warn "V&V tier $tier not implemented yet — skipping"
            continue
        end
        include(file)
    end
end
