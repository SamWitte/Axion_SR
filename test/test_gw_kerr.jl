# Relativistic GW transition rates of the cloud (src/Numerics/gw_kerr.jl) against the
# non-relativistic multipole rates of gw_rates.jl in the α -> 0 limit.
#   julia test/test_gw_kerr.jl
using Test, Suppressor
include(joinpath(@__DIR__, "..", "src", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "src", "state_utils.jl"))
@suppress include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_rates.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "teukolsky.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "kerr_modes.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_kerr.jl"))

# Kerr/NR ratio of a transition, with the NR rate evaluated at the Leaver splitting
function ratio_to_nr(sa, sb, alph, a)
    ma, _ = kerr_mode(sa..., alph, a); mb, _ = kerr_mode(sb..., alph, a)
    Γ, P, modes, ω̃, m̃ = gw_kerr_rate(ma, mb, :tr)
    return Γ / gw_transition_rate(gw_transition_channels(sa, sb), alph, ω̃ / alph), modes, P
end

@testset "GW transitions on Kerr" begin
    @testset "Schwarzschild: Γ_Kerr/Γ_NR = 1 - c α², c constant" begin
        # (channel, coefficient c found numerically; the test checks the α² scaling and the limit)
        for (sa, sb) in (((3, 2, 2), (2, 1, 1)), ((4, 3, 3), (2, 1, 1)), ((3, 1, 1), (2, 1, 1)))
            d = [(1 - ratio_to_nr(sa, sb, α, 0.0)[1]) / α^2 for α in (0.02, 0.04)]
            println("  $(join(sa))->$(join(sb)) a=0: (1 - Γ/Γ_NR)/α² = $(d[1]) (α=0.02), $(d[2]) (α=0.04)")
            @test abs(d[1] / d[2] - 1) < 0.02
            @test 0 < d[1] < 20
        end
    end

    @testset "ℓ structure of 433->322: current quadrupole / mass octupole = 25/2744 in the NR limit" begin
        _, modes, P = ratio_to_nr((4, 3, 3), (3, 2, 2), 0.02, 0.0)
        @test abs(modes[3] / modes[2] / (25 / 2744) - 1) < 0.05
    end

    @testset "spin-orbit term of the |Δm| = 1 (current-type) channels: 1 - Γ/Γ_NR ≈ 3 a α / l_b" begin
        # the lowest channel reproduces the (2,1)-mode factor (1 - 3 q v) of the post-Newtonian waveform
        for (sa, sb, coef) in (((3, 2, 2), (2, 1, 1), 3.0), ((4, 3, 3), (3, 2, 2), 1.5), ((5, 4, 4), (4, 3, 3), 1.0))
            α = 0.02; a = 0.9
            d = 1 - ratio_to_nr(sa, sb, α, a)[1]
            println("  $(join(sa))->$(join(sb)) a=0.9 α=0.02: 1 - Γ/Γ_NR = $d  vs $(coef) a α = $(coef * a * α)")
            @test abs(d / (coef * a * α) - 1) < 0.06
        end
        # odd in a m: the retrograde cloud gets the opposite correction; mirror symmetry is exact
        r_pro = ratio_to_nr((3, 2, 2), (2, 1, 1), 0.02, 0.9)[1]
        r_ret = ratio_to_nr((3, 2, -2), (2, 1, -1), 0.02, 0.9)[1]
        r_mir = ratio_to_nr((3, 2, -2), (2, 1, -1), 0.02, -0.9)[1]
        @test abs((1 - r_pro) + (1 - r_ret)) < 0.1 * abs(1 - r_pro)
        @test r_mir ≈ r_pro rtol = 1e-8
        # mass-type channels have no O(a α) term
        @test abs(ratio_to_nr((4, 3, 3), (2, 1, 1), 0.02, 0.9)[1] - ratio_to_nr((4, 3, 3), (2, 1, 1), 0.02, 0.0)[1]) < 3e-4
    end

    @testset "annihilation 2p x 2p, Schwarzschild: α^14 with ~0.43 x the semi-analytic coefficient" begin
        # The exact rate is ~1e-9 of the bulk partial sums of its integral; it is converged in the grid
        # (ppr, ppw) and the outer cutoff (decay >= 45), identical in Float64 and BigFloat, and the kernel
        # and direct source forms agree. Yoshino & Kodama 2014 (Fig. 5, a = 0) find the same power law with
        # a constant offset below the (484+9π²)/23040 α^14 of Brito+ 2017 (~0.38 read off their figure).
        Cb = (484 + 9π^2) / 46080
        rat = Dict{Float64,Float64}()
        for α in (0.1, 0.07)
            ma, _ = kerr_mode(2, 1, 1, α, 0.0)
            Γ, P, modes, ω̃, m̃ = gw_kerr_rate(ma, ma, :ann; decay=60.0)
            rat[α] = Γ / (Cb * α^14)
            println("  211x211 a=0 α=$α: Γ/(Brito α^14) = $(rat[α])   ℓ=3/ℓ=2 = $(get(modes, 3, 0.0) / modes[2])")
            @test 0.3 < rat[α] < 0.5
            @test get(modes, 3, 0.0) < 1e-12 * modes[2]        # a = 0: odd ℓ̃ does not radiate (parity)
        end
        @test abs(rat[0.07] - rat[0.1]) < 0.05
        # independent of the cutoff once the tail is negligible
        ma, _ = kerr_mode(2, 1, 1, 0.1, 0.0)
        Γ45 = gw_kerr_rate(ma, ma, :ann; decay=45.0)[1]
        @test abs(Γ45 / (Cb * 0.1^14) / rat[0.1] - 1) < 1e-3
    end
end
