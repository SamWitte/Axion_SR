# Relativistic GW correction tables (src/Numerics/gw_kerr_tables.jl, gw_kerr_build.jl) and their
# hook in the evolution cache (gw_build_cache with gw_model=:rel).
#   julia test/test_gw_kerr_tables.jl        (a few minutes: five Teukolsky points)
using Test, Suppressor
include(joinpath(@__DIR__, "..", "src", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "src", "state_utils.jl"))
@suppress include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_rates.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_kerr_build.jl"))

@testset "GW Kerr tables" begin
    @testset "bilinear interpolant in (ln α, a) with edge clamping" begin
        αs = [0.05, 0.1, 0.2, 0.4]; as = [0.0, 0.5, 0.9]
        f(α, a) = 1.5 + 2log(α) - 0.7a + 0.3log(α) * a          # bilinear in (ln α, a): reproduced exactly
        it = gw_interp_from_points([(α, a, f(α, a)) for α in αs for a in as])
        @test it.x ≈ log.(αs) && it.y ≈ as
        for (α, a) in ((0.07, 0.3), (0.13, 0.85), (0.05, 0.0), (0.4, 0.9))
            @test gw_interp2(it, α, a) ≈ f(α, a) rtol = 1e-12
        end
        @test gw_interp2(it, 0.01, 0.3) ≈ gw_interp2(it, 0.05, 0.3)       # clamps outside the grid
        @test gw_interp2(it, 1.0, 0.95) ≈ gw_interp2(it, 0.4, 0.9)
        @test isfinite(gw_interp2(it, NaN, NaN)) && isfinite(gw_interp2(it, -1.0, 0.3))
        @test all(gw_interp_alpha_range(it) .≈ (0.05, 0.4))
        # NaN filling: along α first (lower α preferred), then empty columns from the nearest filled one
        z = [1.0 NaN 3.0; NaN NaN 4.0; 2.0 NaN 5.0]
        @test gw_fill_nan!(z)
        @test z[2, 1] == 1.0 && z[:, 2] == z[:, 1] && z[2, 3] == 4.0
        @test !gw_fill_nan!(fill(NaN, 2, 2))
        # default channel set
        ch = gw_kerr_default_channels(lmax=2, overtones=1)
        @test (:ann, (2, 1, 1), (2, 1, 1)) in ch && (:ann, (2, 1, 1), (3, 2, 2)) in ch && (:ann, (3, 1, 1), (3, 1, 1)) in ch
        @test (:tr, (3, 2, 2), (2, 1, 1)) in ch && (:tr, (3, 2, 2), (3, 1, 1)) in ch && !any(c -> c[1] == :tr && c[2] == (2, 1, 1), ch)
        @test gw_sr_alpha_max(1) ≈ 0.610131 rtol = 1e-4
        @test gw_sr_alpha_max(3) ≈ 1.830394 rtol = 1e-4
        αg, ag = gw_kerr_default_grid((3, 2, 2), (2, 1, 1))
        @test αg[1] ≈ 0.03 && αg[end] ≈ gw_sr_alpha_max(1) && length(αg) == 31 && length(ag) == 12 && ag[end] == 0.998 && 0.95 in ag
    end

    @testset "filling missing points: interpolate, approach the limit, extrapolate a little, never across a threshold" begin
        sa, sb = (3, 2, 2), (2, 1, 1)
        αs = exp.(range(log(0.03), log(0.3), length=8)); as = [0.9, 0.95, 0.998]  # both levels superradiant everywhere
        @test all(gw_kerr_sr_sign(sa, α, a) == 1 && gw_kerr_sr_sign(sb, α, a) == 1 for α in αs, a in as)
        lin(α, a) = 0.2 + 0.1 * log(α) - 0.05a
        # 1. interior gaps of 1-2 nodes: exact for ln ratio linear in ln α
        z = [lin(α, a) for α in αs, a in as]; zt = copy(z)
        z[3:4, 1] .= NaN; z[5, 2] = NaN
        h = gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
        @test h[3, 1] == h[4, 1] == h[5, 2] == 1 && z ≈ zt && count(==(0), h) == 21
        # a gap of 3 is not interpolated along α (spin interpolation, exact here, takes it)
        z = copy(zt); z[3:5, 2] .= NaN
        h = gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
        @test all(h[3:5, 2] .== 4) && z ≈ zt
        # 2. below the lowest computed α a transition tends to ratio 1: ln ratio ∝ α^p, p from the two lowest nodes
        for p in (1.0, 2.0)
            z = [-3.0 * α^p * (1 + a) for α in αs, a in as]; zt = copy(z); z[1:3, :] .= NaN
            h = gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
            @test all(h[1:3, :] .== 2) && z ≈ zt
        end
        z = [-3.0 * α^3 for α in αs, a in as]; z[1:2, :] .= NaN                  # p is clamped to [1, 2]
        gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
        @test z[2, 1] ≈ -3.0 * αs[3]^3 * (αs[2] / αs[3])^2
        # a transition ratio that has levelled off away from 1 is held, not pulled towards 1
        z = fill(log(0.54), 8, 3); z[1:3, :] .= NaN
        h = gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
        @test all(h[1:3, :] .== 2) && all(z[1:3, :] .≈ log(0.54))
        # an annihilation holds its lowest value: "limit" when flat, "held" otherwise
        z = fill(0.3, 8, 3); z[1:2, 1] .= NaN; z[3, 2] = NaN; z[1:2, 2] .= NaN; z[4, 2] = 0.6
        h = gw_kerr_fill_grid!(z, αs, as, :ann, sa, sb)
        @test all(h[1:2, 1] .== 2) && all(z[1:2, 1] .== 0.3)
        @test all(h[1:3, 2] .== 5)
        # 3. one node above the highest computed α is extrapolated; further nodes are not
        z = [lin(α, a) for α in αs, a in as]; zt = copy(z); z[7:8, :] .= NaN
        h = gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
        @test all(h[7, :] .== 3) && z[7, :] ≈ zt[7, :] && all(h[8, :] .== 5) && z[8, :] == z[7, :]
        # ... and not when the step is large (a resolved dip or a steep rise)
        z = [lin(α, a) for α in αs, a in as]; z[6, :] .+= 2.0; z[7:8, :] .= NaN
        h = gw_kerr_fill_grid!(z, αs, as, :tr, sa, sb)
        @test all(h[7, :] .== 5)
        # 4./5. never across a superradiance threshold: (2,1,1) at a = 0.3 stops growing at α ≈ 0.075
        as2 = [0.3]
        sg = [gw_kerr_sr_sign(sb, α, 0.3) for α in αs]
        k = findfirst(==(-1), sg)
        @test k !== nothing && k > 2 && all(sg[1:k-1] .== 1)
        z = [lin(α, 0.3) for α in αs, _ in as2]; z[k-1:k, 1] .= NaN
        h = gw_kerr_fill_grid!(z, αs, as2, :tr, sa, sb)
        @test all(h[k-1:k, 1] .== 5)
        # nothing finite: nothing filled
        @test all(gw_kerr_fill_grid!(fill(NaN, 3, 2), αs[1:3], as[1:2], :tr, sa, sb) .== -1)
        # row selection: noisy rows within GW_KERR_NOISY_ACCEPT and broad rows are used, failures not
        row(st, chk; ratio=0.9) = (alpha=0.1, a=0.5, ratio=ratio, status=st, chk=chk)
        @test gw_kerr_usable(row("ok", NaN)) && gw_kerr_usable(row("broad", NaN)) && gw_kerr_usable(row("noisy", 0.02))
        @test !gw_kerr_usable(row("noisy", 0.2)) && !gw_kerr_usable(row("fail:rate_ErrorException", NaN; ratio=NaN))
        @test !gw_kerr_usable(row("ok", NaN; ratio=-1.0))
    end

    @testset "levels: seeds, tracking between grid nodes, acceptance" begin
        s = (20, 6, 6); a = 0.8206666666666667
        # the seed is the hydrogenic expansion plus the previous offset scaled as α⁴ (one point) or as the power
        # measured from two points; the hydrogenic expansion itself is real
        @test imag(km_omega_guess(s..., 1.7, a)) == 0 && isfinite(real(km_omega_guess(3, 0, 0, 0.1, 0.5)))
        ω1 = km_omega_guess(s..., 0.5, a) - 0.01 * gw_kerr_overtone_gap(s, 0.5)
        @test real(gw_kerr_offset(gw_kerr_seed(ω1, 0.5, 0.6, s, a), 0.6, s, a)) ≈ -0.01 * 1.2^4 rtol = 1e-9
        ω0 = km_omega_guess(s..., 0.4, a) - 0.01 / 1.25^6 * gw_kerr_overtone_gap(s, 0.4)
        @test real(gw_kerr_offset(gw_kerr_seed(ω1, 0.5, 0.6, s, a; ωprev2=ω0, αprev2=0.4), 0.6, s, a)) ≈ -0.01 * 1.2^6 rtol = 1e-9
        @test gw_kerr_seed(nothing, 0.5, 0.6, s, a) === nothing
        @test gw_kerr_overtone_gap(s, 0.3) == km_level_gap(20, 6, 0.3)
        @test gw_kerr_alpha_free(6) ≈ 1.6 atol = 0.01
        @test gw_kerr_alpha_free(3) ≈ 0.95 atol = 0.01
        @test gw_kerr_alpha_free(0) == 0
        # (20,6,6) from α = 1.12 to the top of its grid: the offset grows to half a spacing and is followed
        e0 = km_eigenfrequency(s..., 1.1239, a)
        h = gw_kerr_track(s, a, gw_kerr_push(nothing, 1.1239, e0.ω), 1.8304)
        @test h !== nothing && h[1] == 1.8304
        @test real(gw_kerr_offset(h[2], 1.8304, s, a)) ≈ -0.4905 atol = 2e-3
        # the acceptance of a mode: located and the right overtone
        md, st = gw_kerr_mode_checked(s, 1.8304, a; seed=h[2])
        @test st == "ok" && abs(md.ω - h[2]) < 1e-6 * gw_kerr_overtone_gap(s, 1.8304)
        # without a seed above gw_kerr_alpha_free the overtone cannot be named: not accepted
        @test gw_kerr_mode_checked(s, 1.8304, a)[2] == "root"
        @test gw_kerr_mode_checked(s, 0.3, a)[2] == "ok"
        @test km_radial_nodes(md) == 13
        # a seed one spacing off converges to the neighbouring overtone, which lies within the seed tolerance and
        # is refused by its number of radial nodes (14 for (21,6,6))
        bad = h[2] + gw_kerr_overtone_gap(s, 1.8304)
        @test km_radial_nodes(kerr_mode(s..., 1.8304, a; ω0=bad)[1]) == 14
        @test gw_kerr_mode_checked(s, 1.8304, a; seed=bad)[2] != "ok"
        @test km_radial_nodes(kerr_mode(39, 6, 6, 0.0692, 0.7193)[1]) == 32
        @test km_radial_nodes(kerr_mode(24, 10, 10, 0.0813, 0.998)[1]) == 13
        @test km_radial_nodes(kerr_mode(7, 6, 6, 0.3, 0.9)[1]) == -1                 # pure series: not counted
        # (5,3,3) is not superradiant at this spin above α ~ 0.8 and broadens fast: it is followed while it is
        # a level (width a quarter of the spacing at α = 1.29) with the imaginary part of the seed extrapolated
        sb = (5, 3, 3)
        eb = km_eigenfrequency(sb..., 0.9777, a)
        hb = gw_kerr_track(sb, a, gw_kerr_push(nothing, 0.9777, eb.ω), 1.2919)
        @test hb !== nothing
        zb = gw_kerr_offset(hb[2], 1.2919, sb, a)
        @test real(zb) ≈ -1.0673 atol = 5e-3
        @test imag(zb) ≈ -0.2604 atol = 5e-3
    end

    @testset "radial range of the rate integral: two extended levels" begin
        # (21,10,10) -> (16,8,8) at α = 0.3: the integral was cut where the product of the two radial functions had
        # decayed by e^-6 (ratio 0.531); converged in the range (decay 120 and 200 agree): 0.9465
        α, a = 0.3, 0.01
        ma, _ = gw_kerr_mode_checked((21, 10, 10), α, a); mb, _ = gw_kerr_mode_checked((16, 8, 8), α, a)
        nr = gw_transition_channels((21, 10, 10), (16, 8, 8))
        Γnr = gw_transition_rate(nr, α, (real(ma.ω) - real(mb.ω)) / α)
        r60 = gw_kerr_rate(ma, mb, :tr; lmax=12, horizon=true)
        r120 = gw_kerr_rate(ma, mb, :tr; lmax=12, horizon=true, decay=120.0)
        @test r60[1] / Γnr ≈ 0.9465 atol = 2e-3
        @test r60[1] ≈ r120[1] rtol = 2e-3
        @test abs(r60[8]) < 1e-3                                  # direct and kernel forms of the source agree
        # the range: beyond the more compact level's turning point by decay/q of that level while the other is
        # still extended, then at the sum of the rates
        x = gw_kerr_rmax(ma, mb, 60.0)
        @test x > gw_kerr_turning(ma) > gw_kerr_turning(mb)
        @test gw_kerr_rmax(ma, ma, 60.0) ≈ gw_kerr_turning(ma) + 60 / (2abs(real(ma.L.q))) rtol = 1e-12
    end

    @testset "flat-space limit at the top spin nodes" begin
        # tr (20,6,6) -> (5,3,3) at α = 0.0915: the ratio is the same at every spin to 4e-4. (The stored table of
        # 10 Nov had 0.96194 at a = 0.82, 0.96171 at 0.95 and 4e-13 at 0.998: wrong eigenfrequencies at the top
        # spins, and at every spin a norm of (20,6,6) that the 400-node quadrature had 0.12 % too large.)
        sa, sb = (20, 6, 6), (5, 3, 3); nr = gw_transition_channels(sa, sb)
        rs = [gw_kerr_table_point(:tr, sa, sb, 0.0915, a; nr=nr, lmax=13) for a in (0.8207, 0.998)]
        @test all(r -> r.status == "ok", rs)
        @test rs[1].ratio ≈ 0.96306 atol = 2e-4
        @test rs[2].ratio ≈ rs[1].ratio atol = 4e-4
        @test all(r -> abs(r.Im_omega_a) < 1e-15, rs)
    end

    dir = mktempdir()
    ptr = joinpath(dir, gw_kerr_table_name(:tr, (3, 2, 2), (2, 1, 1)))
    pann = joinpath(dir, gw_kerr_table_name(:ann, (2, 1, 1), (2, 1, 1)))
    @testset "build a small table and read it back" begin
        αs = [0.05, 0.1]; as = [0.01, 0.5]
        gw_kerr_build_channel(:tr, (3, 2, 2), (2, 1, 1); αs=αs, as=as, path=ptr, verbose=false)
        kind, sa, sb, rows = gw_kerr_read_table(ptr)
        @test kind == :tr && sa == (3, 2, 2) && sb == (2, 1, 1) && length(rows) == 4 && all(r -> r.status == "ok", rows)
        r = first(r for r in rows if r.alpha == 0.05 && r.a == 0.01)
        @test 0.95 < r.ratio < 1.0                                   # 1 - c α² - 3 a α / l_b (test_gw_kerr.jl)
        @test r.ratio ≈ r.Gamma_kerr / r.Gamma_nr
        @test r.Gamma_nr ≈ gw_transition_rate(gw_transition_channels((3, 2, 2), (2, 1, 1)), 0.05, r.omega_tilde) rtol = 1e-12
        @test abs(r.Re_omega_a - gw_omega(3, 2, 2, 0.05, 0.01)) < 1e-6      # Leaver vs the O(α⁵) expansion
        @test r.omega_tilde ≈ r.Re_omega_a - r.Re_omega_b rtol = 1e-12
        r2 = first(r for r in rows if r.alpha == 0.05 && r.a == 0.5)
        @test r.ratio - r2.ratio ≈ 3 * 0.49 * 0.05 rtol = 0.15              # spin-orbit term 3 a α / l_b
        # horizon part: computed (default), tiny for a transition, superradiant (< 0) at a = 0.5; direct-form check
        @test all(x -> isfinite(x.Gamma_H) && abs(x.Gamma_H) < 1e-3 * x.Gamma_kerr && abs(x.chk) < 1e-8, rows)
        @test r2.Gamma_H < 0
        # resume: nothing recomputed, no duplicate rows
        gw_kerr_build_channel(:tr, (3, 2, 2), (2, 1, 1); αs=αs, as=as, path=ptr, verbose=false)
        @test length(gw_kerr_read_table(ptr)[4]) == 4
        # one annihilation point
        gw_kerr_build_channel(:ann, (2, 1, 1), (2, 1, 1); αs=[0.1], as=[0.5], path=pann, verbose=false)
        ra = gw_kerr_read_table(pann)[4][1]
        @test ra.status == "ok" && 0.3 < ra.ratio < 0.5                    # ~0.4 x the (484+9π²)/46080 α^14 coefficient
        @test ra.Gamma_H < 0 && abs(ra.Gamma_H) < 0.01 * ra.Gamma_kerr && abs(ra.chk) < 1e-5
        @test ra.Gamma_nr ≈ gw_literature_annihilation((2, 1, 1), (2, 1, 1))[1] * 0.1^14
        @test ra.omega_tilde ≈ 2 * ra.Re_omega_a
    end

    @testset "loader and the evolution hook (gw_model=:rel)" begin
        tabs = gw_kerr_load_tables(dir)
        @test length(tabs.ratio_tr) == 1 && length(tabs.ratio_ann) == 1
        @test Set(keys(tabs.omega)) == Set([(2, 1, 1), (3, 2, 2)])
        @test haskey(tabs.hfrac, (:tr, (2, 1, 1), (3, 2, 2))) && haskey(tabs.hfrac, (:ann, (2, 1, 1), (2, 1, 1)))
        rows = gw_kerr_read_table(ptr)[4]
        it = tabs.ratio_tr[gw_pair_key((3, 2, 2), (2, 1, 1))]
        for r in rows
            @test exp(gw_interp2(it, r.alpha, r.a)) ≈ r.ratio rtol = 1e-10           # grid points reproduced
            @test gw_omega(3, 2, 2, r.alpha, r.a) + gw_interp2(tabs.omega[(3, 2, 2)], r.alpha, r.a) ≈ r.Re_omega_a rtol = 1e-12
        end
        modes = [(s..., 0.0) for s in gw_state_list(4)]
        M = 10.0; α = 0.1; a = 0.5; mu = α / (GNew * M)
        cnr = gw_build_cache(4, modes, mu, M, a)
        crel = withenv("NL_GW_KERR_DIR" => dir) do
            @suppress gw_build_cache(4, modes, mu, M, a; gw_model=:rel)
        end
        @test crel.rel && !cnr.rel
        @test crel.ann_i == cnr.ann_i && crel.tr_i == cnr.tr_i && crel.ann_C == cnr.ann_C
        @test length(cnr.ann_itp) == length(cnr.ann_i) && all(==(0), cnr.ann_itp) && all(==(0), cnr.tr_itp)
        i211 = findfirst(==((2, 1, 1, 0.0)), modes); i322 = findfirst(==((3, 2, 2, 0.0)), modes); i433 = findfirst(==((4, 3, 3, 0.0)), modes)
        t = findfirst(t -> Set((crel.tr_i[t], crel.tr_j[t])) == Set((i211, i322)), eachindex(crel.tr_i))
        k = findfirst(k -> crel.ann_i[k] == i211 && crel.ann_j[k] == i211, eachindex(crel.ann_i))
        @test crel.tr_itp[t] != 0 && crel.ann_itp[k] != 0 && crel.tr_hitp[t] != 0 && crel.ann_hitp[k] != 0
        @test count(!=(0), crel.tr_itp) == 1 && count(!=(0), crel.ann_itp) == 1 && count(!=(0), crel.ω_itp) == 2
        @test gw_has_horizon(crel) && !gw_has_horizon(cnr)
        # at a grid point the hooked rates equal those of the tables, with the Leaver splitting: the levels
        # move at Γ_∞ + Γ_H, the line power uses Γ_∞
        r = first(r for r in rows if r.alpha == α && r.a == a)
        gw_set_omegas_cache!(crel, α, a)
        δ = gw_pair_delta(crel, i322, i211)
        @test δ ≈ r.omega_tilde rtol = 1e-12
        @test gw_tr_rate_inf(crel, t, α, a, abs(δ)) ≈ r.Gamma_kerr rtol = 1e-10
        @test gw_tr_rate(crel, t, α, a, abs(δ)) ≈ r.Gamma_kerr + r.Gamma_H rtol = 1e-10
        ra = gw_kerr_read_table(pann)[4][1]
        @test gw_ann_rate_inf(crel, k, α, a) ≈ ra.Gamma_kerr rtol = 1e-10
        @test gw_ann_rate(crel, k, α, a) ≈ ra.Gamma_kerr + ra.Gamma_H rtol = 1e-10
        # without the horizon tables (NL_GW_KERR_HORIZON=0) only Γ_∞ remains
        crel0 = withenv("NL_GW_KERR_DIR" => dir, "NL_GW_KERR_HORIZON" => "0") do
            @suppress gw_build_cache(4, modes, mu, M, a; gw_model=:rel)
        end
        @test !gw_has_horizon(crel0)
        @test gw_ann_rate(crel0, k, α, a) ≈ ra.Gamma_kerr rtol = 1e-10
        # BH rows: the hole absorbs ω̃ Γ_H u u (energy, units of μ) and m̃ Γ_H u u (angular momentum)
        uu = [1e-3 * (1 + 0.1i) for i in 1:length(modes)]
        nM = length(modes) + 2; spinI = nM - 1; massI = nM
        Xb = vcat(uu, [a, M])
        dub = gw_horizon_bh!(zeros(nM), Xb, crel, α, a, 0.0, spinI, massI)
        eM = 2ra.Re_omega_a * ra.Gamma_H * uu[i211]^2 + r.omega_tilde * r.Gamma_H * uu[i211] * uu[i322]
        eJ = 2 * ra.Gamma_H * uu[i211]^2 + 1 * r.Gamma_H * uu[i211] * uu[i322]
        @test dub[massI] ≈ eM rtol = 1e-8
        @test dub[spinI] ≈ eJ rtol = 1e-8
        @test all(==(0), dub[1:length(modes)])
        @test all(==(0), gw_horizon_bh!(zeros(nM), Xb, cnr, α, a, 0.0, spinI, massI))
        # channels and levels without a table are untouched (433 has no Leaver table: gw_omega for both levels of its pairs)
        gw_set_omegas_cache!(cnr, α, a)
        @test gw_pair_delta(crel, i433, i211) == gw_pair_delta(cnr, i433, i211)
        @test gw_pair_delta(crel, i433, i322) == gw_pair_delta(cnr, i433, i322)
        t2 = findfirst(t -> Set((crel.tr_i[t], crel.tr_j[t])) == Set((i433, i211)), eachindex(crel.tr_i))
        @test gw_tr_rate(crel, t2, α, a, 0.01) == gw_tr_rate(cnr, t2, α, a, 0.01)
        # gw_rhs!: the two caches differ only through the two tabulated channels
        u = [1e-3 * (1 + 0.1i) for i in 1:length(modes)]
        du_nr = gw_rhs!(zeros(length(modes)), u, cnr, α, a)
        du_rel = gw_rhs!(zeros(length(modes)), u, crel, α, a)
        ann, tr = gw_channel_rates(crel, α, a)
        @test sum(du_rel) ≈ -2 * sum(r_ * u[i] * u[j] for (i, j, r_) in ann) rtol = 1e-9     # transitions cancel in the sum
        @test du_rel[i433] ≈ du_nr[i433] rtol = 1e-12
        @test !isapprox(du_rel[i211], du_nr[i211]; rtol=1e-6)
        @test any(x -> x[1] == i322 && x[2] == i211 && x[3] ≈ r.Gamma_kerr + r.Gamma_H, tr)
        # off-grid evaluation is finite and clamps to the table edge
        @test gw_rel_factor(crel, crel.tr_itp[t], 1e-3, 0.999) ≈ exp(gw_interp2(it, 0.05, 0.5))
        # line frequencies: Leaver values where tabulated, consistent pairs for transitions
        ω = zeros(length(modes))
        gw_set_omegas!(ω, crel, M, a, mu)
        @test ω[i322] ≈ r.Re_omega_a && ω[i211] ≈ r.Re_omega_b && ω[i433] == gw_omega(4, 3, 3, α, a)
        @test gw_pair_omegas(crel, i433, i322, false) == (gw_omega(4, 3, 3, α, a), gw_omega(3, 2, 2, α, a))
        @test gw_pair_omegas(crel, i433, i322, true) == (gw_omega(4, 3, 3, α, a), r.Re_omega_a)
    end

    @testset "horizon pumping: no floor on the horizon fraction" begin
        # synthetic table for the same-n pair 322 <-> 311 with Γ_H = -1000 Γ_∞ (the real tables have -10² to -10⁶)
        d2 = mktempdir()
        open(joinpath(d2, gw_kerr_table_name(:tr, (3, 2, 2), (3, 1, 1))), "w") do io
            println(io, "# gw_kerr table  kind=tr  a=3-2-2  b=3-1-1  synthetic (test)")
            for α in (0.1, 0.3), a in (0.5, 0.9)
                wa, wb = gw_omega(3, 2, 2, α, a), gw_omega(3, 1, 1, α, a)
                println(io, join((α, a, 1.0, 1.25, 0.8, wa - wb, wa, wb, 1e-12, 1e-10, -1000.0, 3, 1.0, "ok", 0.0, "2:1.0"), " "))
            end
        end
        α, a, M = 0.2, 0.7, 10.0; mu = α / (GNew * M)
        modes = [(3, 2, 2, 0.0), (3, 1, 1, 0.0)]
        build() = withenv("NL_GW_KERR_DIR" => d2) do
            @suppress gw_build_cache(3, modes, mu, M, a; gw_model=:rel)
        end
        old = GW_HFRAC_FLOOR[]
        try
            @test old == -Inf                                  # default: no floor
            c = build()
            @test gw_kerr_reversible_pairs(gw_kerr_load_tables(d2)) == Set([gw_pair_key((3, 2, 2), (3, 1, 1))])
            t = findfirst(t -> Set((c.tr_i[t], c.tr_j[t])) == Set((1, 2)), eachindex(c.tr_i))
            @test t !== nothing && c.tr_hitp[t] != 0
            gw_set_omegas_cache!(c, α, a)
            δ = gw_pair_delta(c, 1, 2)
            @test δ > 0                                        # 322 is the upper level
            Γinf = gw_tr_rate_inf(c, t, α, a, δ)
            @test gw_hfrac(c, c.tr_hitp[t], α, a) ≈ -1000 rtol = 1e-8
            @test gw_tr_rate(c, t, α, a, δ) ≈ -999 * Γinf rtol = 1e-8      # negative net rate: 311 -> 322
            u = [0.02, 0.03, 0.0, 0.0]                         # rows 3, 4 stand for the spin and mass rows
            r = Γinf * u[1] * u[2]
            du1 = zeros(4); gw_rhs!(du1, u, c, α, a); gw_horizon_bh!(du1, u, c, α, a, 0.0, 3, 4)
            GW_HFRAC_FLOOR[] = -1.0                            # the earlier behaviour: the transition is switched off
            @test gw_tr_rate(c, t, α, a, δ) == 0
            du0 = zeros(4); gw_rhs!(du0, u, c, α, a); gw_horizon_bh!(du0, u, c, α, a, 0.0, 3, 4)
            # (annihilations are common to both and cancel in the differences)
            @test du1[1] - du0[1] ≈ 999r rtol = 1e-6           # the upper, higher-m level gains
            @test du1[2] - du0[2] ≈ -999r rtol = 1e-6          # the lower one loses the same number
            @test du1[4] ≈ -1000 * δ * r rtol = 1e-8           # the hole pays ω̃ per quantum for tide + radiation ...
            @test du1[3] ≈ -1000 * r rtol = 1e-8               # ... and m̃ = 1 unit of angular momentum
            @test du0[4] ≈ -δ * r rtol = 1e-8                  # with the floor it pays only for the radiation
            # energy: each transferred quantum gains δ, so the cloud's gain plus the hole's loss equals the flux to
            # infinity (written with δ: ω_a x - ω_b x would cancel to rounding against the 311 x 311 annihilation)
            @test δ * (du1[1] - du0[1]) + du1[4] ≈ -δ * r rtol = 1e-6
        finally
            GW_HFRAC_FLOOR[] = old
        end
    end

    @testset "spheroidal harmonics and the line pattern" begin
        # spectral ₋₂S_lm against the SpinWeightedSpheroidalHarmonics package (teuk_angular), eigenvalue at c = 0
        xq, wq = gauss(80)
        for (l, m, c) in ((2, 2, 0.0), (2, 2, 0.5), (3, 2, 0.8), (4, 1, 1.2), (2, -1, 0.6), (3, 0, 0.9), (5, 4, 1.5))
            A, lmin, b = gw_swsh(-2, l, m, c)
            c == 0 && @test A ≈ (l + 2) * (l - 1)
            @test abs(sum(b .^ 2) - 1) < 1e-12
            _, Spkg = teuk_angular(-2, l, m, c)
            θs = [0.3, 1.0, 1.7, 2.4, 2.9]
            σ = sign(gw_swsh_eval(-2, l, m, c, cos(θs[1])) / Spkg(θs[1]))
            for θ in θs
                @test gw_swsh_eval(-2, l, m, c, cos(θ)) ≈ σ * Spkg(θ) atol = 1e-7 * (1 + abs(Spkg(θ)))
            end
            @test 2π * sum(wq .* [gw_swsh_eval(-2, l, m, c, x)^2 for x in xq]) ≈ 1 rtol = 1e-10   # unit norm
        end
        # c = 0: the pattern of one ℓ̃ reduces to the non-relativistic gw_pattern (Wigner d's)
        pat(L, m, x) = 2π * (gw_swsh_eval(-2, L, m, 0.0, x)^2 + gw_swsh_eval(-2, L, m, 0.0, -x)^2)
        gwp(L, m, x) = (2L + 1) / 2 * (wigner_d(L, m, 2, x)^2 + wigner_d(L, -m, 2, x)^2)
        for (L, m) in ((2, 2), (3, 1), (4, 2), (2, 0)), x in (-0.9, -0.3, 0.2, 0.75, 1.0)
            @test pat(L, m, x) ≈ gwp(L, m, x) atol = 1e-10
        end
        # tabulated fractions: the small tables above carry lfrac; the pattern averages to 1 over the sphere
        tabs = gw_kerr_load_tables(dir)
        @test haskey(tabs.lfrac, (:tr, gw_pair_key((3, 2, 2), (2, 1, 1))...)) && haskey(tabs.lfrac, (:ann, (2, 1, 1), (2, 1, 1)))
        @test haskey(tabs.lfrac[(:tr, gw_pair_key((3, 2, 2), (2, 1, 1))...)], 2)
        for (kind, sa, sb, α, a) in ((:tr, (3, 2, 2), (2, 1, 1), 0.07, 0.3), (:ann, (2, 1, 1), (2, 1, 1), 0.1, 0.5))
            p = [gw_rel_pattern(tabs, kind, sa, sb, α, a, x) for x in xq]
            @test all(>=(0), p)
            @test sum(wq .* p) / 2 ≈ 1 rtol = 1e-8
        end
        @test gw_rel_pattern(tabs, :tr, (4, 3, 3), (2, 1, 1), 0.1, 0.5, 0.3) === nothing        # no table
        # 211 x 211 at α = 0.1, a = 0.5: on-axis pattern near the quadrupole value 5/2 (the spheroidal ₋₂S_22 at
        # c = aω̃ ≈ 0.1 is slightly more concentrated on the axis; the ℓ̃ = 3, 4 parts are < 1 % here)
        p1 = gw_rel_pattern(tabs, :ann, (2, 1, 1), (2, 1, 1), 0.1, 0.5, 1.0)
        @test 2.4 < p1 < 2.8
        @test gw_rel_pattern(tabs, :ann, (2, 1, 1), (2, 1, 1), 0.1, 0.5, 0.0) < 0.5      # equatorial: 5/16 for pure ℓ̃ = 2
    end
end
