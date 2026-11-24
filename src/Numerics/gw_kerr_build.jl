"""
    gw_kerr_build.jl

Builds the relativistic GW correction tables read by gw_kerr_tables.jl: for one
channel (annihilation a x b or transition a -> b) and a grid of (α, a) points,
the Teukolsky rate of gw_kerr.jl, the non-relativistic rate of gw_rates.jl at the
same (Leaver) frequencies, and their ratio, appended row by row to the channel's
table file (so that a killed job resumes where it stopped).

Grid: by default the (α, a) grid of the superradiance tables (grid_sr_growth.jl:
60 log-spaced α from 0.03 to α_max(m) = 1.3 a_max m / (2 r_+(a_max)), a_max = 0.998,
capped at 2; 40 spins from 0.01 to 0.998), subsampled by `alpha_stride` and
`spin_stride` (the ratio is smooth), with the α range of the pair's more
restrictive level. Below `alpha_min` nothing is computed: the interpolant clamps
to the lowest computed α, which is right because the ratio tends to a constant
(transitions: 1 - c α²; annihilations: the exact leading coefficient over the
flat-space one).

Continuation: within one spin the α loop runs upward and seeds kerr_mode with
a seed extrapolated from the previous eigenfrequencies (gw_kerr_seed; the
hydrogenic expansion alone does not identify the overtone at large α). A mode is
accepted only with a located, validated eigenfrequency of the right overtone
(gw_kerr_mode_checked). Points where the mode or the Teukolsky integration fails are written
with status fail:<reason> and NaN values; the loader fills them from neighbours.

Usage: scripts/build_gw_kerr_tables.jl (one Slurm job per channel).
"""

using Dates

include(joinpath(@__DIR__, "teukolsky.jl"))
include(joinpath(@__DIR__, "kerr_modes.jl"))
include(joinpath(@__DIR__, "gw_kerr.jl"))

const GW_SR_AMAX = 0.998
const GW_SR_NALPHA = 60
const GW_SR_NSPIN = 40

"""α_max of the superradiance table of a level with azimuthal number m (grid_sr_growth.jl)."""
gw_sr_alpha_max(m::Int) = min(2.0, GW_SR_AMAX * m / (2 * (1 + sqrt(1 - GW_SR_AMAX^2))) * 1.3)

"""Extra spin nodes added to the subsampled superradiance grid: the stride-4 grid has 0.90 and 0.998 as its top
nodes, too coarse where the clouds saturate (a → a_c(α) near 1), so 0.95 is added (2026-11-10)."""
const GW_KERR_EXTRA_SPINS = [0.95]

"""The (α, a) grid of the superradiance tables for a level pair, subsampled by the strides, plus `extra_spins`."""
function gw_kerr_default_grid(sa, sb; alpha_stride::Int=2, spin_stride::Int=4, alpha_min=0.03, alpha_max=nothing,
                              extra_spins=GW_KERR_EXTRA_SPINS)
    amax = alpha_max === nothing ? min(gw_sr_alpha_max(abs(sa[3])), gw_sr_alpha_max(abs(sb[3]))) : alpha_max
    αfull = exp10.(LinRange(log10(0.03), log10(amax), GW_SR_NALPHA))
    afull = collect(LinRange(0.01, GW_SR_AMAX, GW_SR_NSPIN))
    αs = unique(vcat(αfull[1:alpha_stride:end], αfull[end]))
    as = sort(unique(vcat(afull[1:spin_stride:end], afull[end], Float64.(extra_spins))))
    return filter(α -> α >= alpha_min * (1 - 1e-12), αs), as
end

"""Non-relativistic annihilation (C, p) exactly as gw_build_cache uses it (with the l = 1 override)."""
function gw_nr_annihilation_coeff(sa, sb; literature_overrides::Bool=true)
    if literature_overrides
        lit = gw_literature_annihilation(sa, sb)
        lit === nothing || return lit
    end
    return gw_annihilation_leading(sa, sb)
end

"""Distance (in ω, M = 1) from level (n,l,m) to the nearest other overtone with the same (l, m), hydrogenic."""
gw_kerr_overtone_gap(s::NTuple{3,Int}, α) = km_level_gap(s[1], s[2], α)

"""Complex offset of a level's frequency from its hydrogenic expansion (km_omega_guess), in level spacings."""
gw_kerr_offset(ω, α, s::NTuple{3,Int}, a) = (ω - km_omega_guess(s..., α, a)) / gw_kerr_overtone_gap(s, α)

# power-law extrapolation of one component of the offset from (v0 at α0, v1 at α1) to α1 r; default power pdef
function gw_kerr_extrap(v0, v1, r01, r, pdef, pmax)
    (v0 === nothing || v0 * v1 <= 0 || abs(v0) < 1e-9) && return v1 * r^pdef
    return v1 * r^clamp(log(v1 / v0) / log(r01), 0.0, pmax)
end

"""
    gw_kerr_seed(ωprev, αprev, α, s, a; ωprev2=nothing, αprev2=nothing) -> seed or nothing

Seed for the eigenfrequency at α from the converged one at α_prev: the hydrogenic expansion at α plus the previous
root's complex offset from its own expansion, extrapolated in units of the level spacing as a power of α, real and
imaginary part separately. The first term the expansion lacks is O(α⁷) in ω, i.e. α⁴ level spacings, which is the
default power of the real part (the imaginary part is kept); with a second earlier root (ωprev2 at αprev2 < αprev)
the powers are measured from the two. (Scaling ω_prev itself by α/α_prev keeps ω/μ fixed while the binding energy
grows as α², and lands on a higher overtone.)
"""
function gw_kerr_seed(ωprev, αprev, α, s::NTuple{3,Int}, a; ωprev2=nothing, αprev2=nothing)
    ωprev === nothing && return nothing
    z1 = gw_kerr_offset(ωprev, αprev, s, a)
    two = ωprev2 !== nothing && αprev2 !== nothing && αprev2 < αprev
    z0 = two ? gw_kerr_offset(ωprev2, αprev2, s, a) : nothing
    r = α / αprev; r01 = two ? αprev / αprev2 : 1.0
    x = gw_kerr_extrap(two ? real(z0) : nothing, real(z1), r01, r, 4.0, 60.0)
    y = gw_kerr_extrap(two ? imag(z0) : nothing, imag(z1), r01, r, 0.0, 60.0)
    return km_omega_guess(s..., α, a) + complex(x, y) * gw_kerr_overtone_gap(s, α)
end

"""
    gw_kerr_track(s, a, hist, α; dz=0.15, dzstep=0.5, mode_tol=1e-8, hmin=1e-3) -> hist' or nothing

Follow level s at spin a from its last known eigenfrequency to α. `hist` = (α1, ω1, α0, ω0): the last two points
(α0, ω0 may be nothing). The level is advanced in steps of ln α, the whole distance first; a step is accepted when
the eigenfrequency is located (km_eigenfrequency error < mode_tol, bound) within dz level spacings of the seed
extrapolated from the two previous points, and the seed itself moves the offset by less than dzstep spacings;
otherwise the step is halved. So where a level drifts fast (large α: a whole spacing between two grid nodes, with
a growing width for levels that are not superradiant) it is followed in steps small enough that the neighbouring
overtone, one spacing away, cannot be mistaken for it. Returns the history ending at (α, ω), or nothing when the
step falls below hmin (the level can no longer be followed: it has become a broad resonance, or the root search
fails).
"""
function gw_kerr_track(s::NTuple{3,Int}, a, hist, α; dz=0.15, dzstep=0.5, mode_tol=1e-8, hmin=1e-3, maxsteps::Int=400)
    α1, ω1, α0, ω0 = hist
    α > α1 || return hist
    h = log(α / α1)
    for it in 1:maxsteps
        αn = h >= log(α / α1) * (1 - 1e-12) ? α : α1 * exp(h)
        seed = gw_kerr_seed(ω1, α1, αn, s, a; ωprev2=ω0, αprev2=α0)
        zs = gw_kerr_offset(seed, αn, s, a)
        good = false; ωn = seed
        if abs(zs - gw_kerr_offset(ω1, α1, s, a)) < dzstep
            e = try
                km_eigenfrequency(s..., αn, a; ω0=seed)
            catch
                nothing
            end
            if e !== nothing && e.converged && e.err < mode_tol && real(e.ω) < αn
                ωn = e.ω
                good = abs(gw_kerr_offset(ωn, αn, s, a) - zs) < dz
            end
        end
        if good
            α0, ω0 = α1, ω1
            α1, ω1 = αn, ωn
            αn == α && return (α1, ω1, α0, ω0)
            h = min(1.5h, log(α / α1))
        else
            h /= 2
            h < hmin && return nothing
        end
    end
    return nothing
end

"""History (α1, ω1, α0, ω0) of a level after a new point (α, ω)."""
gw_kerr_push(hist, α, ω) = hist === nothing ? (α, ComplexF64(ω), nothing, nothing) : (α, ComplexF64(ω), hist[1], hist[2])

"""
    gw_kerr_nearest_overtone(s, ω, α, a) -> n'

The overtone n' (same l, m; n' = l+1 ... n+4) whose hydrogenic expansion (km_omega_guess) is closest to Re ω.
"""
gw_kerr_nearest_overtone(s::NTuple{3,Int}, ω, α, a) =
    argmin(np -> abs(real(ω) - real(km_omega_guess(np, s[2], s[3], α, a))), (s[2] + 1):(s[1] + 4))

"""Largest α at which the hydrogenic expansion alone identifies the overtone of a level with orbital number l:
where the level's offset from the expansion, about 10 α⁴/l³ level spacings (measured for l = 3-10: 0.25, 0.11,
0.054, 0.040, 0.018, 0.0095 × α⁴), stays below 0.3 spacings. Above it the neighbouring overtone can be the root
nearest to the expansion (seen for l = 5, 6 at α = 1.7-2), so only a continuation seed identifies the level."""
gw_kerr_alpha_free(l::Int) = l == 0 ? 0.0 : (0.3 * l^3 / 10)^0.25

"""
    gw_kerr_mode_checked(s, α, a; seed=nothing, tol=0.35, mode_tol=1e-8, free_tol=0.3) -> (KerrMode, status)

kerr_mode for level s = (n,l,m) at (α, a), from `seed` (continuation) or the hydrogenic guess, accepted only if
(1) the eigenfrequency is located: its a-posteriori error estimate (kerr_mode's second return value, in units of
    the level spacing) is below mode_tol, and
(2) it is the right overtone: a seeded root lies within tol level spacings of its seed (a continuation seed is
    accurate; the neighbouring overtone is a whole spacing away); otherwise, and for an unseeded solve, the root
    is accepted if α ≤ gw_kerr_alpha_free(l), its nearest hydrogenic overtone is n and it lies within free_tol
    spacings of the expansion.
(3) for a radial overtone continued by the ODE (n - l - 1 ≥ 2, l ≥ 1, width below 5 % of the spacing) the radial
    function has n - l - 1 nodes (km_radial_nodes): an identification of the overtone that does not rest on the
    seed.
status = "ok", "res" (no located root), "root" (a root that cannot be assigned to this overtone) or "nodes" (a
located root with the wrong number of radial nodes).
(Until 2026-11-21 condition (1) was |continued-fraction residual| < 1e-8 in Float64, which rejected sound roots
at low spin, where the residual is steep, and accepted noise at high spin, where it is flat.)
"""
function gw_kerr_mode_checked(s::NTuple{3,Int}, α, a; seed=nothing, tol=0.35, mode_tol=1e-8, free_tol=0.3)
    gap = gw_kerr_overtone_gap(s, α)
    last = nothing; lastres = Inf; nodes = false
    for ω0 in (seed === nothing ? (nothing,) : (seed, nothing))
        md, res = try
            kerr_mode(s..., α, a; ω0=ω0)
        catch
            (nothing, Inf)
        end
        md === nothing && continue
        last = md; lastres = abs(res)
        abs(res) < mode_tol || continue
        # (3) radial overtones (ODE-continued modes): the radial function has n - l - 1 nodes
        if md.ode !== nothing && s[2] >= 1 && abs(imag(md.ω)) < 0.05 * gap && km_radial_nodes(md) != s[1] - s[2] - 1
            nodes = true
            continue
        end
        ω0 !== nothing && abs(md.ω - ω0) < tol * gap && return md, "ok"
        α <= gw_kerr_alpha_free(s[2]) && gw_kerr_nearest_overtone(s, md.ω, α, a) == s[1] &&
            abs(gw_kerr_offset(md.ω, α, s, a)) < free_tol && return md, "ok"
    end
    last === nothing && error("kerr_mode found no root for $(s) at α=$α a=$a")
    return last, lastres < mode_tol ? (nodes ? "nodes" : "root") : "res"
end

"""
    gw_kerr_table_point(kind, sa, sb, α, a; nr, ω0a, ω0b, T, lmax, ltol, decay, mode_tol) -> NamedTuple

One grid point. `nr` is the NR reference: the transition channel list (kind = :tr) or
(C, p) (kind = :ann). Returns the table row as a NamedTuple with the columns of
GW_KERR_COLS plus the two KerrModes (for seeding the next point). `lost` = :a or :b
writes the failure of a level that the builder could not follow to this α (gw_kerr_track).
status: "ok"; "noisy" (direct and kernel forms of the source disagree by more than chk_tol);
"broad" (a level wider than broad_tol of its level spacing: computed, but not a
quasi-stationary level; the loader does not use it); "fail:<reason>".
"""
function gw_kerr_table_point(kind::Symbol, sa, sb, α, a; nr, ω0a=nothing, ω0b=nothing, T=Float64,
                             lmax::Int=12, ltol=1e-8, decay=60.0, mode_tol=1e-8, horizon::Bool=true,
                             chk_tol=1e-3, broad_tol=0.5, lost::Symbol=:none, kerr_kwargs...)
    t0 = time()
    nan = NaN
    fail(reason, ma=nothing, mb=nothing) = (alpha=α, a=a, Gamma_kerr=nan, Gamma_nr=nan, ratio=nan, omega_tilde=nan,
                                          Re_omega_a=ma === nothing ? nan : real(ma.ω) / α,
                                          Re_omega_b=mb === nothing ? nan : real(mb.ω) / α,
                                          Im_omega_a=ma === nothing ? nan : imag(ma.ω) / α,
                                          Im_omega_b=mb === nothing ? nan : imag(mb.ω) / α,
                                          Gamma_H=nan, nl=0, secs=time() - t0, status="fail:" * reason, chk=nan, lfrac="-", ma=ma, mb=mb)
    lost == :none || return fail("mode_" * string(lost) * "_lost")     # the level could not be followed to this α
    ma = mb = nothing
    try
        ma, sta = gw_kerr_mode_checked(sa, α, a; seed=ω0a, mode_tol=mode_tol)
        sta == "ok" || return fail("mode_a_" * sta, ma, nothing)
        if sa == sb
            mb = ma
        else
            mb, stb = gw_kerr_mode_checked(sb, α, a; seed=ω0b, mode_tol=mode_tol)
            stb == "ok" || return fail("mode_b_" * stb, ma, mb)
        end
    catch e
        return fail("mode_" * string(typeof(e)), ma, mb)
    end
    ωa, ωb = real(ma.ω), real(mb.ω)
    ω̃ = kind == :tr ? ωa - ωb : ωa + ωb
    try
        hi, lo = (kind == :tr && ω̃ < 0) ? (mb, ma) : (ma, mb)
        Γ, P, modes, ω̃abs, m̃, ΓH, modesH, chk = gw_kerr_rate(hi, lo, kind; lmax=lmax, ltol=ltol, decay=decay, T=T, horizon=horizon, kerr_kwargs...)
        Γnr = kind == :tr ? gw_transition_rate(nr, α, ω̃abs / α) : nr[1] * α^nr[2]
        isfinite(Γ) && Γ > 0 && isfinite(Γnr) && Γnr > 0 || return fail("rate_nonfinite", ma, mb)
        # direct and kernel forms of the source disagree: precision floor reached, not usable
        status = (isnan(chk) || abs(chk) < chk_tol) ? "ok" : "noisy"
        # a level wider than broad_tol of the spacing to its neighbours is a resonance (value kept)
        (abs(imag(ma.ω)) > broad_tol * gw_kerr_overtone_gap(sa, α) || abs(imag(mb.ω)) > broad_tol * gw_kerr_overtone_gap(sb, α)) &&
            (status = "broad")
        return (alpha=α, a=a, Gamma_kerr=Γ, Gamma_nr=Γnr, ratio=Γ / Γnr, omega_tilde=ω̃ / α,
                Re_omega_a=ωa / α, Re_omega_b=ωb / α, Im_omega_a=imag(ma.ω) / α, Im_omega_b=imag(mb.ω) / α,
                Gamma_H=ΓH, nl=length(modes), secs=time() - t0, status=status, chk=chk, lfrac=gw_kerr_lfrac(modes, P), ma=ma, mb=mb)
    catch e
        return fail("rate_" * string(typeof(e)), ma, mb)
    end
end

"""Per-ℓ̃ fractions of the power to infinity, "ℓ:frac;..." for fractions above `floor` (for the line inclination patterns)."""
gw_kerr_lfrac(modes, P; floor=1e-6) = (ks = sort([l for l in keys(modes) if modes[l] >= floor * P]);
                                       isempty(ks) ? "-" : join(["$(l):$(round(modes[l] / P, sigdigits=4))" for l in ks], ";"))

gw_kerr_row_string(r) = join((string(r.alpha), string(r.a), string(r.Gamma_kerr), string(r.Gamma_nr), string(r.ratio),
                              string(r.omega_tilde), string(r.Re_omega_a), string(r.Re_omega_b), string(r.Im_omega_a),
                              string(r.Im_omega_b), string(r.Gamma_H), string(r.nl), string(round(r.secs, digits=2)), r.status,
                              string(r.chk), r.lfrac), " ")

"""
    gw_kerr_build_channel(kind, sa, sb; αs, as, path, resume=true, T=Float64, bigfloat_below=0.0,
                          lmax=nothing, ltol=1e-8, decay=60.0, horizon=true, verbose=true, kerr_kwargs...) -> path

Compute the table of one channel on the grid (αs × as) and write/append it to `path`
(default: gw_kerr_table_dir()/gw_kerr_table_name(kind, sa, sb)). With `resume`, points
already present with status ok are skipped. Points with α < bigfloat_below use BigFloat
for the Teukolsky integration. lmax defaults to |m̃| + 10. `horizon` adds the horizon
part of the rate (Gamma_H column; about doubles the cost).
"""
function gw_kerr_build_channel(kind::Symbol, sa, sb; αs=nothing, as=nothing, path=nothing, resume::Bool=true,
                               T=Float64, bigfloat_below=0.0, lmax=nothing, ltol=1e-8, decay=60.0, horizon::Bool=true,
                               alpha_stride::Int=2, spin_stride::Int=4, alpha_min=0.03, alpha_max=nothing,
                               extra_spins=GW_KERR_EXTRA_SPINS, verbose::Bool=true, kerr_kwargs...)
    kind in (:ann, :tr) || error("kind must be :ann or :tr")
    if αs === nothing || as === nothing
        αd, ad = gw_kerr_default_grid(sa, sb; alpha_stride=alpha_stride, spin_stride=spin_stride, extra_spins=extra_spins,
                                      alpha_min=alpha_min, alpha_max=alpha_max)
        αs === nothing && (αs = αd); as === nothing && (as = ad)
    end
    path === nothing && (path = joinpath(gw_kerr_table_dir(), gw_kerr_table_name(kind, sa, sb)))
    mkpath(dirname(path))
    m̃ = kind == :tr ? sa[3] - sb[3] : sa[3] + sb[3]
    lmax === nothing && (lmax = abs(m̃) + 10)
    nr = kind == :tr ? gw_transition_channels(sa, sb) : gw_nr_annihilation_coeff(sa, sb)
    kind == :tr && isempty(nr) && error("no non-relativistic transition channel for $(sa) -> $(sb)")
    done = Set{Tuple{Float64,Float64}}()
    done_ω = Dict{Tuple{Float64,Float64},Tuple{ComplexF64,ComplexF64}}()     # stored roots, to keep seeding on resume
    if resume && isfile(path)
        for r in gw_kerr_read_table(path)[4]
            # ok, noisy and broad rows are final (verdicts that recomputing would repeat); failures are retried
            if r.status in ("ok", "noisy", "broad")
                key = (round(r.alpha, sigdigits=10), round(r.a, sigdigits=10))
                push!(done, key)
                done_ω[key] = (r.alpha * complex(r.Re_omega_a, r.Im_omega_a), r.alpha * complex(r.Re_omega_b, r.Im_omega_b))
            end
        end
    end
    newfile = !(resume && isfile(path))
    io = open(path, newfile ? "w" : "a")
    if newfile
        println(io, "# gw_kerr table  kind=$(kind)  a=$(gw_kerr_level_str(sa))  b=$(gw_kerr_level_str(sb))  " *
                    "generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by gw_kerr_build.jl (M = 1, μ = α; rates in solver units)")
        println(io, "# nr reference: ", kind == :tr ? "gw_transition_channels at the Leaver splitting: " * string(nr) :
                                                      "C alpha^p with (C, p) = " * string(nr))
        println(io, "# options: lmax=$(lmax) ltol=$(ltol) decay=$(decay) T=$(T) bigfloat_below=$(bigfloat_below) horizon=$(horizon) $(kerr_kwargs)")
        println(io, "# columns: ", join(GW_KERR_COLS, " "))
    end
    flush(io)
    total = length(αs) * length(as); n = 0; nfail = 0; t0 = time()
    verbose && println("gw_kerr_build_channel $(kind) $(sa) $(sb): $(length(αs)) α x $(length(as)) a = $(total) points, $(length(done)) done; -> $(path)")
    same = sa == sb
    for a in as
        ha = hb = nothing                    # the last two eigenfrequencies of each level at this spin (gw_kerr_track)
        losta = lostb = false                # a level that could not be followed is not searched for at larger α
        for α in sort(collect(αs))
            n += 1
            key = (round(α, sigdigits=10), round(a, sigdigits=10))
            if key in done
                ha = gw_kerr_push(ha, α, done_ω[key][1]); hb = gw_kerr_push(hb, α, done_ω[key][2])
                continue
            end
            Tuse = α < bigfloat_below ? BigFloat : T
            # follow both levels from their last eigenfrequency to this α (in small steps where they drift fast)
            ta = (ha === nothing || losta) ? nothing : gw_kerr_track(sa, a, ha, α)
            ha !== nothing && ta === nothing && (losta = true)
            tb = same ? ta : ((hb === nothing || lostb) ? nothing : gw_kerr_track(sb, a, hb, α))
            same ? (lostb = losta) : (hb !== nothing && tb === nothing && (lostb = true))
            r = gw_kerr_table_point(kind, sa, sb, α, a; nr=nr, T=Tuse, lmax=lmax, ltol=ltol, decay=decay, horizon=horizon,
                                    ω0a=ta === nothing ? nothing : ta[2], ω0b=tb === nothing ? nothing : tb[2],
                                    lost=losta ? :a : (lostb ? :b : :none), kerr_kwargs...)
            println(io, gw_kerr_row_string(r)); flush(io)
            # next point starts from the followed eigenfrequencies
            a_ok = r.ma !== nothing && !startswith(r.status, "fail:mode_a_")      # level a was accepted at this point
            b_ok = r.mb !== nothing && !startswith(r.status, "fail:mode_")
            ha = ta !== nothing ? ta : (a_ok ? gw_kerr_push(ha, α, r.ma.ω) : ha)
            hb = tb !== nothing ? tb : (b_ok ? gw_kerr_push(hb, α, r.mb.ω) : hb)
            r.status == "ok" || (nfail += 1)
            if verbose
                println("  [$(n)/$(total)] α=$(round(α, sigdigits=5)) a=$(round(a, sigdigits=4)): ratio=$(round(r.ratio, sigdigits=6)) " *
                        "ω̃/μ=$(round(r.omega_tilde, sigdigits=6)) nℓ=$(r.nl) $(round(r.secs, digits=1)) s $(r.status)")
                flush(stdout)
            end
        end
    end
    close(io)
    verbose && println("done: $(total - length(done)) new points, $(nfail) failed or noisy, $(round((time() - t0) / 60, digits=1)) min")
    return path
end
