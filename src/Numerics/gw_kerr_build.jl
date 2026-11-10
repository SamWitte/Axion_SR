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
the previous eigenfrequency scaled by α/α_prev (the hydrogenic seed fails at
α ≳ 0.5). Points where the mode or the Teukolsky integration fails are written
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

# seed for kerr_mode at α from the converged eigenfrequency at α_prev
# Seed for kerr_mode at α from the converged eigenfrequency at α_prev: the hydrogenic expansion at α
# (km_omega_guess, to O(α^5) relative) plus the previous root's deviation from its own expansion, scaled
# as α^7 (the first term the expansion lacks, in units M = 1 where ω ≈ α). (Scaling ω_prev itself by
# α/α_prev keeps ω/μ fixed while the binding energy grows as α², and lands on a higher overtone.)
function gw_kerr_seed(ωprev, αprev, α, s::NTuple{3,Int}, a)
    ωprev === nothing && return nothing
    g0 = km_omega_guess(s..., αprev, a); g1 = km_omega_guess(s..., α, a)
    return g1 + (ωprev - g0) * (α / αprev)^7
end

"""Distance (in ω, M = 1) from level (n,l,m) to the nearest other overtone with the same (l, m), hydrogenic."""
function gw_kerr_overtone_gap(s::NTuple{3,Int}, α)
    n, l, _ = s
    up = 1 / (2n^2) - 1 / (2(n + 1)^2)
    dn = n - 1 > l ? 1 / (2(n - 1)^2) - 1 / (2n^2) : Inf
    return α^3 * min(up, dn)
end

"""
    gw_kerr_nearest_overtone(s, ω, α, a) -> n'

The overtone n' (same l, m; n' = l+1 ... n+4) whose hydrogenic expansion (km_omega_guess) is closest to Re ω.
"""
gw_kerr_nearest_overtone(s::NTuple{3,Int}, ω, α, a) =
    argmin(np -> abs(real(ω) - real(km_omega_guess(np, s[2], s[3], α, a))), (s[2] + 1):(s[1] + 4))

"""
    gw_kerr_mode_checked(s, α, a; seed=nothing, tol=0.5, mode_tol=1e-8) -> (KerrMode, status)

kerr_mode for level s = (n,l,m) at (α, a), from `seed` (continuation) or the hydrogenic guess, with the
root checked against overtone confusion. A seeded root must lie within tol × the gap to the neighbouring
overtones of the seed (a continuation seed is accurate); otherwise, and for an unseeded solve, the root is
accepted if its nearest hydrogenic overtone is n (at large α the expansion itself is off by a sizeable
fraction of the gap, but the overtones shift together). status = "ok", "res" (no converged root) or
"root" (the converged root belongs to another overtone).
"""
function gw_kerr_mode_checked(s::NTuple{3,Int}, α, a; seed=nothing, tol=0.5, mode_tol=1e-8)
    gap = gw_kerr_overtone_gap(s, α)
    last = nothing; lastres = Inf
    for ω0 in (seed === nothing ? (nothing,) : (seed, nothing))
        md, res = try
            kerr_mode(s..., α, a; ω0=ω0)
        catch
            (nothing, Inf)
        end
        md === nothing && continue
        last = md; lastres = abs(res)
        abs(res) < mode_tol || continue
        ω0 !== nothing && abs(real(md.ω) - real(ω0)) < tol * gap && return md, "ok"
        gw_kerr_nearest_overtone(s, md.ω, α, a) == s[1] && return md, "ok"
    end
    last === nothing && error("kerr_mode found no root for $(s) at α=$α a=$a")
    return last, lastres < mode_tol ? "root" : "res"
end

"""
    gw_kerr_table_point(kind, sa, sb, α, a; nr, ω0a, ω0b, T, lmax, ltol, decay, mode_tol) -> NamedTuple

One grid point. `nr` is the NR reference: the transition channel list (kind = :tr) or
(C, p) (kind = :ann). Returns the table row as a NamedTuple with the columns of
GW_KERR_COLS plus the two KerrModes (for seeding the next point).
"""
function gw_kerr_table_point(kind::Symbol, sa, sb, α, a; nr, ω0a=nothing, ω0b=nothing, T=Float64,
                             lmax::Int=12, ltol=1e-8, decay=60.0, mode_tol=1e-8, horizon::Bool=true,
                             chk_tol=1e-3, kerr_kwargs...)
    t0 = time()
    nan = NaN
    fail(reason, ma=nothing, mb=nothing) = (alpha=α, a=a, Gamma_kerr=nan, Gamma_nr=nan, ratio=nan, omega_tilde=nan,
                                          Re_omega_a=ma === nothing ? nan : real(ma.ω) / α,
                                          Re_omega_b=mb === nothing ? nan : real(mb.ω) / α,
                                          Im_omega_a=ma === nothing ? nan : imag(ma.ω) / α,
                                          Im_omega_b=mb === nothing ? nan : imag(mb.ω) / α,
                                          Gamma_H=nan, nl=0, secs=time() - t0, status="fail:" * reason, chk=nan, lfrac="-", ma=ma, mb=mb)
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
        # the direct and kernel forms of the source cancel differently: if their P_∞ disagree, the strongly
        # cancelling (small-α, high-l) annihilation integrals are at the precision floor of T -> not usable
        status = (isnan(chk) || abs(chk) < chk_tol) ? "ok" : "noisy"
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
            # ok and noisy rows are final (noisy is a precision verdict, recomputing gives the same); failures are retried
            if r.status in ("ok", "noisy")
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
    for a in as
        ωa_prev = ωb_prev = nothing; α_prev = nothing
        for α in sort(collect(αs))
            n += 1
            if (round(α, sigdigits=10), round(a, sigdigits=10)) in done
                ωa_prev, ωb_prev = done_ω[(round(α, sigdigits=10), round(a, sigdigits=10))]; α_prev = α
                continue
            end
            Tuse = α < bigfloat_below ? BigFloat : T
            r = gw_kerr_table_point(kind, sa, sb, α, a; nr=nr, T=Tuse, lmax=lmax, ltol=ltol, decay=decay, horizon=horizon,
                                    ω0a=gw_kerr_seed(ωa_prev, α_prev, α, sa, a), ω0b=gw_kerr_seed(ωb_prev, α_prev, α, sb, a),
                                    kerr_kwargs...)
            println(io, gw_kerr_row_string(r)); flush(io)
            if r.status in ("ok", "noisy")          # the modes converged: seed the next α from them
                ωa_prev = r.ma.ω; ωb_prev = r.mb.ω; α_prev = α
                r.status == "noisy" && (nfail += 1)
            else
                nfail += 1                           # keep seeding from the last good root (larger step)
            end
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
