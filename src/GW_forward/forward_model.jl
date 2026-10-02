"""
    forward_model.jl

Forward model of the GW lines a single BH + axion cloud emits today, for BH
parameters drawn from a prior. Each realization:

  1. draws the birth mass M0, birth spin a0, age, Eddington ratio f_Edd,
     distance d and inclination cos(iota) (angle between the BH spin and the
     line of sight) from `BHPriors`;
  2. evolves the cloud from birth to t = age with solve_system (accretion at
     constant Mdot = f_Edd Mdot_Edd(M0), see Core/accretion.jl);
  3. computes every GW line at t = age with gw_lines (frequency including BH
     drift and cloud self-gravity / self-interaction shifts, power, strain).

The birth priors should be broad: SR and accretion move M and a a long way.
Measurements of M and a today enter afterwards as weights on the
realizations (plot_realizations.py --M_obs ... --a_obs ...), so one ensemble
can be reused for different measurements.

Orientation. gw_lines returns the orientation-averaged amplitude
h_avg = sqrt(<A_+^2 + A_x^2>). The amplitude seen at inclination iota is
h(iota) = h_avg sqrt(4 pi (dP/dOmega)(iota) / P). Each line has a definite
azimuthal number (m_a + m_b for annihilations, m_a - m_b for transitions), and
we take dP/dOmega from the spin-weighted harmonics of the leading multipoles,
    4 pi (dP/dOmega) / P = (2L+1)/2 [ d^L_{m,2}(iota)^2 + d^L_{-m,2}(iota)^2 ] :
  - annihilations: only the lowest multipole, L = m_a + m_b. This is an
    approximation: L = m dominates the relativistic 211 x 211 flux at small
    alpha (Yoshino & Kodama 2014), and higher L are neglected. For 211 x 211 it
    gives the familiar h_+ ∝ (1 + cos^2 iota)/2, h_x ∝ cos iota;
  - transitions: every allowed (L, mass/current) channel, weighted by its power
    at the current alpha and Delta omega. Interference between channels with
    different L is neglected. Each channel's pattern is exact.
"""

using Random
using Distributions
using DelimitedFiles
using Printf
using Suppressor
using SpecialFunctions: loggamma
@suppress include(joinpath(@__DIR__, "..", "super_rad.jl"))

"""
    BHPriors(; M, a, log10_age, f_edd, d_kpc, cos_iota)

Distributions for the BH birth mass [M_sun], birth spin, log10(age/yr),
Eddington ratio, distance [kpc] and cos(inclination). Any field can be a
Distributions.jl distribution or a fixed number.
"""
Base.@kwdef struct BHPriors
    M = Uniform(5.0, 20.0)
    a = Uniform(0.0, maxSpin)
    log10_age = Uniform(6.0, 8.0)
    f_edd = 0.0
    d_kpc = truncated(Normal(5.0, 0.5), 0.01, Inf)
    cos_iota = Uniform(-1.0, 1.0)      # isotropic orientation
end

draw(rng, x::Real) = Float64(x)
draw(rng, x::Distribution) = Float64(rand(rng, x))

"""One realization (M, a, age, f_edd, d_kpc, cos_iota) from the priors."""
function sample_system(rng, pr::BHPriors)
    return (M = draw(rng, pr.M), a = draw(rng, pr.a), age = exp10(draw(rng, pr.log10_age)),
            f_edd = draw(rng, pr.f_edd), d_kpc = draw(rng, pr.d_kpc), cos_iota = draw(rng, pr.cos_iota))
end

# ----------------------------------------------------------------------------
# Angular emission pattern
# ----------------------------------------------------------------------------

"""Wigner small-d matrix element d^j_{m'm}(beta), with cos(beta) = c."""
function wigner_d(j::Int, mp::Int, m::Int, c::Float64)
    (abs(m) > j || abs(mp) > j) && return 0.0
    ch, sh = sqrt((1 + c) / 2), sqrt(max(1 - c, 0.0) / 2)
    lf(n) = loggamma(n + 1.0)
    pre = 0.5 * (lf(j + m) + lf(j - m) + lf(j + mp) + lf(j - mp))
    s = 0.0
    for k in max(0, m - mp):min(j + m, j - mp)
        term = exp(pre - lf(j + m - k) - lf(k) - lf(j - k - mp) - lf(k - m + mp))
        s += (isodd(k - m + mp) ? -1 : 1) * term * ch^(2j - 2k + m - mp) * sh^(2k - m + mp)
    end
    return s
end

"""
    gw_pattern(L, m, cos_iota)

4 pi (dP/dOmega)/P of a single (L, m) GW multipole radiated at positive and
negative frequency (sphere average = 1). L = m = 2: (5/16)(1 + 6c^2 + c^4).
"""
gw_pattern(L, m, c) = (2L + 1) / 2 * (wigner_d(L, m, 2, c)^2 + wigner_d(L, -m, 2, c)^2)

const TRANSITION_CHANNELS = Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, Vector{Tuple{Int,Symbol,Float64}}}()

"""
    inclination_factor(line, cos_iota, alpha, a) -> h(iota) / h_avg

`line` is an entry of gw_lines(...).lines; alpha and a are today's values.
"""
function inclination_factor(line, c, α, a)
    sa, sb = line.a, line.b
    if line.kind == :annihilation
        m = sa[3] + sb[3]
        return sqrt(gw_pattern(max(2, m), m, c))
    end
    m = sa[3] - sb[3]
    chans = get!(() -> gw_transition_channels(sa, sb), TRANSITION_CHANNELS, (sa, sb))
    δ = abs(gw_omega(sa..., α, a) - gw_omega(sb..., α, a))
    num, den = 0.0, 0.0
    for (L, kind, K) in chans
        w = (kind == :mass ? K * α^(2 - 2L) : K * α^(4 - 2L)) * δ^(2L + 1)
        num += w * gw_pattern(L, m, c)
        den += w
    end
    return den > 0 ? sqrt(num / den) : 1.0
end

# ----------------------------------------------------------------------------
# One realization
# ----------------------------------------------------------------------------

"""
    forward_model(p; mu, fa, Nmax=3, n_times=4000, window=0.9, rel_floor=1e-12,
                  fdot_frac=0.01, gw_model=:nonrel, self_gravity=true)

Evolve realization `p` (from sample_system) for mu [eV] and fa [GeV] (same
convention as run_Nlevels.jl / gw_lines) and return its GW lines today.

Lines are computed from the part of the trajectory with t >= window × age, and
kept if their strain today is >= rel_floor × the brightest line today. Every
channel is included, also those the solver drops as dynamically irrelevant
(gw_min_rate_per_yr), so faint lines such as 322 x 322 appear.
fdot is the finite difference between t = age and the last saved time
<= (1 - fdot_frac) × age.

Returns (status, M_today, a_today, u_today, lines), with `lines` a vector of
NamedTuples (kind, a, b, f [Hz], fdot [Hz/s], df_cloud [Hz], P [erg/s], h_avg,
h_iota) sorted by h_iota. status is :ok, or :incomplete if the solver stopped
before t = age.
"""
function forward_model(p; mu, fa, Nmax=3, n_times=4000, window=0.9, rel_floor=1e-12,
                       fdot_frac=0.01, gw_model=:nonrel, self_gravity=true)
    t, st, modes, sp, mb = solve_system(mu, fa, p.a, p.M, p.age;
        return_all_info=true, n_times=n_times, Nmax=Nmax, f_edd=p.f_edd, gw_model=gw_model,
        impose_low_cut=1e-100, eq_threshold=1e-100, abstol=1e-30, non_rel=false,
        high_p=true, cheby=true, N_pts_interp=100, N_pts_interpL=100)
    states = permutedims(hcat(st...))            # (levels, times)
    status = isapprox(t[end], p.age; rtol=1e-6) ? :ok : :incomplete
    u_today = sum(states[:, end])

    nt = length(t)
    k0 = min(searchsortedfirst(t, window * t[end]), nt - 1)
    ks = k0:nt
    tw = t[ks]
    res = gw_lines(tw, states[:, ks], modes, sp[ks], mb[ks], mu, Nmax; d_kpc=p.d_kpc,
                   rel_floor=0.0, n_out=length(ks), gw_model=gw_model, min_rate_per_yr=0.0,
                   self_gravity=self_gravity, fa=fa, M0=p.M)
    α = GNew * mb[end] * mu
    kp = max(1, searchsortedlast(res.t, (1 - fdot_frac) * res.t[end]))
    dt_s = (res.t[end] - res.t[kp]) * YEAR_IN_SECONDS
    lines = NamedTuple[]
    for L in res.lines
        L.P[end] > 0 || continue
        hι = L.h0[end] * inclination_factor(L, p.cos_iota, α, sp[end])
        push!(lines, (kind = L.kind, a = L.a, b = L.b, f = L.f[end],
                      fdot = (L.f[end] - L.f[kp]) / dt_s, df_cloud = L.df_cloud[end],
                      P = L.P[end], h_avg = L.h0[end], h_iota = hι))
    end
    if !isempty(lines)
        hmax = maximum(l.h_iota for l in lines)
        filter!(l -> l.h_iota >= rel_floor * hmax, lines)
        sort!(lines, by=l -> -l.h_iota)
    end
    return (status = status, M_today = mb[end], a_today = sp[end], u_today = u_today, lines = lines)
end

# ----------------------------------------------------------------------------
# Ensemble
# ----------------------------------------------------------------------------

line_label(l) = l.kind == :annihilation ? "$(join(l.a))x$(join(l.b))" : "$(join(l.a))->$(join(l.b))"

const SAMPLES_HEADER = "# id M0_Msun a0 age_yr f_edd d_kpc cos_iota status M_today_Msun a_today u_cloud_today n_lines wall_s"
const LINES_HEADER = "# id kind na la ma nb lb mb f_Hz fdot_Hz_per_s df_cloud_Hz P_erg_per_s h_avg h_iota"

"""
    init_outputs(outdir, tag) -> (samples_file, lines_file)

Create (overwrite) <outdir>/samples_<tag>.dat and lines_<tag>.dat with headers.
"""
function init_outputs(outdir, tag)
    mkpath(outdir)
    fs, fl = joinpath(outdir, "samples_$(tag).dat"), joinpath(outdir, "lines_$(tag).dat")
    open(fs, "w") do io; println(io, SAMPLES_HEADER); end
    open(fl, "w") do io; println(io, LINES_HEADER); end
    return fs, fl
end

"""
    run_and_append(p, id, fs, fl; mu, fa, verbose=true, label="", kwargs...)

Forward-model realization `p` and append its row to `fs` and its lines to `fl`
(see run_ensemble). Returns the forward_model output, or nothing on error.
"""
function run_and_append(p, id, fs, fl; mu, fa, verbose=true, label="", kwargs...)
    t0 = time()
    out = try
        forward_model(p; mu=mu, fa=fa, kwargs...)
    catch err
        @warn "realization $id failed" p exception=(err, catch_backtrace())
        nothing
    end
    wall = time() - t0
    open(fs, "a") do io
        if out === nothing
            @printf(io, "%d %.6g %.6g %.6g %.6g %.6g %.6g error NaN NaN NaN 0 %.1f\n",
                    id, p.M, p.a, p.age, p.f_edd, p.d_kpc, p.cos_iota, wall)
        else
            @printf(io, "%d %.6g %.6g %.6g %.6g %.6g %.6g %s %.8g %.8g %.6g %d %.1f\n",
                    id, p.M, p.a, p.age, p.f_edd, p.d_kpc, p.cos_iota, out.status,
                    out.M_today, out.a_today, out.u_today, length(out.lines), wall)
        end
    end
    out === nothing && return nothing
    open(fl, "a") do io
        for l in out.lines
            @printf(io, "%d %s %d %d %d %d %d %d %.12g %.6g %.6g %.6g %.6g %.6g\n",
                    id, l.kind, l.a..., l.b..., l.f, l.fdot, l.df_cloud, l.P, l.h_avg, l.h_iota)
        end
    end
    if verbose
        @printf("[%s] M0=%.3g a0=%.3g age=%.3g yr f_edd=%.3g d=%.3g kpc cos_i=%.2f -> a=%.3g M=%.4g, %d lines (%.0f s)\n",
                label, p.M, p.a, p.age, p.f_edd, p.d_kpc, p.cos_iota, out.a_today, out.M_today,
                length(out.lines), wall)
        for l in out.lines[1:min(3, end)]
            @printf("      %-16s f = %.9g Hz  fdot = %.3g Hz/s  h = %.3g\n", line_label(l), l.f, l.fdot, l.h_iota)
        end
        flush(stdout)
    end
    return out
end

"""
    run_ensemble(priors, N; mu, fa, outdir, tag, seed=1, verbose=true, kwargs...)

Draw N realizations, forward-model each, and append the results to
<outdir>/samples_<tag>.dat (one row per realization) and
<outdir>/lines_<tag>.dat (one row per line), flushing after every realization.
Realizations that throw are recorded with status `error`. kwargs go to
forward_model.
"""
function run_ensemble(pr::BHPriors, N; mu, fa, outdir, tag, seed=1, verbose=true, kwargs...)
    fs, fl = init_outputs(outdir, tag)
    rng = MersenneTwister(seed)
    for id in 1:N
        run_and_append(sample_system(rng, pr), id, fs, fl; mu=mu, fa=fa, verbose=verbose,
                       label="$id/$N", kwargs...)
    end
    return fs, fl
end
