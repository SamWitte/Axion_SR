"""
Forward model of the GW lines a BH + axion cloud emits today: draw birth M, a, age, f_Edd, distance and inclination
from BHPriors, evolve with solve_system to t = age, then gw_lines. Birth priors should be broad; measurements of
M and a today enter afterwards as weights on the realizations (plot_realizations.py).
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

Birth M [M_sun], spin, log10(age/yr), f_Edd, distance [kpc] and cos(iota), each a Distribution or a number.
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

# wigner_d (small-d matrix element) comes from Numerics/gw_kerr_tables.jl (included through super_rad.jl).

"""
    gw_pattern(L, m, cos_iota)

4 pi (dP/dOmega)/P of an (L, +-m) GW multipole (sphere average 1); L = m = 2: (5/16)(1 + 6c^2 + c^4).
"""
gw_pattern(L, m, c) = (2L + 1) / 2 * (wigner_d(L, m, 2, c)^2 + wigner_d(L, -m, 2, c)^2)

const TRANSITION_CHANNELS = Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, Vector{Tuple{Int,Symbol,Float64}}}()

"""
    inclination_factor(line, cos_iota, alpha, a; gw_model=:nonrel) -> h(iota) / h_avg
`line` from gw_lines; alpha, a today. :rel with a Kerr table: gw_rel_pattern; else non-relativistic multipoles
(annihilations: L = m only; transitions: all channels weighted by power, interference neglected).
"""
function inclination_factor(line, c, α, a; gw_model::Symbol=:nonrel)
    sa, sb = line.a, line.b
    if gw_model == :rel
        tabs = gw_kerr_load_tables()
        if tabs !== nothing
            pat = gw_rel_pattern(tabs, line.kind == :annihilation ? :ann : :tr, sa, sb, α, a, c)
            pat === nothing || return sqrt(pat)
        end
    end
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
    forward_model(p; mu, fa, Nmax=3, n_times=4000, window=0.9, rel_floor=1e-12, fdot_frac=0.01, gw_model, self_gravity)
Evolve p (mu [eV], fa [GeV]) to t = age; return (status, M_today, a_today, u_today, lines), lines from t >= window*age
(f [Hz], fdot over the last fdot_frac of age, P [erg/s], h_avg, h_iota), kept if h_iota >= rel_floor * max, sorted.
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
                   self_gravity=self_gravity, fa=fa, M0=p.M, levels=scattering_levels(Nmax))
    α = GNew * mb[end] * mu
    kp = max(1, searchsortedlast(res.t, (1 - fdot_frac) * res.t[end]))
    dt_s = (res.t[end] - res.t[kp]) * YEAR_IN_SECONDS
    lines = NamedTuple[]
    for L in res.lines
        L.P[end] > 0 || continue
        hι = L.h0[end] * inclination_factor(L, p.cos_iota, α, sp[end]; gw_model=gw_model)
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

Forward-model p and append its row to fs and its lines to fl; return the forward_model output (nothing on error).
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

Forward-model N draws into <outdir>/samples_<tag>.dat and lines_<tag>.dat (status error if one throws).
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
