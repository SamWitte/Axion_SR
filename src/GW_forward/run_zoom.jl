# Iterative zoom-in on the birth parameters, conditioned on a measurement of the
# BH mass and spin today.
#
#   julia run_zoom.jl --outdir out --tag zoom --M_obs 10 --a_obs 0.35 \
#                     [--rounds 3 --N_round 100] [--round1_from broad] [prior options]
#
# Parameters x = (M0, a0, log10 age, log10 f_Edd) have a uniform prior box P
# (same flags as run_forward.jl; a parameter with min = max is fixed). Distance
# and cos(iota) are always drawn from their priors (the measurement says nothing
# about them).
#
#   round 1      x ~ P  (or the realizations of an earlier run_forward.jl ensemble
#                drawn from the same prior: --round1_from <tag>)
#   round r > 1  x ~ uniform box B_r: weighted [0.5%, 99.5%] range of every free
#                parameter over all samples so far, padded by pad_frac × its width
#                (at least min_frac × the prior width), clipped to P
#
# All rounds are combined by importance sampling against the mixture of the
# proposals,
#   w(x) = L(x) p(x) / sum_r (N_r / N) q_r(x),
#   L = N(M_today; M_obs, M_obs_sigma) N(a_today; a_obs, a_obs_sigma)
#       × [a_obs_min < a_today < a_obs_max]   (each factor only if given),
# so the weights are a valid posterior even if a box cuts too tightly (round 1
# covers all of P). Realizations that did not reach t = age get w = 0.
#
# Output in <outdir>: samples_<tag>.dat and lines_<tag>.dat (as run_forward.jl,
# all rounds), rounds_<tag>.dat (box of every round) and weights_<tag>.dat
# (id, round, w / max w), rewritten after every round. Plot with
#   python plot_realizations.py <outdir> <tag>
# which picks up weights_<tag>.dat.
using ArgParse
include(joinpath(@__DIR__, "forward_model.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "--outdir";        arg_type = String;  required = true
    "--tag";           arg_type = String;  required = true
    "--round1_from";   arg_type = String;  default = "";   help = "tag of a run_forward.jl ensemble (same outdir, same prior) to use as round 1"
    "--rounds";        arg_type = Int;     default = 3
    "--N_round";       arg_type = Int;     default = 100;  help = "realizations per round"
    "--seed";          arg_type = Int;     default = 1
    "--M_obs";         arg_type = Float64; default = NaN;  help = "measured BH mass today [M_sun]"
    "--M_obs_sigma";   arg_type = Float64; default = 1.0
    "--a_obs";         arg_type = Float64; default = NaN;  help = "measured BH spin today"
    "--a_obs_sigma";   arg_type = Float64; default = 0.05
    "--a_obs_min";     arg_type = Float64; default = NaN;  help = "hard lower bound on the spin today"
    "--a_obs_max";     arg_type = Float64; default = NaN;  help = "hard upper bound on the spin today"
    "--pad_frac";      arg_type = Float64; default = 0.25
    "--min_frac";      arg_type = Float64; default = 0.02
    "--max_wall_sec";  arg_type = Float64; default = 600.0; help = "per-realization solver wall-clock limit (then counted as incomplete)"
    "--Nmax";          arg_type = Int;     default = 3
    "--m_a";           arg_type = Float64; default = NaN;  help = "axion mass [eV] (overrides --alpha)"
    "--alpha";         arg_type = Float64; default = 0.1;  help = "G M_ref m_a"
    "--M_ref";         arg_type = Float64; default = 10.0
    "--f_a";           arg_type = Float64; default = 1e18
    "--M0_min";        arg_type = Float64; default = 5.0
    "--M0_max";        arg_type = Float64; default = 20.0
    "--a0_min";        arg_type = Float64; default = 0.0
    "--a0_max";        arg_type = Float64; default = 0.998
    "--log10_age_min"; arg_type = Float64; default = 5.0
    "--log10_age_max"; arg_type = Float64; default = 8.0
    "--fedd_min";      arg_type = Float64; default = 1e-3; help = "f_Edd log-uniform in [fedd_min, fedd_max]; equal values fix it (0 = no accretion)"
    "--fedd_max";      arg_type = Float64; default = 0.1
    "--d_mean";        arg_type = Float64; default = 5.0
    "--d_sigma";       arg_type = Float64; default = 0.5
    "--cos_iota";      arg_type = Float64; default = NaN
    "--rel_floor";     arg_type = Float64; default = 1e-12
    "--n_times";       arg_type = Int;     default = 4000
end
args = parse_args(s)
all(isnan, (args["M_obs"], args["a_obs"], args["a_obs_min"], args["a_obs_max"])) &&
    error("give --M_obs and/or a spin constraint (--a_obs, --a_obs_min, --a_obs_max)")
haskey(ENV, "NL_MAX_WALL_SEC") || (ENV["NL_MAX_WALL_SEC"] = string(args["max_wall_sec"]))

# ---------------------------------------------------------------------------
# parameter space: x = (M0, a0, log10 age, log10 f_Edd), uniform prior box
# ---------------------------------------------------------------------------
fmin, fmax = args["fedd_min"], args["fedd_max"]
fmax > fmin && fmin <= 0 && error("log-uniform f_Edd needs fedd_min > 0")
names_x = ["M0", "a0", "log10_age", "log10_fedd"]
P_lo = [args["M0_min"], args["a0_min"], args["log10_age_min"], fmax > fmin ? log10(fmin) : 0.0]
P_hi = [args["M0_max"], min(args["a0_max"], maxSpin), args["log10_age_max"], fmax > fmin ? log10(fmax) : 0.0]
fedd_fixed = fmax > fmin ? NaN : fmin
free = findall(P_hi .> P_lo)
println("free parameters: ", names_x[free])

x_of(p) = [p.M, p.a, log10(p.age), isnan(fedd_fixed) ? log10(p.f_edd) : 0.0]
lvol(lo, hi) = sum(log.(hi[free] .- lo[free]))
inbox(x, lo, hi) = all(lo[free] .<= x[free] .<= hi[free])

d_prior = args["d_sigma"] > 0 ? truncated(Normal(args["d_mean"], args["d_sigma"]), 0.01, Inf) : args["d_mean"]
ci_prior = isnan(args["cos_iota"]) ? Uniform(-1.0, 1.0) : args["cos_iota"]
function draw_in_box(rng, lo, hi)
    x = [lo[k] + (hi[k] - lo[k]) * rand(rng) for k in 1:4]
    return (M = x[1], a = x[2], age = exp10(x[3]), f_edd = isnan(fedd_fixed) ? exp10(x[4]) : fedd_fixed,
            d_kpc = draw(rng, d_prior), cos_iota = draw(rng, ci_prior))
end

# ---------------------------------------------------------------------------
# weights
# ---------------------------------------------------------------------------
function loglike(M, a)
    l = 0.0
    isnan(args["M_obs"]) || (l -= 0.5 * ((M - args["M_obs"]) / args["M_obs_sigma"])^2)
    isnan(args["a_obs"]) || (l -= 0.5 * ((a - args["a_obs"]) / args["a_obs_sigma"])^2)
    isnan(args["a_obs_min"]) || a > args["a_obs_min"] || return -Inf
    isnan(args["a_obs_max"]) || a < args["a_obs_max"] || return -Inf
    return l
end

# samples: Vector of (id, round, x, ok, M_today, a_today); boxes: Vector of (lo, hi, N)
function weights(samples, boxes)
    Ntot = sum(b[3] for b in boxes)
    lp = -lvol(P_lo, P_hi)
    lw = map(samples) do sm
        sm.ok || return -Inf
        q = sum(b[3] / Ntot * (inbox(sm.x, b[1], b[2]) ? exp(-lvol(b[1], b[2])) : 0.0) for b in boxes)
        loglike(sm.M, sm.a) + lp - log(q)
    end
    all(isinf, lw) && return zeros(length(lw)), 0.0     # nothing compatible with the measurement
    w = exp.(lw .- maximum(lw))
    return w, sum(w)^2 / sum(w .^ 2)
end

function wquantile(v, w, q)
    i = sortperm(v)
    c = cumsum(w[i]) ./ sum(w)
    return v[i][clamp(searchsortedfirst(c, q), 1, length(v))]
end

function next_box(samples, w, lo_prev, hi_prev)
    lo, hi = copy(P_lo), copy(P_hi)
    for k in free
        v = [sm.x[k] for sm in samples]
        a, b = wquantile(v, w, 0.005), wquantile(v, w, 0.995)
        pad = max(args["pad_frac"] * (b - a), args["min_frac"] * (P_hi[k] - P_lo[k]))
        lo[k], hi[k] = max(P_lo[k], a - pad), min(P_hi[k], b + pad)
    end
    return lo, hi
end

# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
mu = isnan(args["m_a"]) ? args["alpha"] / (GNew * args["M_ref"]) : args["m_a"]
fa = args["f_a"]
outdir, tag = args["outdir"], args["tag"]
fs, fl = init_outputs(outdir, tag)
rng = MersenneTwister(args["seed"])
fwd = (mu=mu, fa=fa, Nmax=args["Nmax"], n_times=args["n_times"], rel_floor=args["rel_floor"])
println("m_a = $mu eV (alpha = $(GNew * args["M_ref"] * mu) at $(args["M_ref"]) M_sun), f_a = $fa")

samples = NamedTuple[]
boxes = Tuple{Vector{Float64}, Vector{Float64}, Int}[]
read_rows(fn) = (m = readdlm(fn, comments=true); ndims(m) == 1 ? reshape(m, 1, :) : m)
sample_row(id, r, p, status, M, a) = (id = id, round = r, x = x_of(p), ok = status == "ok", M = M, a = a)

function write_state(samples, boxes, w, ess)
    open(joinpath(outdir, "rounds_$(tag).dat"), "w") do io
        println(io, "# round N ", join(["$(n)_lo $(n)_hi" for n in names_x], " "))
        for (r, (lo, hi, N)) in enumerate(boxes)
            println(io, r, " ", N, " ", join(["$(lo[k]) $(hi[k])" for k in 1:4], " "))
        end
    end
    open(joinpath(outdir, "weights_$(tag).dat"), "w") do io
        println(io, "# weights for M_obs = $(args["M_obs"]) +- $(args["M_obs_sigma"]), a_obs = $(args["a_obs"]) +- $(args["a_obs_sigma"]), ",
                "a_range = $(args["a_obs_min"]) $(args["a_obs_max"]); ESS = $(ess)")
        println(io, "# id round w")
        for (sm, wi) in zip(samples, w)
            println(io, sm.id, " ", sm.round, " ", wi)
        end
    end
end

next_id = 1
r0 = 1
if !isempty(args["round1_from"])
    src_s = read_rows(joinpath(outdir, "samples_$(args["round1_from"]).dat"))
    cp_lines = read(joinpath(outdir, "lines_$(args["round1_from"]).dat"), String)
    open(fl, "a") do io; write(io, join(filter(l -> !startswith(l, "#"), split(cp_lines, "\n", keepempty=false)) .* "\n")); end
    open(fs, "a") do io
        for i in 1:size(src_s, 1)
            row = src_s[i, :]
            println(io, join(row, " "))
            p = (M = Float64(row[2]), a = Float64(row[3]), age = Float64(row[4]), f_edd = Float64(row[5]))
            x = x_of(p)
            inbox(x, P_lo, P_hi) || @warn "imported realization $(row[1]) lies outside the prior box" x
            push!(samples, sample_row(Int(row[1]), 1, p, String(row[8]), Float64(row[9]), Float64(row[10])))
            global next_id = max(next_id, Int(row[1]) + 1)
        end
    end
    push!(boxes, (copy(P_lo), copy(P_hi), size(src_s, 1)))
    w, ess = weights(samples, boxes)
    write_state(samples, boxes, w, ess)
    println("round 1: imported $(size(src_s, 1)) realizations from $(args["round1_from"]); ESS = $(round(ess, digits=1))")
    r0 = 2
end

lo, hi = copy(P_lo), copy(P_hi)
for r in r0:args["rounds"]
    local w, ess
    if r > 1
        w, ess = weights(samples, boxes)
        if ess < 2
            @warn "ESS = $(round(ess, digits=2)) < 2: keeping the previous box for round $r"
        else
            global lo, hi = next_box(samples, w, lo, hi)
        end
    end
    println("\n=== round $r: box ", join(["$(names_x[k]) [$(round(lo[k], sigdigits=4)), $(round(hi[k], sigdigits=4))]" for k in free], ", "))
    flush(stdout)
    push!(boxes, (copy(lo), copy(hi), 0))
    for n in 1:args["N_round"]
        p = draw_in_box(rng, lo, hi)
        out = run_and_append(p, next_id, fs, fl; fwd..., label="round $r, $n/$(args["N_round"])")
        push!(samples, out === nothing ? sample_row(next_id, r, p, "error", NaN, NaN) :
                                         sample_row(next_id, r, p, String(out.status), out.M_today, out.a_today))
        global next_id += 1
        boxes[end] = (lo, hi, boxes[end][3] + 1)
    end
    w, ess = weights(samples, boxes)
    write_state(samples, boxes, w, ess)
    println("round $r done: $(length(samples)) realizations, ESS = $(round(ess, digits=1))")
    flush(stdout)
end
