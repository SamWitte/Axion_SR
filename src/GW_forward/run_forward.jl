# Random realizations of the GW lines from one BH with uncertain parameters.
#
#   julia run_forward.jl --N 50 --alpha 0.1 --f_a 1e18 --outdir out [prior options]
#
# Priors (all on birth values except distance/orientation; keep M0, a0 broad,
# the measured M and a today are imposed afterwards as weights in
# plot_realizations.py):
#   M0  ~ Uniform(M0_min, M0_max)                              [M_sun]
#   a0  ~ Uniform(a0_min, a0_max), a0_max <= 0.998
#   age ~ log-uniform in [10^log10_age_min, 10^log10_age_max]  [yr]
#   f_Edd: fixed fedd_min if fedd_max <= fedd_min, else log-uniform in [fedd_min, fedd_max]
#   d   ~ Normal(d_mean, d_sigma) truncated to d > 0.01        [kpc]
#   cos(iota) ~ Uniform(-1, 1), or fixed with --cos_iota
# Equal min and max (or a sigma of 0) fix that parameter. The axion mass is
# --m_a [eV], or set by --alpha = G M_ref m_a.
#
# Output (see forward_model.jl): <outdir>/samples_<tag>.dat and lines_<tag>.dat;
# plot with plot_realizations.py.
using ArgParse
include(joinpath(@__DIR__, "forward_model.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "--N";             arg_type = Int;     default = 20;   help = "number of realizations"
    "--seed";          arg_type = Int;     default = 1
    "--outdir";        arg_type = String;  required = true
    "--tag";           arg_type = String;  default = "";   help = "output file tag (default: built from alpha, f_a, seed)"
    "--Nmax";          arg_type = Int;     default = 3
    "--m_a";           arg_type = Float64; default = NaN;  help = "axion mass [eV] (overrides --alpha)"
    "--alpha";         arg_type = Float64; default = 0.1;  help = "G M_ref m_a"
    "--M_ref";         arg_type = Float64; default = 10.0; help = "reference mass for --alpha [M_sun]"
    "--f_a";           arg_type = Float64; default = 1e18; help = "axion decay constant (run_Nlevels.jl convention)"
    "--M0_min";        arg_type = Float64; default = 5.0
    "--M0_max";        arg_type = Float64; default = 20.0
    "--a0_min";        arg_type = Float64; default = 0.0
    "--a0_max";        arg_type = Float64; default = 0.998
    "--log10_age_min"; arg_type = Float64; default = 6.0
    "--log10_age_max"; arg_type = Float64; default = 8.0
    "--fedd_min";      arg_type = Float64; default = 0.0
    "--fedd_max";      arg_type = Float64; default = 0.0
    "--d_mean";        arg_type = Float64; default = 5.0
    "--d_sigma";       arg_type = Float64; default = 0.5
    "--cos_iota";      arg_type = Float64; default = NaN;  help = "fixed cos(inclination); default isotropic"
    "--rel_floor";     arg_type = Float64; default = 1e-12; help = "keep lines with h >= this × brightest today"
    "--n_times";       arg_type = Int;     default = 4000
    "--gw_model";      arg_type = String;  default = "nonrel"
    "--no_cloud_shifts"; action = :store_true
end
args = parse_args(s)

tnormal(m, σ, lo, hi) = σ > 0 ? truncated(Normal(m, σ), lo, hi) : clamp(m, lo, hi)
unif(lo, hi) = hi > lo ? Uniform(lo, hi) : lo
fmin, fmax = args["fedd_min"], args["fedd_max"]
f_edd = fmax > fmin ? (fmin > 0 ? LogUniform(fmin, fmax) : error("log-uniform f_Edd needs fedd_min > 0")) : fmin
lo, hi = args["log10_age_min"], args["log10_age_max"]
priors = BHPriors(
    M = unif(args["M0_min"], args["M0_max"]),
    a = unif(args["a0_min"], min(args["a0_max"], maxSpin)),
    log10_age = hi > lo ? Uniform(lo, hi) : lo,
    f_edd = f_edd,
    d_kpc = tnormal(args["d_mean"], args["d_sigma"], 0.01, Inf),
    cos_iota = isnan(args["cos_iota"]) ? Uniform(-1.0, 1.0) : args["cos_iota"],
)
mu = isnan(args["m_a"]) ? args["alpha"] / (GNew * args["M_ref"]) : args["m_a"]
fa = args["f_a"]
tag = isempty(args["tag"]) ? "ma_$(round(mu, sigdigits=4))_fa_$(fa)_seed_$(args["seed"])" : args["tag"]

println("m_a = $mu eV (alpha = $(GNew * args["M_ref"] * mu) at $(args["M_ref"]) M_sun), f_a = $fa, Nmax = $(args["Nmax"])")
println("priors: ", priors)
fs, fl = run_ensemble(priors, args["N"]; mu=mu, fa=fa, outdir=args["outdir"], tag=tag, seed=args["seed"],
                      Nmax=args["Nmax"], n_times=args["n_times"], rel_floor=args["rel_floor"],
                      gw_model=Symbol(args["gw_model"]), self_gravity=!args["no_cloud_shifts"])
println("wrote $fs\n      $fl")
