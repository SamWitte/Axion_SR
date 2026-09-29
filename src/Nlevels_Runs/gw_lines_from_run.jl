# GW lines (frequency, power, strain vs time) from a finished run_Nlevels.jl output.
#
#   julia gw_lines_from_run.jl --outdir <dir> --fname FullRel_fa_..._Nmax_5.dat [--d_kpc 1] [--rel_floor 1e-8]
#
# Reads Time_/States_/Modes_/Spin_/MassBH_<fname> (pruned States/Modes are fine) and writes,
# next to them (rows = lines, sorted by peak power; columns = the GWTime_ grid, n_out log-spaced times):
#   GWTime_<fname>    t [yr]
#   GWLines_<fname>   kind na la ma nb lb mb peak_P[erg/s] peak_h0 f_start[Hz] f_end[Hz]
#   GWFreq_<fname>    f(t) [Hz], including BH drift and cloud self-gravity + self-interaction shifts
#   GWFreqCloud_<fname>  the cloud-shift part of f(t) [Hz]
#   GWPower_<fname>   P(t) [erg/s]
#   GWStrain_<fname>  h0(t) at d_kpc (orientation-averaged sqrt(A+^2 + Ax^2))
# See gw_lines in ../Numerics/gw_rates.jl for definitions and current limitations.
using DelimitedFiles
using ArgParse
include(joinpath(@__DIR__, "..", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "state_utils.jl"))
include(joinpath(@__DIR__, "..", "Numerics", "gw_rates.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "--outdir";    arg_type = String; required = true
    "--fname";     arg_type = String; required = true; help = "e.g. FullRel_fa_1.0e18_ma_..._Nmax_5.dat"
    "--d_kpc";     arg_type = Float64; default = 1.0
    "--rel_floor"; arg_type = Float64; default = 1e-8; help = "drop lines with peak P below this × brightest"
    "--n_out";     arg_type = Int; default = 2000; help = "number of log-spaced output times"
    "--gw_model";  arg_type = String; default = "nonrel"
    "--no_cloud_shifts"; action = :store_true; help = "omit self-gravity and self-interaction frequency shifts"
end
args = parse_args(s)
dir, fname = args["outdir"], args["fname"]

m = match(r"fa_([0-9.eE+-]+)_ma_([0-9.eE+-]+)_MBH_.*_Nmax_(\d+)", fname)
m === nothing && error("cannot parse f_a / m_a / Nmax from $fname")
fa, mu, Nmax = parse(Float64, m[1]), parse(Float64, m[2]), parse(Int, m[3])

rd(p) = readdlm(joinpath(dir, p * fname))
timeT = vec(rd("Time_")); states = rd("States_"); modes_m = rd("Modes_")
spin = vec(rd("Spin_")); mass = vec(rd("MassBH_"))
ndims(modes_m) == 1 && (modes_m = reshape(modes_m, 1, :))
ndims(states) == 1 && (states = reshape(states, 1, :))
modes = [(round(Int, modes_m[i, 1]), round(Int, modes_m[i, 2]), round(Int, modes_m[i, 3])) for i in 1:size(modes_m, 1)]

tgw, lines = gw_lines(timeT, states, modes, spin, mass, mu, Nmax; d_kpc=args["d_kpc"],
                      rel_floor=args["rel_floor"], n_out=args["n_out"], gw_model=Symbol(args["gw_model"]),
                      self_gravity=!args["no_cloud_shifts"], fa=args["no_cloud_shifts"] ? nothing : fa)
println("$(length(lines)) GW lines (m_a = $mu eV, Nmax = $Nmax, d = $(args["d_kpc"]) kpc)")

meta = [Any[String(L.kind), L.a..., L.b..., maximum(L.P), maximum(L.h0), L.f[1], L.f[end]] for L in lines]
open(joinpath(dir, "GWLines_" * fname), "w") do io
    println(io, "# kind na la ma nb lb mb peak_P_erg_s peak_h0 f_start_Hz f_end_Hz")
    for r in meta
        println(io, join(r, " "))
    end
end
writedlm(joinpath(dir, "GWTime_" * fname), tgw)
writedlm(joinpath(dir, "GWFreq_" * fname), [L.f for L in lines])
writedlm(joinpath(dir, "GWFreqCloud_" * fname), [L.df_cloud for L in lines])
writedlm(joinpath(dir, "GWPower_" * fname), [L.P for L in lines])
writedlm(joinpath(dir, "GWStrain_" * fname), [L.h0 for L in lines])
for L in lines[1:min(10, end)]
    println(rpad("$(L.kind) $(L.a) $(L.kind == :annihilation ? "x" : "->") $(L.b)", 44),
            "  peak h0 = ", round(maximum(L.h0), sigdigits=3), "   f = ", round(L.f[1], sigdigits=8), " -> ", round(L.f[end], sigdigits=8),
            " Hz  (cloud shift range ", round(minimum(L.df_cloud), sigdigits=3), " .. ", round(maximum(L.df_cloud), sigdigits=3), " Hz)")
end
