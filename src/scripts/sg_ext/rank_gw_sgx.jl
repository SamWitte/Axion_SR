# Rank the GW channels of finished extended-SG runs by time-integrated energy flux (NR rates, table "sgx") and list
# those with E/Etot > 1e-4 that have no Kerr (Teukolsky) table: the channels worth building before trusting gw_rel.
#   julia rank_gw_sgx.jl RUN_DIR [RUN_DIR ...]      (adapted from src/scripts/rank_gw_channels.jl, loose .dat files)
using Suppressor, DelimitedFiles
const SRC = "/data/dust/user/samwitte/Axion_SR_sg/src"
include(joinpath(SRC, "Core", "constants.jl"))
include(joinpath(SRC, "state_utils.jl"))
@suppress include(joinpath(SRC, "Numerics", "gw_rates.jl"))
include(joinpath(SRC, "Numerics", "gw_kerr_tables.jl"))
lvl(s) = "$(s[1])-$(s[2])-$(s[3])"
const UNION = Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}},Float64}()
function rank(dir)
    g(k) = joinpath(dir, first(filter(f -> startswith(f, k * "_"), readdir(dir))))
    T = vec(readdlm(g("Time"))); S = vec(readdlm(g("Spin"))); M = vec(readdlm(g("MassBH")))
    modes = readdlm(g("Modes")); U = readdlm(g("States"))
    α0 = parse(Float64, match(r"alpha_([0-9.]+)", dir).captures[1])
    mu = α0 / (GNew * M[1])
    md = [(Int(modes[i, 1]), Int(modes[i, 2]), Int(modes[i, 3])) for i in 1:size(modes, 1)]
    c = gw_build_cache("sgx", md, mu, M[1], S[1]; gw_model=:nonrel, min_rate_per_yr=0.0)
    ks = unique(round.(Int, exp10.(LinRange(0, log10(length(T)), 600))))
    E = Dict{Tuple{Symbol,Int,Int},Float64}(); tprev = T[ks[1]]
    for k in ks
        α = GNew * M[k] * mu; a = S[k]
        ann, tr = gw_channel_rates(c, α, a)
        dt = T[k] - tprev; tprev = T[k]
        for (i, j, r) in ann
            key = (:ann, min(i, j), max(i, j)); E[key] = get(E, key, 0.0) + (c.ω[i] + c.ω[j]) * r * U[i, k] * U[j, k] * dt
        end
        for (hi, lo, r) in tr
            key = (:tr, hi, lo); E[key] = get(E, key, 0.0) + abs(c.ω[hi] - c.ω[lo]) * r * U[hi, k] * U[lo, k] * dt
        end
    end
    tot = sum(values(E)); ord = sort(collect(keys(E)), by=k -> -E[k])
    tabs = gw_kerr_load_tables(); have = tabs === nothing ? Set{String}() : Set(basename.(tabs.files))
    nm(k) = k[1] == :ann ? lvl(md[k[2]]) * " x " * lvl(md[k[3]]) : lvl(md[k[2]]) * " -> " * lvl(md[k[3]])
    println(replace(dir, r".*runs/" => ""), ": top ", join(["$(nm(k)) $(round(E[k] / tot, sigdigits=2))" for k in ord[1:min(6, end)]], ", "))
    for k in ord
        E[k] / tot > 1e-4 || break
        ch = (k[1], md[k[2]], md[k[3]])
        gw_kerr_table_name(ch...) in have || (UNION[ch] = max(get(UNION, ch, 0.0), E[k] / tot))
    end
end
for d in ARGS
    try rank(d) catch e; println(d, ": FAILED ", first(sprint(showerror, e), 200)) end
end
println("\nchannels with E/Etot > 1e-4 in some run and no Kerr table ($(length(UNION))):")
for (ch, f) in sort(collect(UNION), by=x -> -x[2])
    println("$(ch[1]) $(lvl(ch[2])) $(lvl(ch[3]))   $(round(f, sigdigits=2))")
end
