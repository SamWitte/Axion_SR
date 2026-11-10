# Rank the GW channels of finished towers45 leaves (NL_output_Nmax_45.tar.gz) by time-integrated energy flux
# (NR rates as in the runs); print, per leaf, the channels with E/Etot > 1e-4 that lack a Kerr table, and at the
# end the union over leaves. Each archive is extracted to a temporary directory and removed afterwards.
#   julia scripts/rank_gw_channels.jl LEAF [LEAF ...]      LEAF = output/<tag>/BH_x/fa_y/alpha_z (towers45 runs, Nmax 45)
using Suppressor
const SRC = joinpath(@__DIR__, "..")
include(joinpath(SRC, "Core", "constants.jl"))
include(joinpath(SRC, "state_utils.jl"))
@suppress include(joinpath(SRC, "Numerics", "gw_rates.jl"))
include(joinpath(SRC, "Numerics", "gw_kerr_tables.jl"))

lvl(s) = "$(s[1])$(s[2])$(s[3])"
const UNION = Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}},Float64}()
function rank(leaf)
    dir = mktempdir()
    run(`tar xzf $(joinpath(leaf, "NL_output_Nmax_45.tar.gz")) -C $dir`)
    g(k) = first(filter(f -> startswith(f, k * "_"), readdir(dir)))
    T = vec(readdlm(joinpath(dir, g("Time"))))
    S = vec(readdlm(joinpath(dir, g("Spin"))))
    M = vec(readdlm(joinpath(dir, g("MassBH"))))
    modes = readdlm(joinpath(dir, g("Modes")))
    U = readdlm(joinpath(dir, g("States")))
    rm(dir; recursive=true)
    α0 = parse(Float64, match(r"alpha_([0-9.]+)", leaf).captures[1])
    mu = α0 / (GNew * M[1])
    md = [(Int(modes[i, 1]), Int(modes[i, 2]), Int(modes[i, 3])) for i in 1:size(modes, 1)]
    c = gw_build_cache("45_towers45", md, mu, M[1], S[1]; gw_model=:nonrel, min_rate_per_yr=0.0)
    tag = replace(leaf, r".*output/" => "")
    println("$(tag): $(length(md)) levels, $(length(c.ann_i)) ann + $(length(c.tr_i)) tr channels, $(length(T)) times, t_end = $(T[end])")
    ks = unique(round.(Int, exp10.(LinRange(0, log10(length(T)), 400))))
    E = Dict{Tuple{Symbol,Int,Int},Float64}()
    tprev = T[ks[1]]
    for k in ks
        α = GNew * M[k] * mu; a = S[k]
        ann, tr = gw_channel_rates(c, α, a)
        dt = T[k] - tprev; tprev = T[k]
        for (i, j, r) in ann
            key = (:ann, min(i, j), max(i, j))
            E[key] = get(E, key, 0.0) + (c.ω[i] + c.ω[j]) * r * U[i, k] * U[j, k] * dt
        end
        for (hi, lo, r) in tr
            key = (:tr, hi, lo)
            E[key] = get(E, key, 0.0) + abs(c.ω[hi] - c.ω[lo]) * r * U[hi, k] * U[lo, k] * dt
        end
    end
    tot = sum(values(E))
    ord = sort(collect(keys(E)), by=k -> -E[k])
    tabs = gw_kerr_load_tables()
    have = tabs === nothing ? Set{String}() : Set(basename.(tabs.files))
    println("  top 8: ", join(["$(k[1] == :ann ? lvl(md[k[2]]) * "x" * lvl(md[k[3]]) : lvl(md[k[2]]) * "->" * lvl(md[k[3]])) $(round(E[k] / tot, sigdigits=2))" for k in ord[1:min(8, end)]], ", "))
    miss = String[]
    for k in ord
        E[k] / tot > 1e-4 || break
        kind, i, j = k
        ch = (kind, md[i], md[j])
        UNION[ch] = max(get(UNION, ch, 0.0), E[k] / tot)
        gw_kerr_table_name(kind, md[i], md[j]) in have || push!(miss, "$(kind) $(join(md[i], '-')) $(join(md[j], '-')) $(round(E[k] / tot, sigdigits=2))")
    end
    println("  $(length(miss)) channels with E/Etot > 1e-4 lack a Kerr table: ", join(miss, "; "))
    flush(stdout)
end
for leaf in ARGS
    try
        rank(leaf)
    catch e
        println("$(leaf): FAILED $(e)")
    end
end
println("\nUNION of channels with E/Etot > 1e-4 in any leaf ($(length(UNION))):")
for (ch, f) in sort(collect(UNION), by=x -> -x[2])
    println("$(ch[1]) $(join(ch[2], '-')) $(join(ch[3], '-'))   max E/Etot = $(round(f, sigdigits=2))")
end
