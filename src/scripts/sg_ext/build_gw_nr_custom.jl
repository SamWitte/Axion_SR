# Non-relativistic GW annihilation/transition table for an explicit level list (same rows and pruning as
# gw_build_nr_table in Numerics/gw_rates.jl, which only knows the full Nmax basis):
#   julia build_gw_nr_custom.jl LEVELS.txt TAG      -> <repo>/src/rate_sve/gw_nr_rates_Nmax_<TAG>.txt
# LEVELS.txt: one level per line, "n l m" or "n-l-m". Use with run_Nlevels.jl --gw_table_nmax <TAG>.
const SRC = get(ENV, "AXSR_SRC", "/data/dust/user/samwitte/Axion_SR_sg/src")
include(joinpath(SRC, "Core", "constants.jl"))
include(joinpath(SRC, "state_utils.jl"))
include(joinpath(SRC, "Numerics", "gw_rates.jl"))

function build(levels_file, tag; floor_at_alpha1=1e-60, rel_channel=1e-8)
    states = NTuple{3,Int}[]
    for ln in eachline(levels_file)
        s = strip(ln); (isempty(s) || startswith(s, "#")) && continue
        p = parse.(Int, split(replace(s, "-" => " ")))
        push!(states, (p[1], p[2], p[3]))
    end
    states = sort(unique(states))
    pairs = [(states[i], states[j]) for i in 1:length(states) for j in i:length(states)]
    radcache = Dict{NTuple{5,Int}, Float64}(); angcache = Dict{NTuple{5,Int}, Float64}()
    path = gw_table_path(tag); t0 = time(); nrow = 0
    println("GW NR table: $(length(states)) states, $(length(pairs)) pairs -> $(path)")
    open(path, "w") do io
        println(io, "# GW NR rates, custom level list $(basename(levels_file)) ($(length(states)) levels), solver convention (Numerics/gw_rates.jl).")
        println(io, "# ann: rate = C alpha^p.  tr: rate = K alpha^(2-2L mass | 4-2L current) (dw/mu)^(2L+1)")
        println(io, "# type na la ma nb lb mb coeff power kind")
        for (ip, (sa, sb)) in enumerate(pairs)
            C, p = gw_annihilation_leading(sa, sb)
            if C >= floor_at_alpha1
                println(io, join(Any["ann", sa..., sb..., C, p, 0], " ")); nrow += 1
            end
            if sa != sb
                ch = gw_transition_channels(sa, sb; radcache=radcache, angcache=angcache)
                if !isempty(ch)
                    d1 = gw_delta_bound(sa, sb)
                    r1 = [K * d1^(2L + 1) for (L, kind, K) in ch]
                    keep = max(floor_at_alpha1, rel_channel * maximum(r1))
                    for ((L, kind, K), r) in zip(ch, r1)
                        if r >= keep
                            println(io, join(Any["tr", sa..., sb..., K, L, kind == :mass ? 1 : 2], " ")); nrow += 1
                        end
                    end
                end
            end
            ip % 2000 == 0 && (println("  $(ip)/$(length(pairs)) pairs  ($(round(time() - t0, digits=1)) s)"); flush(stdout))
        end
    end
    println("wrote $(nrow) rows in $(round(time() - t0, digits=1)) s")
end
build(ARGS[1], ARGS[2])
