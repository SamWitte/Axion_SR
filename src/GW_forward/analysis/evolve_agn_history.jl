# Evolve an AGN black hole + its axion cloud along an accretion history exported by
# agn_cloud_resonances.export_histories, and write the evolution in the format
# agn_cloud_resonances.run_external reads.
#
#   julia evolve_agn_history.jl --rundir runs --mu 3e-17 --f_a 1e18 [--Nmax 3] [--ids 0 1 2]
#
# Input:  <rundir>/hist_XXXX.txt   (seed mass and spin, dlnM/dt while ON, episode list)
# Output: <rundir>/evol_XXXX.csv   age_yr, M_Msun, a, Mc_n_l_m (= cloud mass in |nlm> / M)
#         <rundir>/evol_XXXX_meta.json   mu_eV, f_a, Nmax, status, ...
#
# The history is cut into its ON and OFF segments and each is one solve_system
# call: while ON, Mdot = dlnM_dt_on * M(t) (acc_dlnM_dt) with Bardeen spin-up from
# the prograde ISCO; while OFF, no accretion. Superradiance, scattering and GW
# emission act throughout. Each call is normalised to the BH mass at its start
# (u = N / (G M^2)); occupations are rescaled from segment to segment, which
# keeps the mass change inside one call small (rates stay at the right alpha).
# track_alpha is on in every segment, so rates follow M(t) within it as well.
# The 211 cap is emax2_mode = :evolving in every segment (current cloud + what the
# current BH can still extract), so chaining does not reset it.
#
# Only coherent accretion (every episode s = +1) is supported: the solver's spin
# is positive, so a retrograde disc that flips the spin cannot be followed yet.
using ArgParse
using Printf
using Suppressor
@suppress include(joinpath(@__DIR__, "..", "..", "super_rad.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "--rundir";   arg_type = String;  required = true
    "--mu";       arg_type = Float64; required = true; help = "boson mass [eV] (must equal Params.mu_eV)"
    "--f_a";      arg_type = Float64; default = 1e18;  help = "decay constant [GeV]"
    "--Nmax";     arg_type = Int;     default = 3
    "--ids";      arg_type = Int;     nargs = '*';     help = "history ids (default: all hist_*.txt in rundir)"
    "--n_out";    arg_type = Int;     default = 200;   help = "output rows per segment (log-spaced in segment time)"
    "--n_times";  arg_type = Int;     default = 2000;  help = "solver save points per segment"
    "--overwrite"; action = :store_true
end
args = parse_args(s)
haskey(ENV, "NL_MAX_WALL_SEC") || (ENV["NL_MAX_WALL_SEC"] = "3600")

"""Read hist_XXXX.txt -> (Dict of header values, episodes matrix [start_yr end_yr s])."""
function read_history(path)
    hdr = Dict{String, String}()
    rows = Vector{Vector{Float64}}()
    for ln in eachline(path)
        if startswith(ln, "#")
            m = match(r"^#\s*(\w+)\s*=\s*(.*)$", ln)
            m === nothing || (hdr[m[1]] = strip(m[2]))
        elseif !isempty(strip(ln))
            push!(rows, parse.(Float64, split(ln)))
        end
    end
    return hdr, isempty(rows) ? zeros(0, 3) : permutedims(hcat(rows...))
end

"""Up to n indices of t (local segment time), log-spaced, plus first and last."""
function thin(t, n)
    length(t) <= n && return collect(eachindex(t))
    tp = t[t .> 0]
    isempty(tp) && return [1, length(t)]
    targets = exp10.(range(log10(minimum(tp)), log10(t[end]), length=n))
    return unique(vcat(1, [clamp(searchsortedfirst(t, x), 1, length(t)) for x in targets], length(t)))
end

function evolve_history(id, rundir, mu, fa, Nmax; n_out=200, n_times=2000)
    hdr, episodes = read_history(joinpath(rundir, @sprintf("hist_%04d.txt", id)))
    age_today = parse(Float64, hdr["age_today_yr"])
    M = parse(Float64, hdr["M_seed_Msun"])
    a = max(parse(Float64, hdr["a_seed"]), 1e-4)          # the solver evolves log(a)
    r_on = parse(Float64, hdr["dlnM_dt_on_per_yr"])
    all(episodes[:, 3] .== 1) || error("history $id has retrograde episodes (s = -1); only coherent accretion is supported")

    edges = sort(unique(vcat(0.0, episodes[:, 1], episodes[:, 2], age_today)))
    edges = edges[edges .<= age_today]
    is_on(t) = any((episodes[:, 1] .<= t) .& (t .< episodes[:, 2]))

    rows = Vector{Vector{Float64}}()
    modes = nothing
    u = nothing                          # occupations normalised to the current segment's start mass
    status = "ok"
    t0w = time()
    nseg = 0
    for (t1, t2) in zip(edges[1:end-1], edges[2:end])
        dur = t2 - t1
        dur > 0 || continue
        on = is_on(0.5 * (t1 + t2))
        nseg += 1
        u0 = u === nothing ? nothing : log.(vcat(u, a, M))
        local tt, st, md, sp, mb
        try
            tt, st, md, sp, mb = @suppress_out solve_system(mu, fa, a, M, dur;
                return_all_info=true, n_times=n_times, Nmax=Nmax, u0_override=u0,
                acc_dlnM_dt=on ? r_on : 0.0, track_alpha=true, emax2_mode=:evolving,
                impose_low_cut=1e-100, eq_threshold=1e-100, abstol=1e-30, non_rel=false,
                high_p=true, cheby=true, N_pts_interp=100, N_pts_interpL=100)
        catch err
            status = @sprintf("error in segment %d (age %.4g-%.4g yr): %s", nseg, t1, t2, sprint(showerror, err))
            break
        end
        modes = md
        Mnorm = M
        for k in thin(tt, n_out)
            Mk = mb[k]
            # cloud mass of level i: N_i mu = u_i G M_norm^2 mu  [M_sun]; written as a fraction of M
            push!(rows, vcat(t1 + tt[k], Mk, sp[k], [st[i][k] * GNew * mu * Mnorm^2 / Mk for i in eachindex(st)]))
        end
        if !isapprox(tt[end], dur; rtol=1e-6)
            status = @sprintf("incomplete: solver stopped at age %.6g yr in segment %d (%.4g-%.4g yr)",
                              t1 + tt[end], nseg, t1, t2)
            break
        end
        Mnew = mb[end]
        u = [st[i][end] for i in eachindex(st)] .* (Mnorm / Mnew)^2   # renormalise to the next start mass
        M, a = Mnew, max(sp[end], 1e-4)
        @printf("  hist %04d seg %3d %s %.4g-%.4g yr: M=%.4g a=%.4f u_max=%.3g\n", id, nseg, on ? "ON " : "off",
                t1, t2, M, a, maximum(u))
        flush(stdout)
    end

    names = ["age_yr", "M_Msun", "a"]
    modes === nothing || append!(names, ["Mc_$(md[1])_$(md[2])_$(md[3])" for md in modes])
    open(joinpath(rundir, @sprintf("evol_%04d.csv", id)), "w") do io
        println(io, join(names, ","))
        for r in rows
            println(io, join([@sprintf("%.10g", x) for x in r], ","))
        end
    end
    open(joinpath(rundir, @sprintf("evol_%04d_meta.json", id)), "w") do io
        @printf(io, "{\"id\": %d, \"mu_eV\": %.10g, \"f_a\": %.6g, \"Nmax\": %d, \"status\": \"%s\", ", id, mu, fa, Nmax,
                replace(status, "\"" => "'"))
        @printf(io, "\"n_segments\": %d, \"M_end_Msun\": %.8g, \"a_end\": %.6g, \"M_today_target_Msun\": %s, \"wall_s\": %.1f}\n",
                nseg, M, a, hdr["M_today_target_Msun"], time() - t0w)
    end
    return status, nseg, M, a, time() - t0w
end

rundir = args["rundir"]
ids = isempty(args["ids"]) ?
      sort([parse(Int, f[6:9]) for f in readdir(rundir) if startswith(f, "hist_") && endswith(f, ".txt")]) :
      args["ids"]
for id in ids
    out = joinpath(rundir, @sprintf("evol_%04d.csv", id))
    if isfile(out) && !args["overwrite"]
        println("exists, skipping: $out")
        continue
    end
    st, nseg, M, a, wall = evolve_history(id, rundir, args["mu"], args["f_a"], args["Nmax"];
                                          n_out=args["n_out"], n_times=args["n_times"])
    @printf("hist %04d: %s, %d segments, M_end = %.4g M_sun, a_end = %.4f, %.0f s\n", id, st, nseg, M, a, wall)
    flush(stdout)
end
