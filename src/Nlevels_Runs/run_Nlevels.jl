using Random
using DelimitedFiles
using Suppressor
@suppress include("../super_rad.jl")
using ArgParse
using Dates

# GNew is now in scope from Core/constants.jl (loaded via super_rad.jl)

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--MassBH"
            arg_type = Float64
            required = true
            help = "BH mass in solar masses"
        "--SpinBH"
            arg_type = Float64
            default  = 0.99
            help     = "Initial BH spin (dimensionless); also used for rate tables on resume"
        "--f_a"
            arg_type = Float64
            required = true
            help     = "Axion decay constant (GeV)"
        "--alpha"
            arg_type = Float64
            required = true
            help     = "alpha = GNew * MassBH * m_a"
        "--Nmax"
            arg_type = Int
            required = true
            help     = "Maximum number of levels included"
        "--tau_max"
            arg_type = Float64
            required = true
            help     = "Integration time (natural units)"
        "--outdir"
            arg_type = String
            required = true
            help     = "Directory in which to write output files"
        "--lm_only"
            arg_type = Bool
            default  = false
            help     = "Use l==m-only rate input file (load_rate_input_Nmax_X_lm.txt)"
        "--gw_model"
            arg_type = String
            default  = "nonrel"
            help     = "GW emission: nonrel (NR annihilation+transition rates), rel (not implemented), off"
        "--f_edd"
            arg_type = Float64
            default  = 0.0
            help     = "Constant accretion: Mdot = f_edd * Mdot_Edd(MassBH), fixed in M_sun/yr for the whole run " *
                       "(Core/accretion.jl; thin disc, Bardeen spin-up, Mdot_Edd = L_Edd/(acc_eta c^2)). 0 = off"
        "--acc_eta"
            arg_type = Float64
            default  = 0.1
            help     = "Radiative efficiency in Mdot_Edd (only used if f_edd > 0)"
        "--no_track_alpha"
            action   = :store_true
            help     = "Keep the SR and scattering rates at the initial alpha (original behaviour; " *
                       "default: they follow M(t), see solve_system track_alpha)"
        "--sr_interp"
            arg_type = String
            default  = "anchored"
            help     = "SR-rate interpolation of the Leaver tables (Numerics/rate_interpolation.jl): anchored " *
                       "(default; monotone cubic along the distance from threshold, zero at the tabulated a_c(alpha)), " *
                       "factorised (2026-11-09/10 default; tag _srfac) or legacy (linear in alpha, as before " *
                       "2026-11-09). legacy output files get an _srlin tag"
        "--save_factor"
            arg_type = Float64
            default  = 5.0
            help     = "Also save a point whenever a relevant level's occupation changed by more than this factor " *
                       "since the last saved point (Inf = off; see solve_system)"
        "--save_dlna"
            arg_type = Float64
            default  = 1e-3
            help     = "Also save a point whenever |d ln a| since the last saved point exceeds this (Inf = off)"
        "--emax2_mode"
            arg_type = String
            default  = "auto"
            help     = "Cap on the 211 cloud (solve_system emax2_mode): auto = none (default since 2026-11-10: no cap), fixed " *
                       "(hard switch at the initial Emax2; the earlier default) or evolving (cap from the current BH " *
                       "state with a smooth taper; the earlier default with accretion)"
        "--t_save_start"
            arg_type = Float64
            default  = 0.0
            help     = "Start of the log-spaced save grid in years (0 = default, max(1, 1e-9 tau_max)); use e.g. 1e-6 to " *
                       "resolve the first growth of fast levels. A resume must use the same value"
        "--self_grav"
            arg_type = Bool
            default  = false
            help     = "Add the Newtonian self-gravity contribution (and its interference with self-interactions) to the " *
                       "scattering rates, from rate_sve/*_LvrHc_SG.dat (Compute_all_rates.jl --self_grav true). " *
                       "Output files get an _SG tag"
        "--si_sign"
            arg_type = Float64
            default  = 1.0
            help     = "Sign of the quartic self-interaction in the SI-SG interference (+1 attractive/axion, -1 repulsive)"
        "--rate_input"
            arg_type = String
            default  = ""
            help     = "Custom channel list (file in rate_sve/ or absolute path) instead of load_rate_input_Nmax_X.txt, " *
                       "e.g. from rate_sve/gen_custom_input.py. Output files get an _rin_<name> tag"
        "--track_alpha_refresh"
            arg_type = Float64
            default  = 0.0
            help     = "Change in ln M that triggers a refresh of the scattering rates under track_alpha " *
                       "(0 = every RHS call, default; the former default 1e-5 puts jumps in the RHS, see solve_system)"
        "--solver"
            arg_type = String
            default  = "radau"
            help     = "radau (RadauIIA5, default; use with --analytic_jac) or trbdf2"
        "--analytic_jac"
            action   = :store_true
            help     = "Analytic Jacobian for the level columns instead of finite differences in every column"
        "--radau_kappa"
            arg_type = Float64
            default  = 0.0
            help     = "radau only: Newton stopping tolerance kappa of RadauIIA5 (0 = the OrdinaryDiffEq default, " *
                       "1/100). Smaller values make a step whose first Newton correction is large take another " *
                       "iteration instead of being accepted unchecked"
        "--maxiters"
            arg_type = Float64
            default  = 5e6
            help     = "Maximum solver iterations before the run stops (then marked incomplete; continue with --resume)"
        "--timescale_saves"
            arg_type = Bool
            default  = true
            help     = "With trbdf2, also save a point whenever the timescale callback fires (default; in strongly " *
                       "oscillating phases that is nearly every step). false saves only the regular grid and " *
                       "save-on-change points; the trajectory is unchanged"
        "--gw_table_nmax"
            arg_type = String
            default  = ""
            help     = "Take the GW channels from gw_nr_rates_Nmax_<this>.txt instead of the run's Nmax " *
                       "(empty = same as Nmax; e.g. 18 or a custom table tag like 45_T1020); " *
                       "levels outside that table get no GW"
        "--regime_split"
            arg_type = Bool
            default  = true
            help     = "Use rate_sve/regime_split.txt: channels whose fourth quantum changes between bound (BH) " *
                       "and escaping (Inf) at alpha_c use their BH half on one side and Inf half on the other " *
                       "(false: each channel keeps its list label at every alpha, as before 2026-10-22)"
        "--keep_levels"
            arg_type = String
            default  = ""
            help     = "Truncate the basis (approximate): \"lm\" keeps l == m levels, optionally plus a list, " *
                       "e.g. \"lm,8-7-6,8-7-5\" (pair with a matching --rate_input). Empty = full basis"
        "--dlnM_sc"
            arg_type = Float64
            default  = 0.05
            help     = "Spacing in ln M of the mass grid on which the scattering rates are tabulated under track_alpha " *
                       "(see solve_system track_alpha_dlnM); use ~0.002 to resolve an SI-SG cancellation in f_a"
        "--resume"
            action   = :store_true
            help     = "Continue from last row of existing Time_/States_/Spin_/MassBH_ .dat"
        "--n_times"
            arg_type = Int
            default  = 1000000
            help     = "Number of log-spaced output times in (t_max*1e-9, t_max); transients are still resolved by the " *
                       "save_factor / save_dlna adaptive saves. 1e6 (default) writes ~14 GB per Nmax 15 run, 2e4 ~0.3 GB"
    end
    return parse_args(s)
end

"""Map (n,l,m) rows in Modes_ to indices in the full Nmax spectrum."""
function _mode_key(row)::Tuple{Int,Int,Int}
    return (round(Int, row[1]), round(Int, row[2]), round(Int, row[3]))
end

"""
Read the last saved point of a run (resume), streaming the .dat files: a States_ row holds
one value per saved time (~1e6), so the history is never loaded into memory.

Archives may be pruned (only active modes kept). Inactive modes are restored
to e_floor / e_init so u0 matches setup_quantum_levels_standard(Nmax, ...).
Returns t_last, u0 (full Nmax spectrum), idx_lvl, modes_full, the Modes_ rows of
the files, the number of saved times, and the last spin and mass.
"""
function load_resume_state(outdir::String, fname::String, e_floor::Float64,
                           Nmax::Int, fa::Float64, alph::Float64, aBH::Float64)
    paths = Dict(
        :time   => joinpath(outdir, "Time_"   * fname),
        :states => joinpath(outdir, "States_" * fname),
        :modes  => joinpath(outdir, "Modes_"  * fname),
        :spin   => joinpath(outdir, "Spin_"   * fname),
        :mass   => joinpath(outdir, "MassBH_" * fname),
    )
    for (k, p) in paths
        isfile(p) || error("--resume requested but missing $k file: $p")
    end

    n_times = _count_rows(paths[:time])
    for k in (:spin, :mass)
        n = _count_rows(paths[k])
        n == n_times || error("$(k) has $n rows but Time has $n_times")
    end
    t_last = parse(Float64, strip(_last_row(paths[:time])))
    a_last = parse(Float64, strip(_last_row(paths[:spin])))
    M_last = parse(Float64, strip(_last_row(paths[:mass])))

    modes_old = readdlm(paths[:modes])
    if ndims(modes_old) == 1
        modes_old = reshape(modes_old, 1, :)
    end

    idx_lvl, _, _, modes_full = setup_quantum_levels_standard(Nmax, fa, M_pl, alph, aBH)
    full_index = Dict{Tuple{Int,Int,Int}, Int}()
    for (i, md) in enumerate(modes_full)
        full_index[(Int(md[1]), Int(md[2]), Int(md[3]))] = i
    end

    # Last value of each States_ row (inactive / pruned modes stay at e_floor).
    E_last = fill(e_floor, idx_lvl)
    n_kept = 0
    for ln in eachline(paths[:states])
        isempty(strip(ln)) && continue
        n_kept += 1
        n_kept <= size(modes_old, 1) || error("States has more rows than Modes ($(size(modes_old, 1)))")
        key = _mode_key(modes_old[n_kept, :])
        haskey(full_index, key) || error("Resume mode $key not in Nmax=$Nmax spectrum")
        E_last[full_index[key]] = max(_last_number(ln), e_floor)
    end
    n_kept == size(modes_old, 1) || error("States rows $n_kept != Modes rows $(size(modes_old, 1))")
    n_kept < idx_lvl && println("  pruned archive: mapped $n_kept/$idx_lvl modes; ",
                                "padded $(idx_lvl - n_kept) inactive with e_floor")

    u0 = log.([E_last; max(a_last, e_floor); max(M_last, e_floor)])
    return t_last, u0, idx_lvl, modes_full, modes_old, n_times, a_last, M_last
end

function _count_rows(path)
    n = 0
    for ln in eachline(path)
        isempty(strip(ln)) || (n += 1)
    end
    return n
end

function _last_row(path)
    last = ""
    for ln in eachline(path)
        isempty(strip(ln)) || (last = ln)
    end
    return last
end

# Last tab/space-separated number of a (long) line.
function _last_number(line::AbstractString)
    s = rstrip(line)
    i = findlast(c -> c == '\t' || c == ' ', s)
    return parse(Float64, i === nothing ? s : s[nextind(s, i):end])
end

"""
Append a resumed segment to the run's .dat files in place. Time_/Spin_/MassBH_ get the
new rows; each States_ row (in the solver's mode order, modes_out) is the old row followed
by the new values, written row by row from the old file, so memory stays at one row.
A mode missing from the old files (pruned archive) gets e_floor at the earlier times.
"""
function append_resumed_segment(outdir, fname, modes_old, n_old, e_floor, t_last,
                                timeT, StatesOut, modes_out, spin, massB)
    start = (!isempty(timeT) && timeT[1] == t_last) ? 2 : 1   # drop the duplicate t_last
    for (prefix, v) in (("Time_", timeT), ("Spin_", spin), ("MassBH_", massB))
        open(joinpath(outdir, prefix * fname), "a") do io
            for x in @view v[start:end]
                println(io, x)
            end
        end
    end

    path = joinpath(outdir, "States_" * fname)
    offsets = Int[]
    open(path) do io
        while !eof(io)
            pos = position(io)
            isempty(strip(readline(io))) || push!(offsets, pos)
        end
    end
    old_row = Dict(_mode_key(modes_old[j, :]) => j for j in 1:size(modes_old, 1))
    tmp = path * ".tmp"
    open(path) do src
        open(tmp, "w") do dst
            for (j, md) in enumerate(modes_out)
                key = (Int(md[1]), Int(md[2]), Int(md[3]))
                if haskey(old_row, key)
                    seek(src, offsets[old_row[key]])
                    write(dst, rstrip(readline(src)))
                else
                    for k in 1:n_old
                        k > 1 && write(dst, '\t')
                        print(dst, e_floor)
                    end
                end
                for x in @view StatesOut[j][start:end]
                    write(dst, '\t')
                    print(dst, x)
                end
                write(dst, '\n')
            end
        end
    end
    mv(tmp, path; force=true)
    writedlm(joinpath(outdir, "Modes_" * fname), modes_out)
end

parsed  = parse_commandline()

MassBH  = parsed["MassBH"]
SpinBH  = parsed["SpinBH"]
f_a     = parsed["f_a"]
alpha   = parsed["alpha"]
Nmax    = parsed["Nmax"]
tau_max = parsed["tau_max"]
outdir  = parsed["outdir"]
lm_only = parsed["lm_only"]
do_resume = parsed["resume"]
f_edd   = parsed["f_edd"]
acc_eta = parsed["acc_eta"]
f_edd >= 0 || error("--f_edd must be >= 0")
self_grav  = parsed["self_grav"]
si_sign    = parsed["si_sign"]
rate_input = isempty(parsed["rate_input"]) ? nothing : parsed["rate_input"]
sr_interp = Symbol(parsed["sr_interp"])
sr_interp in (:anchored, :factorised, :legacy) || error("--sr_interp must be anchored, factorised or legacy")
RateInterpolation.SR_RATE_INTERP[] = sr_interp
solver  = Symbol(parsed["solver"])
solver in (:trbdf2, :radau) || error("--solver must be trbdf2 or radau")

# --keep_levels "lm[,n-l-m,...]": l == m levels plus the listed extras
keep_spec = parsed["keep_levels"]
keep_levels = if isempty(keep_spec)
    nothing
else
    toks = strip.(split(keep_spec, ","))
    keep_lm = "lm" in toks
    # "basisK": every level of the Nmax = K basis (truncation modes included)
    base = Set{NTuple{3,Int}}()
    for t in filter(t -> startswith(t, "basis"), toks)
        _, _, _, mK = setup_quantum_levels_standard(parse(Int, t[6:end]), f_a, M_pl, alpha, SpinBH)
        union!(base, Set((Int(md[1]), Int(md[2]), Int(md[3])) for md in mK))
    end
    extras = Set(Tuple(parse.(Int, split(t, "-"))) for t in toks if t != "lm" && !startswith(t, "basis"))
    (n, l, m) -> (keep_lm && l == m) || (n, l, m) in base || (n, l, m) in extras
end

# Derive axion mass from alpha = GNew * MassBH * m_a
m_a = alpha / (GNew * MassBH)

println("================================================")
println("  run_Nlevels.jl  started:  ", Dates.now())
println("  MassBH  = ", MassBH,  "  M_sun")
println("  SpinBH  = ", SpinBH)
println("  f_a     = ", f_a,     "  GeV")
println("  alpha   = ", alpha,   "  (GNew*M*m_a = ", GNew*MassBH*m_a, ")")
println("  m_a     = ", m_a,     "  eV")
println("  Nmax    = ", Nmax)
println("  tau_max = ", tau_max)
println("  outdir  = ", outdir)
println("  lm_only = ", lm_only)
println("  resume  = ", do_resume)
isempty(keep_spec) || println("  keep_levels = ", keep_spec)
isempty(parsed["gw_table_nmax"]) || println("  gw_table_nmax = ", parsed["gw_table_nmax"])
println("  regime_split = ", parsed["regime_split"])
println("  sr_interp = ", sr_interp)
println("  n_times = ", parsed["n_times"])
println("  self_grav = ", self_grav, self_grav ? "  (si_sign = $(si_sign))" : "")
println("  rate_input = ", rate_input === nothing ? "default (load_rate_input_Nmax_$(Nmax).txt)" : rate_input)
println("  solver  = ", solver, parsed["analytic_jac"] ? " (analytic Jacobian)" : "",
        parsed["radau_kappa"] > 0 ? " (kappa = $(parsed["radau_kappa"]))" : "",
        "  track_alpha_refresh = ", parsed["track_alpha_refresh"])
println("  f_edd   = ", f_edd, f_edd > 0 ? "  (Mdot = $(f_edd * eddington_rate(MassBH; eta=acc_eta)) M_sun/yr, eta = $(acc_eta))" : "  (no accretion)")
println("================================================")

# Output filename -- mirrors single_BH.jl naming convention with Nmax appended.
# With accretion a _fedd_<f_edd> tag goes before _Nmax_ (the *Nmax_N.dat globs and
# the gw_lines_from_run.jl filename parser still match); without it the name is unchanged.
fedd_tag = f_edd > 0 ? "_fedd_$(f_edd)" : ""
# _SG (self-gravity on), _rin_<list> (custom channel list) and _srlin (legacy SR-rate interpolation) tags also
# go before _Nmax_, so SI-only and SG runs, default and custom channel lists, or the two SR interpolations never
# overwrite each other.
sg_tag  = self_grav ? "_SG" : ""
rin_tag = rate_input === nothing ? "" : "_rin_" * splitext(basename(rate_input))[1]
srlin_tag = sr_interp === :legacy ? "_srlin" : sr_interp === :factorised ? "_srfac" : ""
fname = "FullRel_fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH)$(fedd_tag)$(sg_tag)$(rin_tag)$(srlin_tag)_Nmax_$(Nmax).dat"

spin_path = joinpath(outdir, "Spin_" * fname)
if isfile(spin_path) && !do_resume
    println("Output already exists, skipping: ", spin_path)
    exit(0)
end

mkpath(outdir)

# Fixed numerical parameters -- taken from single_BH.jl
impose_low_cut   = 1e-100
return_all_info  = true
n_times          = parsed["n_times"]
eq_threshold     = 1e-100
stop_on_a        = 0.0
abstol           = 1e-30
N_pts_interp     = 100
N_pts_interpL    = 100
cheby            = true
non_rel          = false
high_p           = true

e_floor = 1.0 / (GNew * MassBH^2 * M_to_eV)

u0_override = nothing
t_start = 0.0
t_last = 0.0
modes_old = nothing
n_old = 0

if do_resume
    t_last, u0_override, n_modes, modes_full, modes_old, n_old, spin_last, mass_last =
        load_resume_state(outdir, fname, e_floor, Nmax, f_a, alpha, SpinBH)
    t_start = t_last
    if keep_levels !== nothing
        # Truncated basis: the kept levels (basis order, as solve_system does).
        kidx = findall(md -> keep_levels(Int(md[1]), Int(md[2]), Int(md[3])), modes_full)
        u0_override = [u0_override[kidx]; u0_override[end-1:end]]
        n_modes = length(kidx)
    end
    println("Resume from t=", t_last, "  (tau_max=", tau_max, ")")
    println("  n_modes=", n_modes, "  n_times_old=", n_old)
    println("  u0_len=", length(u0_override), "  spin_last=", spin_last,
            "  mass_last=", mass_last)
    if t_last >= tau_max
        println("Already at/past tau_max; nothing to do.")
        exit(0)
    end
    flush(stdout)
end

timeT, StatesOut, modes_out, spin, massB = @time solve_system(
    m_a, f_a, SpinBH, MassBH, tau_max;
    impose_low_cut  = impose_low_cut,
    return_all_info = return_all_info,
    n_times         = n_times,
    eq_threshold    = eq_threshold,
    abstol          = abstol,
    non_rel         = non_rel,
    debug           = true,
    high_p          = high_p,
    N_pts_interp    = N_pts_interp,
    N_pts_interpL   = N_pts_interpL,
    Nmax            = Nmax,
    cheby           = cheby,
    lm_only         = lm_only,
    u0_override     = u0_override,
    t_start         = t_start,
    gw_model        = Symbol(parsed["gw_model"]),
    f_edd           = f_edd,
    acc_eta         = acc_eta,
    track_alpha     = !parsed["no_track_alpha"],
    track_alpha_dlnM = (0.0025, parsed["dlnM_sc"]),
    save_factor     = parsed["save_factor"],
    save_dlna       = parsed["save_dlna"],
    self_grav       = self_grav,
    si_sign         = si_sign,
    rate_input      = rate_input,
    track_alpha_refresh = parsed["track_alpha_refresh"],
    solver          = solver,
    analytic_jac    = parsed["analytic_jac"],
    keep_levels     = keep_levels,
    maxiters        = parsed["maxiters"],
    timescale_saves = parsed["timescale_saves"],
    t_save_start    = parsed["t_save_start"] > 0 ? parsed["t_save_start"] : nothing,
    emax2_mode      = parsed["emax2_mode"] == "auto" ? nothing : Symbol(parsed["emax2_mode"]),
    gw_table_nmax   = isempty(parsed["gw_table_nmax"]) ? nothing : parsed["gw_table_nmax"],
    regime_split    = parsed["regime_split"],
    radau_kappa     = parsed["radau_kappa"] > 0 ? parsed["radau_kappa"] : nothing,
)

println("StatesOut size: ", size(StatesOut))
println("Modes:  ", modes_out)
println("Final:  t=", timeT[end], "  spin=", spin[end], "  massB=", massB[end])

# Always write the full spectrum Modes_ from the solver (not a pruned subset).
if do_resume
    append_resumed_segment(outdir, fname, modes_old, n_old, e_floor, t_last,
                           timeT, StatesOut, modes_out, spin, massB)
else
    writedlm(joinpath(outdir, "Time_"   * fname), timeT)
    writedlm(joinpath(outdir, "States_" * fname), StatesOut)
    writedlm(joinpath(outdir, "Modes_"  * fname), modes_out)
    writedlm(joinpath(outdir, "Spin_"   * fname), spin)
    writedlm(joinpath(outdir, "MassBH_" * fname), massB)
end

println("Saved output to: ", outdir)
println("  ", fname)

# Stopped early (wall-clock limit): mark the run incomplete so finish_nmax_output.sh
# leaves the loose .dat files for --resume instead of pruning/archiving them.
incomplete_marker = joinpath(outdir, "NL_incomplete_Nmax_$(Nmax)")
if timeT[end] < tau_max * (1 - 1e-12)
    write(incomplete_marker, "stopped at t=$(timeT[end]) of tau_max=$(tau_max); continue with --resume\n")
    println("INCOMPLETE: stopped at t=", timeT[end], " of ", tau_max, "; rerun with --resume to continue")
else
    rm(incomplete_marker; force=true)
end
println("Finished: ", Dates.now())
