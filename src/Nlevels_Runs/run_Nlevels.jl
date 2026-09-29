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
            help     = "Axion decay constant (eV)"
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
        "--resume"
            action   = :store_true
            help     = "Continue from last row of existing Time_/States_/Spin_/MassBH_ .dat"
    end
    return parse_args(s)
end

"""Map (n,l,m) rows in Modes_ to indices in the full Nmax spectrum."""
function _mode_key(row)::Tuple{Int,Int,Int}
    return (round(Int, row[1]), round(Int, row[2]), round(Int, row[3]))
end

"""
Load prior trajectory and rebuild log-space u0 from the last saved time.

Archives may be pruned (only active modes kept). Inactive modes are restored
to e_floor / e_init so u0 matches setup_quantum_levels_standard(Nmax, ...).
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

    time_old   = vec(readdlm(paths[:time]))
    states_old = readdlm(paths[:states])          # (n_kept, n_times)
    modes_old  = readdlm(paths[:modes])
    spin_old   = vec(readdlm(paths[:spin]))
    mass_old   = vec(readdlm(paths[:mass]))

    n_kept, n_times = size(states_old)
    if length(time_old) != n_times
        error("Time length $(length(time_old)) != States cols $n_times")
    end
    if length(spin_old) != n_times || length(mass_old) != n_times
        error("Spin/MassBH length mismatch vs States cols $n_times")
    end
    if ndims(modes_old) == 1
        modes_old = reshape(modes_old, 1, :)
    end
    if size(modes_old, 1) != n_kept
        error("Modes rows $(size(modes_old,1)) != States rows $n_kept")
    end

    idx_lvl, _, _, modes_full = setup_quantum_levels_standard(Nmax, fa, M_pl, alph, aBH)
    full_index = Dict{Tuple{Int,Int,Int}, Int}()
    for (i, md) in enumerate(modes_full)
        full_index[(Int(md[1]), Int(md[2]), Int(md[3]))] = i
    end

    # Expand pruned States to full spectrum (inactive = e_floor for all times).
    states_full = fill(e_floor, idx_lvl, n_times)
    n_mapped = 0
    for j in 1:n_kept
        key = _mode_key(modes_old[j, :])
        haskey(full_index, key) || error("Resume mode $key not in Nmax=$Nmax spectrum")
        i = full_index[key]
        @inbounds for t in 1:n_times
            states_full[i, t] = Float64(states_old[j, t])
        end
        n_mapped += 1
    end
    if n_kept < idx_lvl
        println("  pruned archive: mapped $n_mapped/$idx_lvl modes; ",
                "padded $(idx_lvl - n_mapped) inactive with e_floor")
    elseif n_kept != idx_lvl
        error("States has $n_kept modes but Nmax=$Nmax expects $idx_lvl")
    end

    t_last = Float64(time_old[end])
    E_last = [max(Float64(states_full[j, end]), e_floor) for j in 1:idx_lvl]
    a_last = max(Float64(spin_old[end]), e_floor)
    M_last = max(Float64(mass_old[end]), e_floor)
    u0 = log.([E_last; a_last; M_last])

    return time_old, states_full, modes_full, spin_old, mass_old, t_last, u0, idx_lvl
end

function concat_trajectories(time_old, states_old, spin_old, mass_old,
                             time_new, states_new, spin_new, mass_new)
    # Drop duplicate t_last if the new segment starts at the same time.
    start = 1
    if !isempty(time_new) && !isempty(time_old) && time_new[1] == time_old[end]
        start = 2
    end
    if start > length(time_new)
        return time_old, states_old, spin_old, mass_old
    end

    time_out = vcat(time_old, time_new[start:end])
    spin_out = vcat(spin_old, spin_new[start:end])
    mass_out = vcat(mass_old, mass_new[start:end])

    n_modes = size(states_old, 1)
    # states_new from solve_system is Vector of Vectors (one per mode)
    states_out = similar(states_old, n_modes, length(time_out))
    for j in 1:n_modes
        new_row = states_new[j][start:end]
        states_out[j, :] = vcat(states_old[j, :], new_row)
    end
    return time_out, states_out, spin_out, mass_out
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

# Derive axion mass from alpha = GNew * MassBH * m_a
m_a = alpha / (GNew * MassBH)

println("================================================")
println("  run_Nlevels.jl  started:  ", Dates.now())
println("  MassBH  = ", MassBH,  "  M_sun")
println("  SpinBH  = ", SpinBH)
println("  f_a     = ", f_a,     "  eV")
println("  alpha   = ", alpha,   "  (GNew*M*m_a = ", GNew*MassBH*m_a, ")")
println("  m_a     = ", m_a,     "  eV")
println("  Nmax    = ", Nmax)
println("  tau_max = ", tau_max)
println("  outdir  = ", outdir)
println("  lm_only = ", lm_only)
println("  resume  = ", do_resume)
println("================================================")

# Output filename -- mirrors single_BH.jl naming convention with Nmax appended
fname = "FullRel_fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH)_Nmax_$(Nmax).dat"

spin_path = joinpath(outdir, "Spin_" * fname)
if isfile(spin_path) && !do_resume
    println("Output already exists, skipping: ", spin_path)
    exit(0)
end

mkpath(outdir)

# Fixed numerical parameters -- taken from single_BH.jl
impose_low_cut   = 1e-100
return_all_info  = true
n_times          = 1000000
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
time_old = nothing
states_old = nothing
modes_old = nothing
spin_old = nothing
mass_old = nothing

if do_resume
    time_old, states_old, modes_old, spin_old, mass_old, t_last, u0_override, n_modes =
        load_resume_state(outdir, fname, e_floor, Nmax, f_a, alpha, SpinBH)
    t_start = t_last
    println("Resume from t=", t_last, "  (tau_max=", tau_max, ")")
    println("  n_modes=", n_modes, "  n_times_old=", length(time_old))
    println("  u0_len=", length(u0_override), "  spin_last=", spin_old[end],
            "  mass_last=", mass_old[end])
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
)

if do_resume
    # StatesOut is Vector-of-Vectors (one per mode); convert concat to Matrix for writedlm
    timeT, states_mat, spin, massB = concat_trajectories(
        time_old, states_old, spin_old, mass_old,
        timeT, StatesOut, spin, massB,
    )
    # writedlm Matrix writes rows; match cold-start Vector-of-Vectors layout
    StatesOut = [states_mat[j, :] for j in 1:size(states_mat, 1)]
end
# Always write the full spectrum Modes_ from the solver (not a pruned subset).
modes_write = modes_out

println("StatesOut size: ", size(StatesOut))
println("Modes:  ", modes_out)
println("Final:  t=", timeT[end], "  spin=", spin[end], "  massB=", massB[end])

writedlm(joinpath(outdir, "Time_"   * fname), timeT)
writedlm(joinpath(outdir, "States_" * fname), StatesOut)
writedlm(joinpath(outdir, "Modes_"  * fname), modes_write)
writedlm(joinpath(outdir, "Spin_"   * fname), spin)
writedlm(joinpath(outdir, "MassBH_" * fname), massB)

println("Saved output to: ", outdir)
println("  ", fname)
println("Finished: ", Dates.now())
