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
            help     = "Initial BH spin (dimensionless)"
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
    end
    return parse_args(s)
end

parsed  = parse_commandline()

MassBH  = parsed["MassBH"]
SpinBH  = parsed["SpinBH"]
f_a     = parsed["f_a"]
alpha   = parsed["alpha"]
Nmax    = parsed["Nmax"]
tau_max = parsed["tau_max"]
outdir  = parsed["outdir"]

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
println("================================================")

# Output filename -- mirrors single_BH.jl naming convention with Nmax appended
fname = "FullRel_fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH)_Nmax_$(Nmax).dat"

# Skip if all five output files already exist (makes re-runs idempotent)
if isfile(joinpath(outdir, "Spin_" * fname))
    println("Output already exists, skipping: ", joinpath(outdir, "Spin_" * fname))
    exit(0)
end

# Create output directory (no-op if it already exists)
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
cheby            = false
non_rel          = false
high_p           = true

timeT, StatesOut, idx_lvl, spin, massB = @time solve_system(
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
)

println("StatesOut size: ", size(StatesOut))
println("Modes:  ", idx_lvl)
println("Final:  t=", timeT[end], "  spin=", spin[end], "  massB=", massB[end])

writedlm(joinpath(outdir, "Time_"   * fname), timeT)
writedlm(joinpath(outdir, "States_" * fname), StatesOut)
writedlm(joinpath(outdir, "Modes_"  * fname), idx_lvl)
writedlm(joinpath(outdir, "Spin_"   * fname), spin)
writedlm(joinpath(outdir, "MassBH_" * fname), massB)

println("Saved output to: ", outdir)
println("  ", fname)
println("Finished: ", Dates.now())
