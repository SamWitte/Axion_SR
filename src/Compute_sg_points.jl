include("solve_sr_rates.jl")
include("Core/constants.jl")
include("state_utils.jl")
using ArgParse
using DelimitedFiles

# Self-interaction + self-gravity rate of a list of scattering channels at an explicit list of alpha values, with the
# settings of Compute_all_rates.jl (a = 0.95, Leaver/Heun bound states, BHlmax 2, Gauss-Legendre angles, BH recoil)
# and three changes that make large channel lists with high-n levels tractable and accurate at low alpha:
#   --rel_nmax N   levels with n <= N use the relativistic (Leaver/Heun) bound state as in the existing tables, levels
#                  with n > N the hydrogenic one (gf_radial NR_states). The relativistic solve is slow and fragile at
#                  n >~ 10, where the hydrogenic state is accurate; for low n it is needed (hydrogenic 322, 211 change
#                  the self-gravity rate of 322 322 -> 211 + inf by 2x at alpha = 0.2).
#   --alpha_h c    when the fourth quantum is bound, the step of its radial integration (h = r+/h_mve in r*) is capped
#                  at c/alpha, i.e. >= 2 pi/c points per near-horizon wavelength. The default step grows like
#                  n_max^2/alpha^2 and is not converged at low alpha (commit 2767135; pilot of 2026-11-12:
#                  433 433 -> 10-6-6 off by 65% at alpha 0.05, 322 322 -> 544 by 12%). 0 = default step only.
#   --cache DIR    gf_radial caches (eigenvalue searches and relativistic radial solves, GF_CACHE_DIR), shared on disk.
# Never touches rate_sve/*_LvrHc_*.dat: each channel goes to <out>/<S1>_<S2>_<S3>_<S4>.pts with rows
#   alpha  G_SI  G_SG  G_x  kappa  seconds  nr_mask  h_fac  inf_flag
# (rates at f_a = M_pl as in the *_LvrHc_SG.dat tables; nr_mask = bit i set if level i is hydrogenic; h_fac = step
# refinement relative to the default; inf_flag = 1 if the hydrogenic + fine-structure energies put the fourth quantum
# above mu at this alpha). Rows already present with the same nr_mask and at least this h_fac are skipped, so a
# stopped task list can simply be rerun.
#   julia Compute_sg_points.jl --tasks tasks.txt --alpha_list 0.05,0.1,0.2 --out DIR --cache CACHE
#   julia Compute_sg_points.jl --precompute levels.txt --alpha_list ... --cache CACHE [--part 3/64]
# Task lines: "S1 S2 S3 S4 [rel_nmax=N] [alpha_h=c] [h_mve_fac=X] [alphas=a,b,c]" (per-task overrides).

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--tasks"
            arg_type = String
            default = ""
        "--precompute"   # level list (n-l-m per line): fill the cache for every level and alpha, then exit
            arg_type = String
            default = ""
        "--part"         # with --precompute: "i/N" does every N-th (level, alpha) pair starting at i
            arg_type = String
            default = "1/1"
        "--out"
            arg_type = String
            default = ""
        "--cache"
            arg_type = String
            default = ""
        "--alpha_list"
            arg_type = String
            default = ""
        "--rel_nmax"
            arg_type = Int
            default = 8
        "--alpha_h"
            arg_type = Float64
            default = 0.12
        "--h_mve_fac"
            arg_type = Float64
            default = 1.0
        "--ang_method"
            arg_type = String
            default = "gl"
        "--self_grav"
            arg_type = Bool
            default = true
        "--sg_recoil"
            arg_type = Bool
            default = true
        "--BHlmax"
            arg_type = Int
            default = 2
        "--alpha_min"   # only for the NptsCh mapping of Compute_all_rates.jl
            arg_type = Float64
            default = 0.03
        "--debug"
            arg_type = Bool
            default = false
    end
    return parse_args(s)
end

const PA = parse_commandline()
GF_CACHE_DIR[] = PA["cache"]
const A_SPIN = 0.95
const M_UNIT = 1.0
const NTOT_SAFE = 30000
const NPTS_BND = 2000

function have_rows(fOUT)
    rows = Tuple{Float64, Int, Float64}[]
    isfile(fOUT) || return rows
    for ln in eachline(fOUT)
        p = split(ln)
        length(p) >= 9 || continue
        push!(rows, (round(parse(Float64, p[1]), digits=10), Int(parse(Float64, p[7])), parse(Float64, p[8])))
    end
    return rows
end

function run_channel(S1, S2, S3, S4, alphas, rel_nmax, alpha_h, hfac; rpts=1500, rmaxT=100, Nang=5000000)
    a = A_SPIN; M = M_UNIT; prec = 200; cvg_acc = 1e-20
    NptsCh_Min = 40; NptsCh_Max = 70; iterC = 10; Lcheb = 4; der_acc = 1e-20
    n1, l1, m1 = parse_state_from_args(S1)
    n2, l2, m2 = parse_state_from_args(S2)
    n3, l3, m3 = parse_state_from_args(S3)
    min_m = minimum([m1, m2, m3])
    alpha_max = min(a * min_m / (2 * (1 + sqrt(1 - a^2))) * 1.03, 1.7)
    to_inf = S4 != "BH"
    nr_st = (n1 > rel_nmax, n2 > rel_nmax, n3 > rel_nmax)
    mask = Int(nr_st[1]) + 2 * Int(nr_st[2]) + 4 * Int(nr_st[3])
    rp = 1 + sqrt(1 - a^2)
    fOUT = joinpath(PA["out"], join((S1, S2, S3, S4), "_") * ".pts")
    done = have_rows(fOUT)
    for alpp in alphas
        alpp > alpha_max && continue
        mu = alpp / (M * GNew)
        maxN = maximum([n1 n2 n3])
        rmax_1 = (2.0 .^(2.0 .* 3 .- 2 .* (1 .+ 3)) .* gamma(2 .+ 2 .* 3) ./ (0.03 .^2) ./ factorial(2 .* 3 - 1) .* 7.0)
        rmax_ratio = (2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ alpp.^2 ./ factorial(big(2 .* maxN - 1)) .* 7.0) ./ rmax_1
        h_def = Float64(0.2 / rmax_ratio)
        e4 = (ergL(n1, l1, m1, mu, M, a; full=true) + ergL(n2, l2, m2, mu, M, a; full=true) - ergL(n3, l3, m3, mu, M, a; full=true)) * GNew * M
        inf_flag = real(e4) > alpp ? 1 : 0
        h_mve = h_def * hfac
        if alpha_h > 0 && (!to_inf || inf_flag == 0)
            h_mve = max(h_mve, rp * alpp / alpha_h)
        end
        h_fac = round(h_mve / h_def, sigdigits=5)
        any(r -> r[1] == round(alpp, digits=10) && r[2] == mask && r[3] >= 0.999 * h_fac, done) && continue
        NptsCh = Int(round(NptsCh_Min + (NptsCh_Max - NptsCh_Min) * clamp((alpp - PA["alpha_min"]) / (alpha_max - PA["alpha_min"]), 0, 1)))
        t0 = time()
        res = try
            gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=NPTS_BND, eps_fac=1e-3, Ntot_safe=NTOT_SAFE,
                      Nang=Nang, NON_REL=false, NR_states=nr_st, h_mve=h_mve, to_inf=to_inf, rmaxT=rmaxT, run_leaver=true, NptsCh=NptsCh,
                      cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=PA["debug"], der_acc=der_acc, Lcheb=Lcheb, use_heunc=true,
                      BHlmax=PA["BHlmax"], self_grav=PA["self_grav"], sg_recoil=PA["sg_recoil"], ang_method=Symbol(PA["ang_method"]))
        catch e
            println("FAILED ", join((S1, S2, S3, S4), " "), " alpha ", alpp, ": ", first(sprint(showerror, e), 300))
            flush(stdout)
            continue
        end
        dt = time() - t0
        g_si, g_sg, g_x = PA["self_grav"] ? (res.si, res.sg, res.cross) : (Float64(res), 0.0, 0.0)
        kappa = (g_si > 0 && g_sg > 0) ? g_x / (2 * sqrt(g_si * g_sg)) : 0.0
        open(fOUT, "a") do f
            println(f, join((alpp, g_si, g_sg, g_x, kappa, round(dt, digits=1), mask, h_fac, inf_flag), "\t"))
        end
        println("pt ", join((S1, S2, S3, S4), " "), "  alpha ", alpp, "  SI ", g_si, "  SG ", g_sg, "  kappa ", round(kappa, digits=4),
                "  ", round(dt, digits=1), " s  mask=", mask, " hfac=", h_fac)
        flush(stdout)
        GC.gc()     # many of these processes share a node; without this their heaps grow until the node runs out of memory
    end
end

default_alphas = isempty(PA["alpha_list"]) ? Float64[] : parse.(Float64, split(PA["alpha_list"], ","))

if !isempty(PA["precompute"])
    isempty(PA["cache"]) && error("--precompute needs --cache")
    i0, np = parse.(Int, split(PA["part"], "/"))
    levels = [parse_state_from_args(String(strip(ln))) for ln in eachline(PA["precompute"]) if !isempty(strip(ln)) && !startswith(strip(ln), "#")]
    jobs = [(lv, al) for lv in levels for al in default_alphas]
    for (k, ((n, l, m), alpp)) in enumerate(jobs)
        (k - i0) % np == 0 || continue
        mu = alpp / (M_UNIT * GNew)
        t0 = time()
        try
            gf_cached("erg", (n, l, m, mu, M_UNIT, A_SPIN, NTOT_SAFE)) do
                find_im_part(mu, M_UNIT, A_SPIN, n, l, m; Ntot_force=NTOT_SAFE, for_s_rates=true, return_both=true)
            end
        catch e
            println("erg FAILED ", (n, l, m), " alpha ", alpp, ": ", first(sprint(showerror, e), 200))
        end
        t1 = time()
        if n <= PA["rel_nmax"]
            try
                gf_cached("wf", (n, l, m, mu, M_UNIT, A_SPIN, NPTS_BND, NTOT_SAFE, true)) do
                    solve_radial(mu, M_UNIT, A_SPIN, n, l, m; rpts=NPTS_BND, return_erg=true, Ntot_safe=NTOT_SAFE, use_heunc=true, debug=false, pre_compute_erg=nothing)
                end
            catch e
                println("wf FAILED ", (n, l, m), " alpha ", alpp, ": ", first(sprint(showerror, e), 200))
            end
        end
        println("pre ", (n, l, m), " alpha ", alpp, "  erg ", round(t1 - t0, digits=1), " s  wf ", round(time() - t1, digits=1), " s")
        flush(stdout)
        empty!(GF_CACHE_MEM)
    end
    println("PRECOMPUTE_DONE ", PA["part"])
    exit(0)
end

isempty(PA["tasks"]) && error("give --tasks or --precompute")
isempty(PA["out"]) && error("--tasks needs --out")
mkpath(PA["out"])
for ln in eachline(PA["tasks"])
    p = split(ln)
    (isempty(p) || startswith(p[1], "#")) && continue
    rel_nmax = PA["rel_nmax"]; alpha_h = PA["alpha_h"]; hfac = PA["h_mve_fac"]; alphas = default_alphas
    for tok in p[5:end]
        k, v = split(tok, "=")
        k == "rel_nmax" && (rel_nmax = parse(Int, v))
        k == "alpha_h" && (alpha_h = parse(Float64, v))
        k == "h_mve_fac" && (hfac = parse(Float64, v))
        k == "alphas" && (alphas = parse.(Float64, split(v, ",")))
    end
    try
        run_channel(String(p[1]), String(p[2]), String(p[3]), String(p[4]), alphas, rel_nmax, alpha_h, hfac)
    catch e
        println("FAILED channel ", join(p[1:4], " "), ": ", first(sprint(showerror, e), 300))
        flush(stdout)
    end
    length(GF_CACHE_MEM) > 400 && empty!(GF_CACHE_MEM)     # the disk copies stay
end
println("TASKS_DONE ", PA["tasks"])
