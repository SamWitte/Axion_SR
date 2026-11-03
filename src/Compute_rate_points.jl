include("solve_sr_rates.jl")
include("Core/constants.jl")
include("state_utils.jl")
using ArgParse

# Compute selected alpha grid points of one scattering channel, on the same grid and with the same
# settings as Compute_all_rates.jl, optionally with hydrogenic bound states (--non_rel true: the three
# bound states from radial_bound_NR; the energies, the fourth quantum and the BH multipoles unchanged) and
# Gauss-Legendre angular overlaps (--ang_method gl).
# Writes "alpha rate seconds" for the computed points to --out; never touches rate_sve/*_LvrHc_.dat.
#   julia Compute_rate_points.jl --S1 533 --S2 533 --S3 40-6-6 --S4 BH --alpha_idx 4,7,9 --non_rel true --out rate_sve/nr_test/x.dat
#   julia Compute_rate_points.jl --tasks tasks.txt --non_rel true --ang_method gl --out rate_sve/nr_patch

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--alpha_min"
            arg_type = Float64
            default = 0.03
        "--alpha_pts"
            arg_type = Int
            default = 15
        "--alpha_idx" # comma-separated 1-based grid indices; empty = all
            arg_type = String
            default = ""
        "--S1"
            arg_type = String
            default = "211"
        "--S2"
            arg_type = String
            default = "211"
        "--S3"
            arg_type = String
            default = "722"
        "--S4"
            arg_type = String
            default = "BH"
        "--run_leaver"
            arg_type = Bool
            default = true
        "--use_heunc"
            arg_type = Bool
            default = true
        "--BHlmax"
            arg_type = Int
            default = 2
        "--non_rel"
            arg_type = Bool
            default = false
        "--ang_method" # angular overlaps: "mc" (Monte Carlo, as Compute_all_rates.jl) or "gl" (Gauss-Legendre, exact)
            arg_type = String
            default = "mc"
        "--out"
            arg_type = String
            default = ""
        "--tasks" # file with lines "S1 S2 S3 S4 i1,i2,..." computed in one process; each channel is written to
                  # <out_dir>/<S1>_<S2>_<S3>_<S4>_NRGL.dat (--out then names the directory)
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

function main(S1, S2, S3, S4, idx_str, fOUT; kpts=14, rpts=1500, rmaxT=100, Nang=5000000, Npts_Bnd=2000)
    # settings of Compute_all_rates.jl
    a = 0.95
    M = 1.0
    Ntot_safe=30000
    prec=200
    cvg_acc=1e-20
    NptsCh_Min=40
    NptsCh_Max=70
    iterC=10
    Lcheb=4
    der_acc=1e-20
    debug=true

    alpha_min, alpha_pts = parsed_args["alpha_min"], parsed_args["alpha_pts"]
    NON_REL = parsed_args["non_rel"]
    run_leaver, use_heunc, BHlmax = parsed_args["run_leaver"], parsed_args["use_heunc"], parsed_args["BHlmax"]
    ang_method = Symbol(parsed_args["ang_method"])

    n1, l1, m1 = parse_state_from_args(S1)
    n2, l2, m2 = parse_state_from_args(S2)
    n3, l3, m3 = parse_state_from_args(S3)
    min_m = minimum([m1, m2, m3])
    alpha_max = a .* min_m ./ (2 .* (1 .+ sqrt.(1 .- a.^2))) .* 1.03
    if alpha_max > 1.7
       alpha_max = 1.7
    end
    alpha_list = LinRange(alpha_min, alpha_max, alpha_pts)
    NptsCh_list = Int.(round.(LinRange(NptsCh_Min, NptsCh_Max, alpha_pts)))
    idx = isempty(idx_str) ? collect(1:alpha_pts) : parse.(Int, split(idx_str, ","))
    to_inf = S4 != "BH"

    mkpath(dirname(fOUT))
    println(S1, "\t", S2, "\t", S3, "\t", S4, "\tNON_REL = ", NON_REL, "\tang_method = ", ang_method, "\tpoints ", idx)
    for i in idx
        alpp = alpha_list[i]
        mu = alpp ./ (M * GNew)
        maxN = maximum([n1 n2 n3])
        rmax_1 = (2.0 .^(2.0 .* 3 .- 2 .* (1 .+ 3)) .* gamma(2 .+ 2 .* 3) ./ (0.03 .^2) ./ factorial(2 .* 3 - 1) .* 7.0)
        rmax_ratio = (2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ alpp.^2 ./ factorial(big(2 .* maxN - 1)) .* 7.0) ./ rmax_1
        h_mve = (0.2) ./ rmax_ratio
        t0 = time()
        rate = gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, eps_fac=1e-3, Ntot_safe=Ntot_safe, Nang=Nang, NON_REL=NON_REL, h_mve=h_mve, to_inf=to_inf, rmaxT=rmaxT, run_leaver=run_leaver, NptsCh=NptsCh_list[i], cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=debug, der_acc=der_acc, Lcheb=Lcheb, use_heunc=use_heunc, BHlmax=BHlmax, ang_method=ang_method)
        dt = time() - t0
        println("alpha \t", alpp, "\t", rate, "\t", round(dt, digits=1), " s")
        open(fOUT, "a") do f
            println(f, alpp, "\t", rate, "\t", round(dt, digits=1))
        end
    end
end

if isempty(parsed_args["tasks"])
    main(parsed_args["S1"], parsed_args["S2"], parsed_args["S3"], parsed_args["S4"], parsed_args["alpha_idx"], parsed_args["out"])
else
    for ln in eachline(parsed_args["tasks"])
        p = split(ln)
        (isempty(p) || startswith(p[1], "#")) && continue
        try
            main(String(p[1]), String(p[2]), String(p[3]), String(p[4]), length(p) > 4 ? String(p[5]) : "", joinpath(parsed_args["out"], join(p[1:4], "_") * "_NRGL.dat"))
        catch e
            println("FAILED ", join(p[1:4], " "), ": ", e)
        end
    end
end
