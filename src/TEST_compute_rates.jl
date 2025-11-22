include("solve_sr_rates.jl")
include("Constants.jl")
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin

            
        "--alpha_min"
            arg_type = Float64
            default = 0.03
            
        "--alpha_pts"
            arg_type = Int
            default = 50
            
        "--S1"
            arg_type = String
            default = "541"
            
        "--S2"
            arg_type = String
            default = "844"
            
        "--S3"
            arg_type = String
            default = "432"
            
        "--S4"
            arg_type = String
            default = "Inf"
            
        "--ftag"
            arg_type = String
            default = "_"
  
        "--run_leaver"
            arg_type = Bool
            default = true
            

        "--check_err" # check if there could be an error in the existing file, overwrite if there is...
            arg_type = Bool
            default = false
            
            
    end
    return parse_args(s)
end


parsed_args = parse_commandline()

alpha_min = parsed_args["alpha_min"];
alpha_pts = parsed_args["alpha_pts"];
inf_nr=true


S1 = parsed_args["S1"]
S2 = parsed_args["S2"]
S3 = parsed_args["S3"]
S4 = parsed_args["S4"]

ftag = parsed_args["ftag"];

run_leaver = parsed_args["run_leaver"];
check_err = parsed_args["check_err"];

print(S1, "\t", S2, "\t", S3, "\t", S4, "\n")


function main(;kpts=14, rpts=1000, rmaxT=100, Nang=50000, Npts_Bnd=300)
    a = 0.95
    M = 1.0
    Ntot_safe=20000
    NON_REL = false
    
    prec=200
    cvg_acc=1e-8
    NptsCh=80
    iterC=20
    Lchby=4
    debug=true
    der_acc=1e-20
    
    min_m = minimum([parse(Int, S1[end]), parse(Int, S2[end]), parse(Int, S3[end])])
    alpha_max = a .* min_m ./ (2 .* (1 .+ sqrt.(1 .- a.^2))) .* 1.03
    
    alpha_list = LinRange(alpha_min, alpha_max, alpha_pts)
    # alpha_list = [0.1]

    State1 = [parse(Int, c) for c in S1]
    n1 = State1[1]
    l1 = State1[2]
    m1 = State1[3]
    
    State2 = [parse(Int, c) for c in S2]
    n2 = State2[1]
    l2 = State2[2]
    m2 = State2[3]
    
    State3 = [parse(Int, c) for c in S3]
    n3 = State3[1]
    l3 = State3[2]
    m3 = State3[3]
        
    
    if S4 == "BH"
       to_inf = false
    else
       to_inf = true
    end
    
    # alp = 0.1202191971028369
    alp = 0.26456991246737593
    mu = alp ./ (M * GNew)
            
    maxN = maximum([n1 n2 n3])
    rmax_1 = (2.0 .^(2.0 .* 3 .- 2 .* (1 .+ 3)) .* gamma(2 .+ 2 .* 3) ./ (0.03 .^2) ./ factorial(2 .* 3 - 1) .* 7.0)

    rmax_ratio = (2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ alp.^2 ./ factorial(2 .* maxN - 1) .* 7.0) ./ rmax_1
    h_mve = (0.2) ./ rmax_ratio
    # h_mve /= 10.0
    
    use_heunc=true
    output_sve = @time gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts,  rmaxT=rmaxT, run_leaver=run_leaver, NptsCh=NptsCh, cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=debug,  Npts_Bnd=Npts_Bnd, Lcheb=Lchby, der_acc=der_acc, use_heunc=use_heunc, Nang=Nang, eps_fac=1e-3)


end


main()
