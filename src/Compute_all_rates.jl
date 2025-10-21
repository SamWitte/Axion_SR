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
            default = 30
            
        "--S1"
            arg_type = String
            default = "322"
            
        "--S2"
            arg_type = String
            default = "322"
            
        "--S3"
            arg_type = String
            default = "211"
            
        "--S4"
            arg_type = String
            default = "Inf"
            
        "--ftag"
            arg_type = String
            default = "_"
            
        "--run_leaver"
            arg_type = Bool
            default = false
            
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

function main(;kpts=14, rpts=1000, rmaxT=100, Nang=200000, Npts_Bnd=2000)
    a = 0.95
    M = 10.0
    Ntot_safe=20000
    NON_REL = false
    
    prec=200
    cvg_acc=5e-3
    NptsCh=100
    iterC=40
    
    debug=true
    
    
    min_m = minimum([parse(Int, S1[end]), parse(Int, S2[end]), parse(Int, S3[end])])
    alpha_max = a .* min_m ./ (2 .* (1 .+ sqrt.(1 .- a.^2))) .* 1.03
    
    # alpha_list = 10 .^LinRange(log10(alpha_min), log10(alpha_max), alpha_pts)
    alpha_list = LinRange(alpha_min, alpha_max, alpha_pts)
    # alpha_list = [0.1]
    # alpha_list = [0.45]
    output_sve = zeros(alpha_pts)
    
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
    
    fOUT = "rate_sve/"*S1*"_"*S2*"_"*S3*"_"*S4*ftag*".dat"
    overwrite_file = false
    if isfile(fOUT)&&check_err
        file_load = readdlm(fOUT)
        alpha_test = file_load[1, 1]
        computed_out = file_load[1, 2]
        mu = alpha_test ./ (M * GNew)
        h_mve = 0.1
        testVal = gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd,  eps_fac=1e-3, Ntot_safe=Ntot_safe, Nang=Nang, NON_REL=NON_REL, h_mve=h_mve, to_inf=to_inf, rmaxT=rmaxT, run_leaver=run_leaver, NptsCh=NptsCh, cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=debug)
        if (abs.(computed_out .- testVal) ./ computed_out) .> 0.4
            overwrite_file = true
            println("Test failed.... alpha, old / new \t ", alpha_test, "\t", computed_out, "\t", testVal)
            rm(fOUT)
        end
    end
    if !isfile(fOUT)||overwrite_file
        
        for i in 1:alpha_pts
            mu = alpha_list[i] ./ (M * GNew)
            
            maxN = maximum([n1 n2 n3])
            rmax_1 = (2.0 .^(2.0 .* 3 .- 2 .* (1 .+ 3)) .* gamma(2 .+ 2 .* 3) ./ (0.03.^2) ./ factorial(2 .* 3 - 1) .* 7.0)
            rmax_ratio = (2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ alpha_list[i].^2 ./ factorial(2 .* maxN - 1) .* 7.0) ./ rmax_1
            h_mve = (0.1) ./ rmax_ratio
           
            output_sve[i] = gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, eps_fac=1e-3, Ntot_safe=Ntot_safe, Nang=Nang, NON_REL=NON_REL, h_mve=h_mve, to_inf=to_inf, rmaxT=rmaxT, run_leaver=run_leaver, NptsCh=NptsCh, cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=debug)
        
            print("alpha \t", alpha_list[i], "\t", output_sve[i], "\n")
        end
       
        alpha_list = alpha_list[output_sve .> 0.0]
        output_sve = output_sve[output_sve .> 0.0]
        writedlm(fOUT, hcat(alpha_list, output_sve))
    end
    
end


main()
