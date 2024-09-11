include("solve_sr_rates.jl")
include("Constants.jl")
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    
        "--alpha_max"
            arg_type = Float64
            default = 1.0
            
        "--alpha_min"
            arg_type = Float64
            default = 0.03
            
        "--alpha_pts"
            arg_type = Int
            default = 20
            
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
            
            
    end
    return parse_args(s)
end


parsed_args = parse_commandline()

alpha_max = parsed_args["alpha_max"];
alpha_min = parsed_args["alpha_min"];
alpha_pts = parsed_args["alpha_pts"];
inf_nr=false

S1 = parsed_args["S1"]
S2 = parsed_args["S2"]
S3 = parsed_args["S3"]
S4 = parsed_args["S4"]

function main(;kpts=14, rpts=50000, rmaxT=100, Nang=200000, Npts_Bnd=20000)
    a = 0.9
    M = 10.0
    
    alpha_list = 10 .^LinRange(log10(alpha_min), log10(alpha_max), alpha_pts)
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
    
    lF_min = 0
    
    if S4 == "BH"
        for i in 1:alpha_pts
            # print("alpha \t", alpha_list[i], "\n")
            mu = alpha_list[i] ./ (M * GNew)
            output_sve[i] = s_rate_bnd(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; kpts=kpts, rpts=rpts, rmaxT=rmaxT, inf_nr=inf_nr, Nang=Nang, Npts_Bnd=Npts_Bnd, debug=false, include_cont=true, bnd_thresh=1e-5)
            print("alpha \t", alpha_list[i], "\t", output_sve[i], "\n")
        end
    elseif S4 == "Inf"
        for i in 1:alpha_pts
            print("alpha \t", alpha_list[i], "\n")
            mu = alpha_list[i] ./ (M * GNew)
            output_sve[i] = s_rate_inf(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min; rpts=rpts, rmaxT=rmaxT,  sve_for_test=false, inf_nr=inf_nr, Npts_Bnd=Npts_Bnd, Nang=Nang)
        end
    else
        print("what is S4??? ", S4, "\n")
    end
    
    writedlm("rate_sve/"*S1*"_"*S2*"_"*S3*"_"*S4*ftag*".dat", hcat(alpha_list, output_sve))
    
end


main()
