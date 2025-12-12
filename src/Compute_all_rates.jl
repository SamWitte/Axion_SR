include("solve_sr_rates.jl")
include("Core/constants.jl")
include("state_utils.jl")
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
            default = true
            
        "--use_heunc"
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
use_heunc = parsed_args["use_heunc"];
check_err = parsed_args["check_err"];

print(S1, "\t", S2, "\t", S3, "\t", S4, "\n")

function main(;kpts=14, rpts=1000, rmaxT=100, Nang=5000000, Npts_Bnd=500)
    a = 0.95
    M = 1.0
    Ntot_safe=20000
    NON_REL = false
    
    prec=200
    cvg_acc=1e-10
    NptsCh_Min=40
    NptsCh_Max=70
    iterC=10
    Lcheb=4
    der_acc=1e-20
    
    debug=false
    overwrite_file = false
    
    flag_file = "fishy_rates.dat"
    
    # Parse state strings using utility functions
    n1, l1, m1 = parse_state_from_args(S1)
    n2, l2, m2 = parse_state_from_args(S2)
    n3, l3, m3 = parse_state_from_args(S3)

    min_m = minimum([m1, m2, m3])
    alpha_max = a .* min_m ./ (2 .* (1 .+ sqrt.(1 .- a.^2))) .* 1.03

    alpha_list = LinRange(alpha_min, alpha_max, alpha_pts)
    NptsCh_list = Int.(round.(LinRange(NptsCh_Min, NptsCh_Max, alpha_pts)))

    output_sve = zeros(alpha_pts)

    # Handle special case where states are equal and m values sum to > 9
    if (n1 == n2) && (l1 == l2) && (m1 == m2) && (m1 + m2 > 9)
        n3 = l1 + l2 + 1
        l3 = l1 + l2
        m3 = m1 + m2
    end

    if S4 == "BH"
        to_inf = false
    else
        to_inf = true
    end

    # Create output filename using state strings (preserves format from command line)
    fOUT = "rate_sve/"*S1*"_"*S2*"_"*S3*"_"*S4*ftag*".dat"
    if !isfile(fOUT)||overwrite_file
        
        for i in 1:alpha_pts
            alpp = alpha_list[i]
            mu = alpp ./ (M * GNew)
            
            maxN = maximum([n1 n2 n3])
            rmax_1 = (2.0 .^(2.0 .* 3 .- 2 .* (1 .+ 3)) .* gamma(2 .+ 2 .* 3) ./ (0.03 .^2) ./ factorial(2 .* 3 - 1) .* 7.0)
            
            rmax_ratio = (2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ alpp.^2 ./ factorial(big(2 .* maxN - 1)) .* 7.0) ./ rmax_1

           
            h_mve = (0.2) ./ rmax_ratio
           
            output_sve[i] = gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, eps_fac=1e-3, Ntot_safe=Ntot_safe, Nang=Nang, NON_REL=NON_REL, h_mve=h_mve, to_inf=to_inf, rmaxT=rmaxT, run_leaver=run_leaver, NptsCh=NptsCh_list[i], cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=debug, der_acc=der_acc, Lcheb=Lcheb, use_heunc=use_heunc)
            
            if (i > 5)&&(i < (alpha_pts - 4))&&(output_sve[i] > 0.0)
                itp = LinearInterpolation(log10.(alpha_list[1:i-1]), log10.(output_sve[1:i-1]), extrapolation_bc=Line())
                out_guess = 10 .^itp(log10.(alpha_list[i])) ./ output_sve[i]
                if (out_guess > 1e2)||(out_guess < 1e-2)
                    linein = S1*" "*S2*" "*S3*" "*S4
                    lines = isfile(flag_file) ? readlines(flag_file) : String[]
                    
                    if linein in lines
                        println("File already flagged....")
                    else
                        # Append to file
                        open(flag_file, "a") do f
                            write(f, linein * "\n")
                        end
                        println("Rate adddd to fishy files.")
                    end
                end
                
            end
        
            print("alpha \t", alpha_list[i], "\t", output_sve[i], "\n")
        end
       
        alpha_list = alpha_list[output_sve .> 0.0]
        output_sve = output_sve[output_sve .> 0.0]
        writedlm(fOUT, hcat(alpha_list, output_sve))
    end
    
end


main()
