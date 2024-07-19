using ArgParse
include("solve_sr_rates.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    
        "--n1"
            arg_type = Int
            default = 2
            
        "--l1"
            arg_type = Int
            default = 1
            
        "--m1"
            arg_type = Int
            default = 1
            
        "--n2"
            arg_type = Int
            default = 2
            
        "--l2"
            arg_type = Int
            default = 1
            
        "--m2"
            arg_type = Int
            default = 1
            
        "--n3"
            arg_type = Int
            default = 3
            
        "--l3"
            arg_type = Int
            default = 2
            
        "--m3"
            arg_type = Int
            default = 2
            
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

alphL = range(0.03, 1.0, 20);
rpts=50000

n1=parsed_args["n1"];
l1=parsed_args["l1"];
m1=parsed_args["m1"];
n2=parsed_args["n2"];
l2=parsed_args["l2"];
m2=parsed_args["m2"];
n3=parsed_args["n3"];
l3=parsed_args["l3"];
m3=parsed_args["m3"];

a = 0.9
M = 22.2

function runner()
    outArr = zeros(length(alphL), 3);
    for i in 1:length(alphL)
    	mu = alphL[i] ./ (GNew * M)
        outstd, outNew = test_projection_scatter(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, rmaxT=70)
        # print(outstd, "\t", outNew)
        outArr[i, :] = [alphL[i] outstd outNew]
        print(outArr[i, :], "\n")
    end

    writedlm("test_store/Xi_rel_check_Erg100_$(n1)$(l1)$(m1)_$(n2)$(l2)$(m2)_$(n3)$(l3)$(m3)_.dat", outArr)

end

runner()
