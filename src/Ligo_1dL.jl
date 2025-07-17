using ArgParse
include("stat_1dL.jl")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    
        "--dataname"
            arg_type = String
            default = "LIGO_test"
    
        "--ax_mass"
            arg_type = Float64
            default = 1e-13
            
        "--fa_max"
            arg_type = Float64
            default = 1e20
            
        "--fa_min"
            arg_type = Float64
            default = 1e11

        "--tau_max_override"
            arg_type = Float64
            default = 4.5e7
            
        "--delt_M"
            arg_type = Float64
            default = 0.05
            
        "--Ftag"
            arg_type = String
            default = "_"

        "--non_rel"
            arg_type = Bool
            default = false
            
        "--Nmax"
            arg_type = Int
            default = 3
            
        "--cheby"
            arg_type = Bool
            default = false
            
        "--Nsamples"
            arg_type = Int
            default = 1
            
        "--numsamples_perwalker"
            arg_type = Int
            default = 1000
            
        "--burnin"
            arg_type = Int
            default = 500
            
        "--use_kde"
            arg_type = Bool
            default = true

    end

    return parse_args(s)
end

parsed_args = parse_commandline()

dataname = parsed_args["dataname"];
ax_mass = parsed_args["ax_mass"];
fa_min = parsed_args["fa_min"];
fa_max = parsed_args["fa_max"];

Ftag = parsed_args["Ftag"];

non_rel = parsed_args["non_rel"]
Nmax = parsed_args["Nmax"]
cheby = parsed_args["cheby"]
tau_max = parsed_args["tau_max_override"]
Nsamples = parsed_args["Nsamples"]
numsamples_perwalker = parsed_args["numsamples_perwalker"]
burnin = parsed_args["burnin"]
delt_M = parsed_args["delt_M"]
use_kde = parsed_args["use_kde"]

print("Deets...\n\n")
println("Datafile: ", dataname)
println("F name: ", Ftag)
println("Nmax: ", Nmax)
println("Use Cheby (if not, Leaver): ", cheby)
println("Nonrel? : ", non_rel)
println("Tau max : ", tau_max)
println("Use KDE : ", use_kde)

println("Ax mass : ", ax_mass)
println("Fa min : ", fa_min)
println("Fa max : ", fa_max)

#### if file exists, don't run! ...
dont_over_run = true



# Random.seed!(2084339465932781399)

data = open(readdlm, "BH_data/"*dataname*".dat")
   
Fname = "LIGO_"*dataname*"_TauMax_"*string(round(tau_max, sigdigits=2))
Fname *= "_M_ax_"*string(round(ax_mass, sigdigits=4))
Fname *= "_Nmax_$(Nmax)_"
if cheby
    Fname *= "_cheby_"
end
 
if non_rel
    Fname *= "_NonRel_"
else
    Fname *= "_FullRel_"
end

if use_kde
    Fname *= "_KDE_"
else
    Fname *= "_GA_"
end


check_exists = "output_mcmc/"*Fname*"_mcmc.dat"
if (!dont_over_run || !isfile(check_exists))
    time0=Dates.now()
    @inbounds @fastmath profileL_func_minimize(data, ax_mass, Fname, Nsamples, fa_min=fa_min, fa_max=fa_max, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, numsamples_perwalker=numsamples_perwalker, delt_M=delt_M, burnin=burnin, use_kde=use_kde)

    time1=Dates.now()
    print("\n\n Run time: ", time1-time0, "\n")
end
