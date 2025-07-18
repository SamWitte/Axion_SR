using Glob
using DelimitedFiles
using ArgParse
using KernelDensity, StatsBase
using QuadGK
using Interpolations

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
            
        "--Nwalkers"
            arg_type = Int
            default = 4
            
        "--targetP"
            arg_type = Float64
            default = 0.95
            
        "--use_kde"
            arg_type = Bool
            default = true
            
        "--kde_test"
            arg_type = Bool
            default = false

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
Nwalkers = parsed_args["Nwalkers"]

delt_M = parsed_args["delt_M"]

targetP = parsed_args["targetP"]
kde_test = parsed_args["kde_test"]
use_kde = parsed_args["use_kde"]


Fname = "LIGO_"*dataname*"_TauMax_"*string(round(tau_max, sigdigits=2))
Fname *= "_M_ax_*"
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

fileList = Glob.glob("output_mcmc/"*Fname*"_mcmc.dat")


output_lim = zeros(length(fileList), 2)
for i in 1:length(fileList)
    # grab mass
    pattern = r"M_ax_(-?\d+(\.\d+)?([eE][+-]?\d+)?)"
    max = match(pattern, fileList[i])
    value_str = max.captures[1]
    value = parse(Float64, value_str)
    
    # get limit
    input_data = vec(readdlm(fileList[i]))
    
    # Compute KDE
    kd = kde(input_data)
    # Create linear interpolation from the KDE
    itp = interpolate(kd.density, BSpline(Linear()))
    x_itp = range(first(kd.x), last(kd.x), length=length(kd.x))
    f_interp = extrapolate(scale(itp, x_itp), Flat())  # safely extrapolate beyond bounds

    # kd.x = points where density is evaluated
    # kd.density = density values
    function kde_cdf(x)
        integral, _ = quadgk(t -> f_interp(t), minimum(kd.x), x)
        return integral
    end
    
    if kde_test
        faL = LinRange(log10.(fa_min), log10.(fa_max), 50)
        writedlm("output_mcmc/KDE_"*string(i)*".dat", cat(faL, f_interp(faL), dims=2))
    end
    
    f_low = log10.(fa_min)
    f_high = log10.(fa_max)
    f_test = mean([f_low, f_high])
    found_it = false
    npt=0
    while !found_it
        cdfV = kde_cdf(f_test)
        if abs.(cdfV - targetP) < 0.001
            found_it = true
        else
            if cdfV > targetP
                f_high = f_test
                f_test = mean([f_low, f_high])
            else
                f_low = f_test
                f_test = mean([f_low, f_high])
            end
            npt += 1
            
            if npt > 50
                println("Fail \t ", f_low, " ", f_high, " ", f_test, " ", cdfV)
                found_it = true
            end
        end
    
    end
    output_lim[i,:] = [value, f_test]

end

Fname = "LIGO_"*dataname*"_TauMax_"*string(round(tau_max, sigdigits=2))
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

sorted_out = sort(eachrow(output_lim), by = x -> x[1]) |> collect
sorted_out = vcat(sorted_out...)  # convert back to a Matrix
writedlm("output_mcmc/Lim_"*Fname*".dat", sorted_out)
