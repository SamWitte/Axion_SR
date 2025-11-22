using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Dates
using Interpolations
# Include constants only if not already defined (avoid double-include when called from AxionSR module)
if !@isdefined(GNew)
    include("Core/constants.jl")
end
include("Core/evolution_helpers.jl")
include("Numerics/rate_computation.jl")
include("solve_sr_rates.jl")
include("load_rates.jl")
using Printf
using .RateComputation

function super_rad_check(M_BH, aBH, massB, f_a; tau_max=1e4, alpha_max_cut=100.0, debug=false, impose_low_cut=0.01, stop_on_a=0, eq_threshold=1e-100, abstol=1e-30, non_rel=true, high_p=true, N_pts_interp=100, N_pts_interpL=100, Nmax=3, cheby=true, spinone=false)

    alph = GNew .* M_BH .* massB

    if alph .> alpha_max_cut
        return aBH, M_BH
    elseif alph .< impose_low_cut
        return aBH, M_BH
    end

    # Use unified solve_system with spinone parameter
    if !spinone
        final_spin, final_BH = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, impose_low_cut=impose_low_cut, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel, high_p=high_p, N_pts_interp=N_pts_interp, N_pts_interpL=N_pts_interpL, Nmax=Nmax, cheby=cheby, spinone=false)
    else
        final_spin, final_BH = solve_system(massB, nothing, aBH, M_BH, tau_max, spinone=true, n_times=10000)
    end

    return final_spin, final_BH

end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu
    
    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(abs.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)) # abs just in case...
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function isapproxsigfigs(a, b, precision)
    return round(a, sigdigits=precision) == round(b, sigdigits=precision)
end

include("solve_system_unified.jl")
