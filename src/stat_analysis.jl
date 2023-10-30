using Random
using Distributions
import AffineInvariantMCMC
using DelimitedFiles
include("super_rad.jl")

Random.seed!(1234)


M_BH = 1.0
aBH = 0.98
massB = 1e-11
f_a = 1e18
tau_max=1e8
alpha_max_cut = 0.2
lg_m_low=-13
lg_m_high=-10
lg_f_high=19
lg_f_low=11



# data = open(readdlm, "data_in/TestData.dat")
data = open(readdlm, "BH_data/Doddy_Vals.dat")

function prior(theta, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    log_m, log_f = theta
    
    if (lg_m_low < log_m < lg_m_high)&&(lg_f_low < log_f < lg_f_high)
        return 0.0
    end
    return -Inf
end

function log_probability(theta, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high; tau_max=1e4, alpha_max_cut=0.2)

    lp = prior(theta, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    if !isfinite.(lp)
        return -Inf
    end

    return lp + log_likelihood(theta, data, tau_max=tau_max, alpha_max_cut=alpha_max_cut)
end

function log_likelihood(theta, data; tau_max=1e4, alpha_max_cut=0.2)
    
    log_m, log_f = theta
    MassBH_c = data[:, 1]
    MassBH_err = data[:, 2]
    SpinBH_c = data[:, 3]
    SpinBH_err = data[:, 4]
    
    num_data = length(MassBH_c)
    sum_loglike = 0.0
    for i in 1:num_data
        mass_found = false
        MassBH = nothing
        d_mass = Normal(MassBH_c[i], MassBH_err[i])
        while !mass_found
            MassBH = rand(d_mass,1)[1]
            if MassBH > 0
                mass_found = true
            end
        end
        d_spin = Normal(SpinBH_c[i], SpinBH_err[i])
        SpinBH = rand(d_spin,1)[1]
        if SpinBH > 0.998
            SpinBH = 0.998
        end
        if SpinBH < 0.0
            SpinBH = 0.0
        end
#        print(MassBH, "\t", SpinBH, "\t", 10 .^log_m ,"\n" )
        final_spin = super_rad_check(MassBH, SpinBH, 10 .^log_m, 10 .^log_f, tau_max=tau_max, alpha_max_cut=alpha_max_cut, debug=false)
#        print("Final spin \t", final_spin, "\t", SpinBH, "\t",SpinBH_err[i],   "\n")
        sum_loglike += -0.5 * (SpinBH_c[i] - final_spin).^2 / SpinBH_err[i].^2
    end
#    print("loglike \t", sum_loglike, "\n")
    return sum_loglike
end

function initialize_walkers(numwalkers, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    len_mass = abs.(lg_m_high - lg_m_low )
    len_fa = abs.(lg_f_high .- lg_f_low)
    
    x0 = rand(2, numwalkers)
    x0[1, :] .*= len_mass
    x0[1, :] .+= lg_m_low
    x0[2, :] .*= len_fa
    x0[2, :] .+= lg_f_low
    return x0
end

function mcmc_func_minimize(data; lg_m_low=-20, lg_m_high=-18, lg_f_high=19, lg_f_low=18, tau_max=1e4, alpha_max_cut=0.2)
    numdims = 2
    numwalkers = 100
    thinning = 10
    numsamples_perwalker = 10000
    burnin = 1000

    function llhood(x)
        return log_probability(x, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    end

    x0 = initialize_walkers(numwalkers, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    print("starting burn-in...\n")
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, x0, burnin, 1)
    print("starting main run...\n")
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
    flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)

    writedlm("output_mcmc/Test_mcmc.dat", flatchain')
    writedlm("output_mcmc/Test_likevals.dat", flatllhoodvals')
    
end
        
mcmc_func_minimize(data, lg_m_low=lg_m_low, lg_m_high=lg_m_high, lg_f_high=lg_f_high, lg_f_low=lg_f_low, tau_max=tau_max, alpha_max_cut=alpha_max_cut)
