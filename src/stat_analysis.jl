using Random
using Distributions
import AffineInvariantMCMC
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates
using MCMCDiagnosticTools
# using PyCall





function prior(theta, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    log_m, log_f = theta
    
    if (lg_m_low < log_m < lg_m_high)&&(lg_f_low < log_f < lg_f_high)
        return 0.0
    end
    return -Inf
end

function log_probability(theta, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high; tau_max=1e4, alpha_max_cut=0.2, use_input_table=true, solve_322=true, impose_low_cut=1e-100, max_mass_matrix=nothing, input_data="Masha", solve_n4=false, solve_n5=false, stop_on_a=0.0, eq_threshold=1e-100, abstol=1e-30, non_rel=true)

    lp = prior(theta, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    if !isfinite.(lp)
        print("Returning minus Inf..... \n")
        return -Inf
    end

    return lp + log_likelihood(theta, data, tau_max=tau_max, alpha_max_cut=alpha_max_cut, use_input_table=use_input_table, solve_322=solve_322, impose_low_cut=impose_low_cut, max_mass_matrix=max_mass_matrix, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel)
end

function log_likelihood(theta, data; tau_max=1e4, alpha_max_cut=0.2, use_input_table=true, solve_322=true, impose_low_cut=1e-100, max_mass_matrix=nothing, input_data="Masha", solve_n4=false, solve_n5=false, stop_on_a=0.0, eq_threshold=1e-100, abstol=1e-30, non_rel=true, debug=false)
    
    log_m, log_f = theta
    sum_loglike = 0.0
    
    if use_input_table
        MassBH_c = data[:, 1]
        MassBH_errU = data[:, 2]
        MassBH_errD = data[:, 3]
        SpinBH_c = data[:, 4]
        SpinBH_errU = data[:, 5]
        SpinBH_errD = data[:, 6]
        age = data[:, 7]
        Mc = data[:, 8]
        orbitT = data[:, 9]
        
        num_data = length(MassBH_c)
        
        for i in 1:num_data
        
            val_found = false
            p_or_neg = Int(round(rand()))
            MassBH = nothing
            d_mass = nothing
            if p_or_neg == 0
                d_mass = Normal(MassBH_c[i], MassBH_errU[i])
            else
                d_mass = Normal(MassBH_c[i], MassBH_errD[i])
            end
            
            # print("looking for mass \n")
            while !val_found
                MassBH = rand(d_mass,1)[1]
                if (p_or_neg == 0) && (MassBH >= MassBH_c[i])
                    val_found = true
                elseif (p_or_neg == 1) && (MassBH <= MassBH_c[i]) && (MassBH > 0)
                    val_found = true
                end
            end
            
            # print("looking for spin \n")
            SpinBH = nothing
            p_or_neg = Int(round(rand()))
            SpinBH = nothing
            val_found = false
            if p_or_neg == 0
                if maxSpin > SpinBH_c[i]
                    SpinBH = rand(Uniform(SpinBH_c[i], maxSpin))
                else
                    SpinBH = maxSpin
                end
            else
                d_spin = Normal(SpinBH_c[i], SpinBH_errD[i])
                while !val_found
                    SpinBH = rand(d_spin,1)[1]
                    if (SpinBH <= SpinBH_c[i]) && (SpinBH >= 0.0)
                        val_found = true
                    end
                end
            end
                

            if SpinBH > maxSpin
                SpinBH = maxSpin
            end
            if SpinBH < 0.0
                SpinBH = 0.0
            end

            if age[i] < tau_max
                maxtime = age[i]
            else
                maxtime = tau_max
            end
                
            alph = GNew .* MassBH .* 10 .^log_m #
            print("Mass and fa \t", 10 .^log_m, "\t", 10 .^log_f, "\n")
            if debug
                print(MassBH, "\t", SpinBH, "\n")
            end
            final_spin, final_mass = @time super_rad_check(MassBH, SpinBH, 10 .^log_m, 10 .^log_f, tau_max=maxtime, alpha_max_cut=alpha_max_cut, debug=false, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel)
            print("Init/Final spin \t", SpinBH, "\t", final_spin, "\n\n")
            ### Likelihood part
            
            
            # masha cut
            # a_max = 4 .* alph ./ (1 .+ 4 .* alph.^2)
            
            allwd_err = 0.01
            if (input_data == "Masha")
                a_max, MBH_max = super_rad_check(MassBH, SpinBH, 10 .^log_m, 1e19, tau_max=maxtime, alpha_max_cut=alpha_max_cut, debug=false, solve_322=false, impose_low_cut=impose_low_cut, input_data="Doddy", solve_n4=false, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel)
                
                frac_SD = (SpinBH_c[i] - a_max) * 0.2
                if ((final_spin .- a_max) .> frac_SD)
                    final_spin = SpinBH_c[i]
                end
            end
            
            if final_spin > SpinBH_c[i]
                sum_loglike += -0.5 * (SpinBH_c[i] - final_spin).^2 / SpinBH_errU[i].^2
            else
                sum_loglike += -0.5 * (SpinBH_c[i] - final_spin).^2 / SpinBH_errD[i].^2
            end
            
#            if final_mass > MassBH_c[i]
#                sum_loglike += -0.5 * (MassBH_c[i] - final_mass).^2 / MassBH_errU[i].^2
#            else
#                sum_loglike += -0.5 * (MassBH_c[i] - final_mass).^2 / MassBH_errD[i].^2
#            end
            
        end
    else
        # assume data format: [Indx DataEntry Label], Label=0 [velocity disp] Label=1 Galaxy Mass
        # Indx = data[:, 1]
        # DataEntry = data[:, 2]
        # DataLabel = data[:, 3]
        ##### LOG 10 Values ########
        MassBH_c = data[:, 1]
        MassBH_errU = data[:, 2]
        MassBH_errD = data[:, 3]
        
    
        num_data = length(MassBH_c)
        
        for i in 1:num_data
            val_found = false
            p_or_neg = Int(round(rand()))
            MassBH = nothing
            d_mass = nothing
            if p_or_neg == 0
                d_mass = Normal(MassBH_c[i], MassBH_errU[i])
            else
                d_mass = Normal(MassBH_c[i], MassBH_errD[i])
            end
            
            while !val_found
                MassBH = 10 .^rand(d_mass,1)[1]
                if (p_or_neg == 0) && (MassBH >= 10 .^MassBH_c[i])
                    val_found = true
                elseif (p_or_neg == 1) && (MassBH <= 10 .^MassBH_c[i]) && (MassBH > 0)
                    val_found = true
                end
            end
            
            prior_spins(a) = ones(size(a))  # agnostic spin prior.
            min_spin = one_d_spin_fixed_mass(MassBH .* Ms, prior_spins, max_mass_matrix; return_all=false)
                                
            spinBH_sample = rand() .* (maxSpin .- min_spin) .+ min_spin
            
            # alph = GNew .* MassBH .* 10 .^log_m #
            # day_to_inVeV = 24.0 * 60 * 60 / 6.58e-16

            final_spin, final_mass = super_rad_check(MassBH, spinBH_sample, 10 .^log_m, 10 .^log_f, tau_max=tau_max, alpha_max_cut=alpha_max_cut, debug=false, solve_322=solve_322, impose_low_cut=alpha_min_cut, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, abstol=abstol, non_rel=non_rel)
            loglike = tde_like(final_mass, final_spin, max_mass_matrix; plot=false)

            sum_loglike += loglike
            
        end
        
    
    end
    return sum_loglike
end

function spin_prior(x)
    return 1.0
end

function mass_prior(x, log10_M, sig)
    return exp.(-(x .- log10_M).^2 ./ (2 .* sig.^2)) ./ (sqrt.(2 .* pi) .* sig)
end

function BH_mass_from_obs(val, label)
    if label == 1
        MBH = 10 .^ (7.43 .+ 1.61 .+ log10.(val ./ 3e10))
        return MBH
    else
        MBH = 10 .^ (7.87 .+ 4.384 .* log10.(val ./ 160.0))
        return MBH
    end
end

function initialize_walkers(numwalkers, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    len_mass = abs.(lg_m_high - lg_m_low )
    len_fa = abs.(lg_f_high .- lg_f_low)
    mu_std = log10.(0.2 ./ (data[1, 1] .* GNew))
    
    x0 = rand(2, numwalkers)
    x0[1, :] .*= len_mass
    x0[1, :] .+= lg_m_low
    # x0[1, :] .= mu_std
    x0[2, :] .*= len_fa
    x0[2, :] .+= lg_f_low
    return x0
end

function mcmc_func_minimize(data, Fname; lg_m_low=-20, lg_m_high=-18, lg_f_high=19, lg_f_low=18, tau_max=1e4, alpha_max_cut=0.2, alpha_min_cut=1e-100, use_input_table=true, solve_322=true, numwalkers=10, thinning=1, numsamples_perwalker=2000, burnin=500, max_mass_matrix=nothing, input_data="Masha", solve_n4=false, solve_n5=false, stop_on_a=0.0, eq_threshold=1e-100, abstol=1e-30, non_rel=true)
    numdims = 2


    function llhood(x)
        return log_probability(x, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high, tau_max=tau_max, alpha_max_cut=alpha_max_cut, use_input_table=use_input_table, solve_322=solve_322, impose_low_cut=alpha_min_cut, max_mass_matrix=max_mass_matrix, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel)
    end

    # Random.seed!(1215740967981745887)
    x0 = initialize_walkers(numwalkers, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    # x0 = [-12.11097910149243 -10.946032301911332 -11.186738089283743 -11.15637998808181; 13.015578527313194 17.574769391364217 15.064923984685148 13.622376649826602]
    print("Init walkers \t", x0, "\n")
    
    print("Starting burn-in...\n")
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, x0, burnin, 1)
    print("Starting main run...\n")
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
    
    print("Finished main run...\n")
    flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)
    rHvals = rhat(flatchain')
    essVals = ess(flatchain')
    
    print("Convergence... (rhat - 1) \t ", rHvals - 1.0, "\n")
    print("Convergence... ess \t ", essVals, "\n")
    
    writedlm("output_mcmc/"*Fname*"_mcmc.dat", flatchain')
#    writedlm("output_mcmc/"*Fname*"_likevals.dat", flatllhoodvals')
    
end
   
