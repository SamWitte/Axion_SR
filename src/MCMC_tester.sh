alpha_max_cut=10.0 
alpha_min_cut=0.001 
tau_max_override=5.0e7 
Ftag="_NR_v1_"
one_BH="Cyg_X1"
numwalkers=6
numsamples_perwalker=2000
burnin=400
stop_on_a=0.9
non_rel=true
Nmax=3
cheby=false
julia MCMC.jl --alpha_max_cut $alpha_max_cut --alpha_min_cut $alpha_min_cut --tau_max_override $tau_max_override --Ftag $Ftag --one_BH $one_BH --numwalkers $numwalkers --numsamples_perwalker $numsamples_perwalker --stop_on_a $stop_on_a --burnin $burnin --non_rel $non_rel --Nmax $Nmax --cheby $cheby
