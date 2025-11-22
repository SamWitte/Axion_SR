dataname="GW231123_NRSur"
tau_max_override=1e7
delt_M=0.05
Nmax=5
Ftag="_"
targetP=0.9
julia Ligo_limit.jl --dataname $dataname --tau_max_override $tau_max_override --delt_M $delt_M --Nmax $Nmax  --Ftag $Ftag --kde_test false --use_kde true --targetP $targetP
