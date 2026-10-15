SG rate tables used for the corrected 5-level evolution (same format as rate_sve/*_LvrHc_SG.dat:
alpha, Gamma_SI, Gamma_SG, Gamma_x, kappa).

* 211_211_322_BH, 322_322_211_Inf, 322_322_544_BH, 544_544_322_Inf: exact one-graviton exchange
  (= ../tables/*_LvrHc_SGrel.dat).
* the other 8 dipole channels of the network: Newtonian kernel + BH recoil
  (Compute_all_rates.jl --self_grav true --sg_recoil true --ftag _LvrHc_ --alpha_min 0.03 --alpha_pts 14 --ang_method gl).
* all other channels: unchanged Newtonian tables from rate_sve/ (no dipole, recoil has no effect).

Note: the file names match the `*_LvrHc_SG.dat` pattern in .gitignore; use `git add -f` to commit them.
