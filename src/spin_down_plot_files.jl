include("super_rad.jl")

fa_vals = [1e18, 1e16, 1e14, 1e13, 1e12]
ma_list = [2e-12, 2.5e-12, 3e-12, 3.5e-12, 4e-12, 5e-12]
MassBH = 22.2
SpinBH = 0.99
tau_max = 4.8e6
impose_low_cut = 1e-100
high_p=true
N_pts_interp=30
N_pts_interpL=30
n_times = 1000000


for i in 1:length(fa_vals)
    out_m2 = zeros(length(ma_list), 4)
    out_all = zeros(length(ma_list), 4)

    for j in 1:length(ma_list)
        alph = GNew * MassBH * ma_list[j]
        print("\n")
        print(alph, "\t", fa_vals[i], "\n\n")
        print("N = 3 \n")
        final_spin_3, final_BH = @time solve_system(ma_list[j], fa_vals[i], SpinBH, MassBH, tau_max, debug=false, input_data="Me", solve_n4=false, solve_n5=false, non_rel=false, max_m_2=true, high_p=high_p, N_pts_interp=N_pts_interp, n_times=n_times, N_pts_interpL=N_pts_interpL)
        print("N = 3A \n")
        final_spin_3A, final_BH = @time solve_system(ma_list[j], fa_vals[i], SpinBH, MassBH, tau_max, debug=false, input_data="Me", solve_n4=false, solve_n5=false, non_rel=false, max_m_2=false, high_p=high_p, N_pts_interp=N_pts_interp, n_times=n_times, N_pts_interpL=N_pts_interpL)
        
        print("N = 4 \n")
        final_spin_4, final_BH = @time solve_system(ma_list[j], fa_vals[i], SpinBH, MassBH, tau_max, debug=false, input_data="Me", solve_n4=true, solve_n5=false, non_rel=false, max_m_2=true, high_p=high_p, N_pts_interp=N_pts_interp, n_times=n_times, N_pts_interpL=N_pts_interpL)
        print("N = 4A \n")
        final_spin_4A, final_BH = @time solve_system(ma_list[j], fa_vals[i], SpinBH, MassBH, tau_max, debug=false, input_data="Me", solve_n4=true, solve_n5=false, non_rel=false, max_m_2=false, high_p=high_p, N_pts_interp=N_pts_interp, n_times=n_times, N_pts_interpL=N_pts_interpL)
        
        print("N = 5 \n")
        final_spin_5, final_BH = @time solve_system(ma_list[j], fa_vals[i], SpinBH, MassBH, tau_max, debug=false, input_data="Me", solve_n4=true, solve_n5=true, non_rel=false, max_m_2=true,  impose_low_cut=impose_low_cut, high_p=high_p, N_pts_interp=N_pts_interp, n_times=n_times, N_pts_interpL=N_pts_interpL)
        print("N = 5A \n")
        final_spin_5A, final_BH = @time solve_system(ma_list[j], fa_vals[i], SpinBH, MassBH, tau_max, debug=false, input_data="Me", solve_n4=true, solve_n5=true, non_rel=false, max_m_2=false, impose_low_cut=impose_low_cut, high_p=high_p, N_pts_interp=N_pts_interp, n_times=n_times, N_pts_interpL=N_pts_interpL)
        
        out_m2[j, :] .= [alph, final_spin_3, final_spin_4, final_spin_5]
        out_all[j, :] .= [alph, final_spin_3A, final_spin_4A, final_spin_5A]
        
    end
    
    writedlm("test_store/a_final_info_$(fa_vals[i])_SpinBH_$(SpinBH)_tauM_$(tau_max)_m2SD_.dat", out_m2)
    writedlm("test_store/a_final_info_$(fa_vals[i])_SpinBH_$(SpinBH)_tauM_$(tau_max)_.dat", out_all)
end

