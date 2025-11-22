include("solve_sr_rates.jl")
using DelimitedFiles

N_mu = 5000
Mbh_range = 10 .^ LinRange(log10.(3), log10.(50e9), 5000)
mu_range = 10 .^ LinRange(-22, -10, N_mu)

a_spin = 0.99


del_thresh = 0.01
tscale = 5e7
ftag = "_tau_5e7_"
outarr = zeros(N_mu, 2)


for i in 1:N_mu

    fa_guess = 1e17
    fa_max = 1e22
    fa_min = 1e9

    
    idx_bh = 0
    FoM = 0.0
    
    for j in 1:length(Mbh_range)
        alp = GNew .* Mbh_range[j] .* mu_range[i]
        eq_211 = 2.5e-1 * (0.01 ./ alp).^3 .* sqrt.(a_spin ./ 0.9) .* (fa_guess ./ 1e22).^2
        ahold = 4 .* alp .* (2.0 - 1.0) ./ ((2.0 .- 1.0).^2 .+ 4 .* alp.^2)
        tau = 1.0 ./ (2 .*  sr_rates(2, 1, 1, mu_range[i], Mbh_range[j], a_spin) .* alp)  .* 6.58e-16 ./ 3.15e7
        if tau < 0
            tau = 1e100
        end
        # println(alp, "\t", Mbh_range[j], "\t", ahold, "\t ", 180 .* tau)
        if (ahold < a_spin)&&(180 .* tau .< tscale)
            if (eq_211 ./ tau) .> FoM
                FoM = (eq_211 ./ tau)
                idx_bh = j
            end
        end
    end
    
    if idx_bh == 0
        outarr[i, :] .= [mu_range[i], 1e22]
        continue
    end
    
    found_f = false
    cnt = 0
    while !found_f
        
        alp = GNew .* Mbh_range[idx_bh] .* mu_range[i]
        eq_211 = 2.5e-1 * (0.01 ./ alp).^3 .* sqrt.(a_spin ./ 0.9) .* (fa_guess ./ 1e15).^2
        tau = 1.0 ./ (2 .*  sr_rates(2, 1, 1, alp ./ (GNew .* Mbh_range[idx_bh]),Mbh_range[idx_bh], a_spin) .* alp)  .* 6.58e-16 ./ 3.15e7
        if tau < 0
            tau = 1e100
        end
        ahold = 4 .* alp .* (2.0 - 1.0) ./ ((2.0 .- 1.0).^2 .+ 4 .* alp.^2)
        g_time = log.(eq_211 ./ 1e-80) .* tau
        sd_time = (a_spin .- ahold) ./ (eq_211 ./ tau)
        tot_time = g_time .+ sd_time
        if abs.(tscale .- tot_time) ./ tscale .< del_thresh
            found_f = true
            continue
        end
        
        if tot_time .< tscale
            fa_max = fa_guess
            fa_guess = 10 .^ ( (log10.(fa_guess) .+ log10.(fa_min)) ./ 2)
        else
            fa_min = fa_guess
            fa_guess = 10 .^ ( (log10.(fa_guess) .+ log10.(fa_max)) ./ 2)
        end
        
        cnt +=1
        if cnt > 100
            println(cnt, "\t", fa_guess, "\t", tot_time)
            found_f = true
        end
    end
    
    outarr[i, :] .= [mu_range[i], fa_guess]
end

writedlm("test_store/Projected_Eg_"*ftag*".dat", outarr)
