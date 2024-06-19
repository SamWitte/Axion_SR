
function load_rate_coeffs(mu, M, a, f_a; non_rel=true, input_data="Me", solve_n4=true)
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4
    
    Drate = Dict()
    
    Drate["211_322^GW"] = 0.0 * alph^16
    if non_rel
        
        if input_data == "Doddy"
            Drate["211_211^322^BH"] = 0.0 # format gamma_{A x B}^{C x D}
            Drate["322_322^211^Inf"] = 0.0
            Drate["211_211^GW"] = 0.0
            Drate["322_211^GW"] = 0.0
            Drate["211_211_211^Inf"] = 0.0
            Drate["322_322^GW"] = 0.0
            
            # n = 4
            if solve_n4
                Drate["211_411^322^BH"] = 0.0
                Drate["411_411^322^BH"] = 0.0
                Drate["322_411^211^Inf"] = 0.0
                Drate["211_211^422^BH"] = 0.0
                Drate["411_422^211^Inf"] = 0.0
                Drate["411_411^422^BH"] = 0.0
                Drate["211_422^433^BH"] = 0.0
                Drate["422_322^211^Inf"] = 0.0
                Drate["433_433^211^Inf"] = 0.0
                Drate["322_433^211^Inf"] = 0.0
                Drate["211_322^433^BH"] = 0.0
            end
            
            
        else
            Drate["211_211^322^BH"] = 4.0e-7 .* alph^11 .* faFac * rP
            Drate["322_322^211^Inf"] = 1.0e-8 * alph^8 .* faFac
            Drate["211_211^GW"] = 1.0e-2 * alph^14
            Drate["322_211^GW"] = 5.0e-6 * alph^10
            Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
            Drate["322_322^GW"] = 3.0e-8 * alph^18
            
            # n = 4
            if solve_n4
                Drate["211_411^322^BH"] = 2.5e-8 * alph^11 * faFac * rP
                Drate["411_411^322^BH"] = 9.8e-11 * alph^11 * faFac * rP
                Drate["322_411^211^Inf"] = 2.8e-9 * alph^8 * faFac
                Drate["211_211^422^BH"] = 1.5e-7 * alph^11 * faFac * rP
                Drate["411_422^211^Inf"] = 2.2e-9 * alph^8 * faFac
                Drate["411_411^422^BH"] = 2.3e-7 * alph^7 * faFac * rP
                Drate["211_422^433^BH"] = 1.1e-9 * alph^7 * faFac * rP
                Drate["422_322^211^Inf"] = 1.6e-8 * alph^8 * faFac
                Drate["433_433^211^Inf"] = 9.2e-11 * alph^8 * faFac
                Drate["322_433^211^Inf"] = 2.6e-9 * alph^8 * faFac
                Drate["211_322^433^BH"] = 9.1e-8 * alph^11 * faFac * rP
            end
            
        end
        
    else
        Drate["211_211^GW"] = 1.0e-2 * alph^14
        Drate["322_211^GW"] = 5.0e-6 * alph^10
        Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        Drate["322_322^GW"] = 3.0e-8 * alph^18
        
        
        # rates computed for fixed rP(a=0.9)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.9^2))
        
        # load files and get rates
        data = open(readdlm, "rate_sve/211_211_322_BH.dat")
        itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
        Drate["211_211^322^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
        
        data = open(readdlm, "rate_sve/322_322_211_Inf.dat")
        itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
        Drate["322_322^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
        
        if solve_n4
            data = open(readdlm, "rate_sve/211_411_322_BH.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_411^322^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            
            data = open(readdlm, "rate_sve/411_411_322_BH.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["411_411^322^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
        
            
            data = open(readdlm, "rate_sve/322_411_211_Inf.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["322_411^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/211_211_422_BH.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_211^422^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/411_422_211_Inf.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["411_422^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/411_411_422.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["411_411^422^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/211_422_433_BH.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_422^433^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/422_322_211_Inf.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["422_322^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/433_433_211_Inf.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["433_433^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/322_433_211_Inf.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["322_433^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/211_322_433_BH.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_322^433^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac

        end
        
    end
    
    return Drate
end
