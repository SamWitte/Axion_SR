
function load_rate_coeffs(mu, M, a, f_a; non_rel=true, input_data="Me")
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4
    
    Drate()
    
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
            
            
        else
            Drate["211_211^322^BH"] = 4.0e-7 .* alph^11 .* faFac * rP
            Drate["322_322^211^Inf"] = 1.0e-8 * alph^8 .* faFac
            Drate["211_211^GW"] = 1.0e-2 * alph^14
            Drate["322_211^GW"] = 5.0e-6 * alph^10
            Drate["211_211_211^Inf"] = 1.5e-8 * alph^21
            Drate["322_322^GW"] = 3.0e-8 * alph^18
            
            # n = 4
            Drate["211_411^322^BH"] = 2.5e-8 * alph^11 * faFac * rP
            Drate["411_411^322^BH"] = 9.8e-11 * alph^11 * faFac * rP
            Drate["322_411^211^Inf"] = 2.8e-9 * alph^8 * faFac
            Drate["211_211^422^BH"] = 1.5e-7 * alph^11 * faFac * rP
            Drate["411_422^211^Inf"] = 5.0e-12 * alph^8 * faFac
            Drate["411_411^422^BH"] = 2.3e-7 * alph^7 * faFac * rP
            Drate["211_422^433^BH"] = 1.1e-9 * alph^7 * faFac * rP
            Drate["422_322^211^Inf"] = 4.0e-9 * alph^8 * faFac
            Drate["433_433^211^Inf"] = 9.2e-11 * alph^8 * faFac
            Drate["322_433^211^Inf"] = 2.6e-9 * alph^8 * faFac
            Drate["211_322^433^BH"] = 9.1e-8 * alph^11 * faFac
            
        end
        
    else
        # rates computed for fixed rP(a=0.9)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.9^2))
        
        # load files
        
        
        # interpolate
        
        
        print("Not done yet... \n")
    end
    
    return Drate
end
