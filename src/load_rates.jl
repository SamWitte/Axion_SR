using Glob

function load_rate_coeffs(mu, M, a, f_a, Nmax, SR_rates; non_rel=true)
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4
    
    rate_list = readdlm("rate_sve/load_rate_input_Nmax_$(Nmax).txt")
    cnt = 1
    
    kill_lvls = []
    for nn in 1:Nmax, l in 1:(nn - 1),  m in 1:l
        if SR_rates[cnt] <= 0
            push!(kill_lvls, [nn, l, m])
        end
    end
    keep_idx = []
    for i in 1:length(rate_list[:,1])
        if !((rate_list[i, 1] in kill_lvls)||(rate_list[i, 2] in kill_lvls)||(rate_list[i, 3] in kill_lvls))
            push!(keep_idx, i)
        end
    end
    rate_list = rate_list[keep_idx, :]
    
    Drate = Dict()
    
    Drate["211_322^GW"] = 0.0 * alph^16
    
    
    ### EDITING ##
    
    if non_rel
        include_m1 = true
        include_m2 = true
        include_m3 = true
        if include_m1
            Drate["211_211^322^BH"] = 4.2e-7 .* alph^11 .* faFac * rP
            Drate["322_322^211^Inf"] = 1.1e-8 * alph^8 .* faFac
            Drate["211_211^GW"] = 1.0e-2 * alph^14
            Drate["322_211^GW"] = 5.0e-6 * alph^10
            Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
            
            Drate["211_311^322^BH"] = 3.1e-10 .* alph^7 .* faFac * rP
            Drate["311_311^211^Inf"] = 5.1e-8 .* alph^8 * faFac
            Drate["311_322^211^Inf"] = 1.2e-8 .* alph^8 * faFac
            Drate["311_311^322^BH"] = 1.62e-10 .* alph^7 .* faFac * rP
            
        end
        Drate["322_322^GW"] = 3.0e-8 * alph^18
            
        # n = 4
        if Nmax >= 4
            if include_m1
                Drate["211_411^322^BH"] = 2.5e-8 * alph^11 * faFac * rP
                Drate["322_411^211^Inf"] = 3.7e-9 * alph^8 * faFac
                Drate["211_211^422^BH"] = 1.5e-7 * alph^11 * faFac * rP
                Drate["411_422^211^Inf"] = 2.2e-9 * alph^8 * faFac
                Drate["411_411^322^BH"] = 1.7e-11 * alph^11 * faFac * rP ### Disagree
                Drate["411_411^422^BH"] = 2.2e-11 * alph^7 * faFac * rP
                Drate["211_422^433^BH"] = 7.83e-11 * alph^7 * faFac * rP ### Disagree
                
                # New ones
                Drate["411_411^211^Inf"] = 1.7e-9 * alph^8 * faFac
                Drate["411_433^211^Inf"] = 1.1e-10 * alph^8 * faFac
                Drate["422_422^211^Inf"] = 1.6e-9 * alph^8 * faFac
                Drate["422_433^211^Inf"] = 6.1e-10 * alph^8 * faFac
                Drate["211_411^422^BH"] = 3.2e-11 * alph^7 * faFac * rP
                Drate["411_422^433^BH"] = 2.3e-11 * alph^7 * faFac * rP
                
                Drate["211_311^422^BH"] = 2.7e-7 .* alph^11 .* faFac * rP
                Drate["311_311^422^BH"] =  1.7e-11 .* alph^11 .* faFac * rP
                Drate["311_322^433^BH"] = 7.0e-8 .* alph^11 .* faFac * rP
                Drate["311_411^211^Inf"] = 1.9e-8 .* alph^8 * faFac
                Drate["311_422^211^Inf"] = 7.0e-9 .* alph^8 * faFac
                Drate["311_433^211^Inf"] = 2.2e-10 .* alph^8 * faFac
                Drate["311_411^322^BH"] = 1.9e-10 .* alph^7 .* faFac * rP
                Drate["311_411^422^BH"] = 3.8e-13 .* alph^7 .* faFac * rP
                Drate["311_422^433^BH"] = 7.7e-12 .* alph^7 .* faFac * rP
                
                
                    
                Drate["422_322^211^Inf"] = 1.6e-8 * alph^8 * faFac
                Drate["433_433^211^Inf"] = 9.2e-11 * alph^8 * faFac
                Drate["322_433^211^Inf"] = 2.6e-9 * alph^8 * faFac
                Drate["211_322^433^BH"] = 9.1e-8 * alph^11 * faFac * rP
                Drate["322_411^433^BH"] = 3.8e-11 * alph^7 * faFac * rP
            end
           
                
            if Nmax >= 5
                if include_m1
                    Drate["211_211^522^BH"] = 7.5e-8 * alph^11 * faFac * rP
                    Drate["322_411^533^BH"] = 2.0e-8 * alph^11 * faFac * rP
                    Drate["411_411^522^BH"] = 9.0e-11 * alph^11 * faFac * rP #
                    
                    Drate["211_311^522^BH"] = 1.0e-7 .* alph^11 .* faFac * rP
                    Drate["211_411^522^BH"] = 9.9e-8 .* alph^11 .* faFac * rP
                    Drate["211_322^533^BH"] = 3.1e-8 .* alph^11 .* faFac * rP
                    Drate["211_422^533^BH"] = 1.1e-7 .* alph^11 .* faFac * rP
                    Drate["211_433^544^BH"] = 1.1e-9 .* alph^11 .* faFac * rP
                    Drate["211_511^322^BH"] = 2.9e-8 .* alph^11 .* faFac * rP
                    Drate["211_511^422^BH"] = 2.1e-10 .* alph^11 .* faFac * rP
                    Drate["211_522^433^BH"] = 6.5e-8 .* alph^11 .* faFac * rP
                    Drate["211_511^522^BH"] = 6.6e-12 .* alph^7 .* faFac * rP
                    Drate["211_522^533^BH"] = 2.6e-11 .* alph^7 .* faFac * rP
                    Drate["211_533^544^BH"] = 4.6e-13 .* alph^7 .* faFac * rP
                    Drate["311_311^522^BH"] = 2.6e-12 .* alph^11 .* faFac * rP
                    Drate["311_322^533^BH"] = 2.9e-8 .* alph^11 .* faFac * rP
                    Drate["311_411^522^BH"] = 8.1e-10 .* alph^11 .* faFac * rP
                    Drate["311_422^533^BH"] = 2.1e-9 .* alph^11 .* faFac * rP
                    Drate["311_433^544^BH"] = 1.5e-8 .* alph^11 .* faFac * rP
                    Drate["311_522^211^Inf"] = 3.9e-9 * alph^8 * faFac
                    Drate["311_533^211^Inf"] = 7.6e-11 * alph^8 * faFac
                    Drate["311_544^211^Inf"] = 5.3e-11 * alph^8 * faFac
                    Drate["322_511^211^Inf"] = 1.7e-9 * alph^8 * faFac
                    Drate["411_511^211^Inf"] = 4.3e-14 * alph^8 * faFac
                    Drate["411_522^211^Inf"] = 8.1e-10 * alph^8 * faFac
                    Drate["411_533^211^Inf"] = 5.8e-16 * alph^8 * faFac
                    Drate["411_544^211^Inf"] = 4.2e-14 * alph^8 * faFac
                    Drate["422_511^211^Inf"] = 2.0e-10 * alph^8 * faFac
                    Drate["422_522^211^Inf"] = 3.7e-9 * alph^8 * faFac
                    Drate["422_533^211^Inf"] = 6.2e-14 * alph^8 * faFac
                    Drate["422_544^211^Inf"] = 2.7e-12 * alph^8 * faFac
                    Drate["433_511^211^Inf"] = 2.6e-16 * alph^8 * faFac
                    Drate["433_522^211^Inf"] = 2.7e-10 * alph^8 * faFac
                    Drate["433_533^211^Inf"] = 1.9e-10 * alph^8 * faFac
                    Drate["433_544^211^Inf"] = 2.4e-11 * alph^8 * faFac
                    
                    Drate["511_511^422^BH"] = 4.0e-13 .* alph^11 .* faFac * rP
                    Drate["511_522^433^BH"] = 1.6e-12 .* alph^11 .* faFac * rP
                    Drate["511_511^522^BH"] = 4.1e-12 .* alph^7 .* faFac * rP
                    Drate["511_522^533^BH"] = 1.5e-12 .* alph^7 .* faFac * rP
                    Drate["511_533^544^BH"] = 5.2e-13 .* alph^7 .* faFac * rP
                    
                    
                end
                if include_m2
                    Drate["322_322^544^BH"] = 1.9e-9 * alph^11 * faFac * rP
                    Drate["322_422^544^BH"] = 1.3e-11 * alph^11 * faFac * rP #
                    Drate["322_522^544^BH"] = 3.4e-12 * alph^7 * faFac * rP
                    Drate["422_422^544^BH"] = 2.3e-9 * alph^11 * faFac * rP #
                    Drate["422_522^544^BH"] = 3.7e-14 * alph^7 * faFac * rP
                    Drate["522_522^544^BH"] = 2.7e-13 * alph^7 * faFac * rP #
                    
                    Drate["422_544^322^Inf"] = 2.5e-11 * alph^8 * faFac
                    Drate["433_544^322^Inf"] = 7.8e-10 * alph^8 * faFac
                    Drate["422_533^322^Inf"] = 1.2e-9 * alph^8 * faFac
                    Drate["433_533^322^Inf"] = 2.8e-9 * alph^8 * faFac
                    Drate["433_522^322^Inf"] = 6.3e-10 * alph^8 * faFac
                    Drate["422_522^322^Inf"] = 1.6e-9 * alph^8 * faFac
                    
                    Drate["522_544^322^Inf"] = 2.2e-11 * alph^8 * faFac
                    Drate["533_544^322^Inf"] = 1.8e-10 * alph^8 * faFac
                    Drate["544_544^322^Inf"] = 4.3e-11 * alph^8 * faFac
                    Drate["522_533^322^Inf"] = 4.4e-10 * alph^8 * faFac
                    Drate["533_533^322^Inf"] = 3.1e-10 * alph^8 * faFac
                    Drate["522_522^322^Inf"] = 1.6e-10 * alph^8 * faFac
                end
            end
        end
            
    elseif !non_rel
        
        Drate["211_211^GW"] = 1.0e-2 * alph^14
        Drate["322_211^GW"] = 5.0e-6 * alph^10
        # Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        Drate["322_322^GW"] = 3.0e-8 * alph^18
        
        
        # rates computed for fixed rP(a=0.9)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.9^2))
        
        dirN = "rate_sve/"
        ftag = "_LvrHc_"
        
        for i in 1:length(rate_list[:,1])
            nm_tag = string(rate_list[i, 1]) * "_" * string(rate_list[i, 2]) * "^" * string(rate_list[i, 3]) * "^" * string(rate_list[i, 4])
            fileT = dirN * string(rate_list[i, 1]) * "_" * string(rate_list[i, 2]) * "_" * string(rate_list[i, 3]) * "_" * string(rate_list[i, 4]) * ftag * ".dat"
            if isfile(fileT)
                data = open(readdlm, fileT)
                data = data[data[:,2] .!= 0.0, :]
                
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                rate_out = 10 .^itp(log10.(alph)) .* faFac
                if string(rate_list[i, 4]) == "BH"
                    rate_out *= rP_ratio
                end
                Drate[nm_tag] = rate_out
            end
        end
    end

    return Drate
end

function key_to_indx(keyN, Nmax)
    if keyN[9:10] == "GW"
        totN = 3
    else
        totN = 4
    end
    outPix = zeros(Int, totN)
    sgn = zeros(totN)
    
    for i in 1:totN
        if i == 1
            nme = keyN[1:3]
            sgn[i] = -1.0
        elseif i == 2
            nme = keyN[5:7]
            sgn[i] = -1.0
        elseif i == 3
            if totN == 4
                nme = keyN[9:11]
                if keyN[8] == "_"
                    sgn[i] = -1.0
                else
                    sgn[i] = 1.0
                end
            else
                nme = "GW"
                sgn[i] = 1.0
            end
        else
            nme = keyN[13:end]
            sgn[i] = 1.0
        end
        
        
        if nme == "BH"
            outPix[i] = -1.0
            continue
        elseif (nme == "Inf")||(nme == "GW")
            outPix[i] = 0.0
            continue
        else
            outPix[i] = get_state_idx(nme, Nmax)
        end
            
    end
    return outPix, sgn
end

function get_state_idx(str_nlm, Nmax)
    cnt = 1
    out_idx = -1

    for nn in 1:Nmax, l in 1:(nn - 1),  m in 1:l

        if str_nlm == string(nn)*string(l)*string(m)
            out_idx = cnt
            found = true
            break
        end
        cnt += 1
    end


    if out_idx == -1
        seen_truncation_modes = Set()
        for nn in 1:Nmax, l in 1:(nn - 1)
            m_new = (2 * l)
            if m_new < Nmax
                continue # already included
            else
                truncation_key = (m_new + 1, m_new, m_new)
                if !(truncation_key in seen_truncation_modes)
                    if str_nlm == string(m_new + 1)*string(m_new)*string(m_new)
                        out_idx = cnt
                        break
                    end
                    push!(seen_truncation_modes, truncation_key)
                    cnt += 1
                end
            end
        end
    end
    return out_idx
end
