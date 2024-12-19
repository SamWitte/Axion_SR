using Glob

function load_rate_coeffs(mu, M, a, f_a; non_rel=true, input_data="Me", solve_n4=true, solve_n5=true, amin_211=nothing, amin_322=nothing, amin_433=nothing)
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4
    
    Drate = Dict()
    
    Drate["211_322^GW"] = 0.0 * alph^16
    
    
    if (input_data == "Doddy")
    
        Drate["211_211^322^BH"] = 0.0 # format gamma_{A x B}^{C x D}
        Drate["322_322^211^Inf"] = 0.0
        Drate["211_211^GW"] = 0.0
        Drate["322_211^GW"] = 0.0
        Drate["211_211_211^Inf"] = 0.0
        Drate["322_322^GW"] = 0.0
        
      
    elseif non_rel
        include_m1 = true
        include_m2 = true
        include_m3 = true
        if (alph < 0.2) && (a > 0.95)
            nothing;
        else
            ### check if m=1 SR
            if isnothing.(amin_211)
                wI_test = find_im_part(mu, M, a, 2, 1, 1;  for_s_rates=true, QNM=false, Ntot_force=5000)
                if wI_test < 0.0
                    include_m1 = false
                end
            else
                if a < amin_211
                    include_m1 = false
                end
            end
        
            ### check if m=2 SR
            if isnothing.(amin_322)
                wI_test = find_im_part(mu, M, a, 3, 2, 2;  for_s_rates=true, QNM=false, Ntot_force=5000)
                if wI_test < 0.0
                    include_m2 = false
                end
            else
                if a < amin_322
                    include_m2 = false
                end
            end
            
            ### check if m=3 SR
            if isnothing.(amin_433)
                wI_test = find_im_part(mu, M, a, 4, 3, 3;  for_s_rates=true, QNM=false, Ntot_force=5000)
                if wI_test < 0.0
                    include_m3 = false
                end
            else
                if a < amin_433
                    include_m3 = false
                end
            end
        end
        
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
        if solve_n4
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
           
                
            if solve_n5
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
        Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        Drate["322_322^GW"] = 3.0e-8 * alph^18
        
        
        # rates computed for fixed rP(a=0.95)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.95^2))
        
        include_m1 = true
        include_m2 = true
        include_m3 = true
        # careful skip?
        if (alph < 0.2) && (a > 0.95)
            nothing;
        else
            ### check if m=1 SR
            if isnothing.(amin_211)
                wI_test = find_im_part(mu, M, a, 2, 1, 1;  for_s_rates=true, QNM=false, Ntot_force=5000)
                if wI_test < 0.0
                    include_m1 = false
                end
            else
                if a < amin_211
                    include_m1 = false
                end
            end
        
            ### check if m=2 SR
            if isnothing.(amin_322)
                wI_test = find_im_part(mu, M, a, 3, 2, 2;  for_s_rates=true, QNM=false, Ntot_force=5000)
                if wI_test < 0.0
                    include_m2 = false
                end
            else
                if a < amin_322
                    include_m2 = false
                end
            end
            
            ### check if m=3 SR
            if isnothing.(amin_433)
                wI_test = find_im_part(mu, M, a, 4, 3, 3;  for_s_rates=true, QNM=false, Ntot_force=5000)
                if wI_test < 0.0
                    include_m3 = false
                end
            else
                if a < amin_433
                    include_m3 = false
                end
            end
        end
        
        dirN = "rate_sve/"
        ftag = "_GF_v2_"
        RateLoad_List = glob(dirN * "*" * ftag * "*.dat")
                         
                         

        # load files and get rates
        #######################
        
        for i in 1:length(RateLoad_List)
            st_idx = findfirst(dirN, RateLoad_List[i])[end] + 1
            cut_idx = findfirst("_GF_", RateLoad_List[i])[1] - 1
            n1 = RateLoad_List[i][st_idx:end][1:1]
            n2 = RateLoad_List[i][st_idx:end][5:5]
            n3 = RateLoad_List[i][st_idx:end][9:9]
            maxn = maximum([n1 n2 n3])
            
            
            # first check if rate is at correct b
            if (maxn == "3")||((maxn == "4")&&(solve_n4))||((maxn == "5")&&(solve_n5))
                # now check if alpha is too large for rate
                if (occursin("1_", RateLoad_List[i]))&&(!include_m1)
                    continue
                elseif (occursin("2_", RateLoad_List[i]))&&(!include_m2)
                    continue
                elseif (occursin("3_", RateLoad_List[i]))&&(!include_m3)
                    continue
                else
                    data = open(readdlm, RateLoad_List[i])
                    data = data[data[:,2] .!= 0.0, :]
                    # data[data[:,2] .!= 0.0, 2] .= 1e-100
                    itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                    rate_out = 10 .^itp(log10.(alph)) .* faFac
                    if occursin("BH", RateLoad_List[i])
                        rate_out *= rP_ratio
                    end
                    Drate[RateLoad_List[i][st_idx:cut_idx]] = rate_out
                end
            end
        end
    end

    return Drate
end

function key_to_indx(keyN; solve_n4=false, solve_n5=false)
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

        
            
        if solve_n4
            if nme == "211"
                outPix[i] = 1
            elseif nme == "311"
                outPix[i] = 2
            elseif nme == "322"
                outPix[i] = 3
            elseif nme == "411"
                outPix[i] = 4
            elseif nme == "422"
                outPix[i] = 5
            elseif nme == "433"
                outPix[i] = 6
            elseif nme == "BH"
                outPix[i] = -1
            elseif (nme == "GW")||(nme == "Inf")
                outPix[i] = 0
            else
                if solve_n5
                    if nme == "511"
                        outPix[i] = 7
                    elseif nme == "522"
                        outPix[i] = 8
                    elseif nme == "533"
                        outPix[i] = 9
                    elseif nme == "544"
                        outPix[i] = 10
                    else
                        print("WHAT!? \t ", nme, "\n")
                    end
                else
                    print("WHAT!? \t ", nme, "\n")
                end
            end
        else
            if nme == "211"
                outPix[i] = 1
            elseif nme == "311"
                outPix[i] = 2
            elseif nme == "322"
                outPix[i] = 3
            elseif nme == "BH"
                outPix[i] = -1
            elseif (nme == "GW")||(nme == "Inf")
                outPix[i] = 0
            else
                print("WHAT!? \t ", nme, "\n")
            end
        end
    end
    return outPix, sgn
end


