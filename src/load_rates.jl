
function load_rate_coeffs(mu, M, a, f_a; non_rel=true, input_data="Me", solve_n4=true, solve_n5=false)
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
        Drate["211_211^322^BH"] = 4.2e-7 .* alph^11 .* faFac * rP
        Drate["322_322^211^Inf"] = 1.1e-8 * alph^8 .* faFac
        Drate["211_211^GW"] = 1.0e-2 * alph^14
        Drate["322_211^GW"] = 5.0e-6 * alph^10
        Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        Drate["322_322^GW"] = 3.0e-8 * alph^18
            
        # n = 4
        if solve_n4
            Drate["211_411^322^BH"] = 2.5e-8 * alph^11 * faFac * rP
            Drate["411_411^322^BH"] = 1.7e-11 * alph^11 * faFac * rP ### Disagree
            Drate["322_411^211^Inf"] = 3.7e-9 * alph^8 * faFac
            Drate["211_211^422^BH"] = 1.5e-7 * alph^11 * faFac * rP
            Drate["411_422^211^Inf"] = 2.2e-9 * alph^8 * faFac
            Drate["411_411^422^BH"] = 2.2e-11 * alph^7 * faFac * rP
            Drate["211_422^433^BH"] = 7.83e-11 * alph^7 * faFac * rP ### Disagree
            Drate["422_322^211^Inf"] = 1.6e-8 * alph^8 * faFac
            Drate["433_433^211^Inf"] = 9.2e-11 * alph^8 * faFac
            Drate["322_433^211^Inf"] = 2.6e-9 * alph^8 * faFac
            Drate["211_322^433^BH"] = 9.1e-8 * alph^11 * faFac * rP
            Drate["322_411^433^BH"] = 3.8e-11 * alph^7 * faFac * rP
                
            if solve_n5
                Drate["322_322^544^BH"] = 1.9e-9 * alph^11 * faFac * rP
                Drate["322_411^533^BH"] = 2.0e-8 * alph^11 * faFac * rP
                Drate["322_422^544^BH"] = 1.3e-11 * alph^11 * faFac * rP #
                Drate["411_411^522^BH"] = 9.0e-11 * alph^11 * faFac * rP #
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
            
    elseif !non_rel
        
        Drate["211_211^GW"] = 1.0e-2 * alph^14
        Drate["322_211^GW"] = 5.0e-6 * alph^10
        Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        Drate["322_322^GW"] = 3.0e-8 * alph^18
        
        
        # rates computed for fixed rP(a=0.98)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.98^2))
        
        include_m1 = true
        include_m2 = true
        include_m3 = true
        # careful skip?
        if (alph < 0.2) && (a > 0.95)
            nothing;
        else
            ### check if m=1 SR
            wI_test = find_im_part(mu, M, a, 2, 1, 1;  for_s_rates=true, QNM=false, Ntot_force=5000)
            if wI_test < 0.0
                include_m1 = false
            end
        
            ### check if m=2 SR
            wI_test = find_im_part(mu, M, a, 3, 2, 2;  for_s_rates=true, QNM=false, Ntot_force=5000)
            if wI_test < 0.0
                include_m2 = false
            end
            
            ### check if m=3 SR
            wI_test = find_im_part(mu, M, a, 4, 3, 3;  for_s_rates=true, QNM=false, Ntot_force=5000)
            if wI_test < 0.0
                include_m3 = false
            end
        end
        
        dirN = "rate_sve/"
        ftag = "_GF_v2_"
        RateLoad_List = ["211_211_322_BH", "322_322_211_Inf", "211_411_322_BH", "411_411_322_BH",
                         "322_411_211_Inf", "211_211_422_BH", "411_422_211_Inf", "411_411_422_BH",
                         "211_422_433_BH", "422_322_211_Inf", "433_433_211_Inf", "322_433_211_Inf",
                         "211_322_433_BH", "322_411_433_BH", "322_322_544_BH", "322_411_533_BH",
                         "322_422_544_BH", "411_411_522_BH", "322_522_544_BH", "422_422_544_BH",
                         "422_522_544_BH", "522_522_544_BH", "422_544_322_Inf", "433_544_322_Inf",
                         "422_533_322_Inf", "433_533_322_Inf", "422_522_322_Inf", "522_544_322_Inf",
                         "533_544_322_Inf", "544_544_322_Inf", "522_533_322_Inf", "533_533_322_Inf",
                         "522_522_322_Inf"]

        # load files and get rates
        #######################
        
        for i in 1:length(RateLoad_List)
            
            n1 = RateLoad_List[i][1:1]
            n2 = RateLoad_List[i][5:5]
            n3 = RateLoad_List[i][9:9]
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
                    data = open(readdlm, dirN * RateLoad_List[i] * ftag * ".dat")
                    data = data[data[:,2] .!= 0.0, :]
                    itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                    rate_out = 10 .^itp(log10.(alph)) .* faFac
                    if occursin("BH", RateLoad_List[i])
                        rate_out *= rP_ratio
                    end
                    Drate[RateLoad_List[i]] = rate_out
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
            elseif nme == "322"
                outPix[i] = 2
            elseif nme == "411"
                outPix[i] = 3
            elseif nme == "422"
                outPix[i] = 4
            elseif nme == "433"
                outPix[i] = 5
            elseif nme == "BH"
                outPix[i] = -1
            elseif (nme == "GW")||(nme == "Inf")
                outPix[i] = 0
            else
                if solve_n5
                    if nme == "522"
                        outPix[i] = 6
                    elseif nme == "533"
                        outPix[i] = 7
                    elseif nme == "544"
                        outPix[i] = 8
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
            elseif nme == "322"
                outPix[i] = 2
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
