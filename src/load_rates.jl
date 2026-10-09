using Glob
include("state_utils.jl")
include(joinpath(@__DIR__, "Numerics", "gw_rates.jl"))

# Leaver rate files (rate_sve/*_LvrHc_.dat), read once per session: solve_system
# with track_alpha evaluates the scattering rates on a grid of BH masses.
if !@isdefined(LVR_FILE_CACHE)
    const LVR_FILE_CACHE = Dict{String, Any}()
end

# Channels requested with self_grav=true but lacking an SG table (warned about once per session).
if !@isdefined(SG_MISSING_WARNED)
    const SG_MISSING_WARNED = Ref(false)
end

# Scattering (BH / Inf) rate coefficients. GW annihilations and transitions are
# not in this dictionary: they come from Numerics/gw_rates.jl (gw_build_cache /
# gw_rhs!) inside solve_system.
#
# self_grav = true adds the Newtonian self-gravity contribution to every channel with an SG table
# (rate_sve/<...>_LvrHc_SG.dat, written by Compute_all_rates.jl --self_grav true). The SI and SG
# amplitudes interfere, so
#     Gamma(f_a) = Gamma_SI (M_pl/f_a)^4 + si_sign * 2 kappa sqrt(Gamma_SI Gamma_SG) (M_pl/f_a)^2 + Gamma_SG,
# with kappa = cos(relative phase) = +-1 tabulated per channel and si_sign = +1 for an attractive
# quartic (axion cosine potential), -1 for a repulsive one. Only available with non_rel = false.
#
# rate_input: custom channel list (file name in rate_sve/ or absolute path, same format as
# load_rate_input_Nmax_X.txt, e.g. written by rate_sve/gen_custom_input.py) used instead of the default
# list for this Nmax. Every state in it must be one of the Nmax modes (standard or truncation).
function load_rate_coeffs(mu, M, a, f_a, Nmax, SR_rates; non_rel=true, lm_only=false, self_grav=false, si_sign=1.0, rate_input=nothing)
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4

    input_suffix = lm_only ? "_lm" : ""
    rate_file = if rate_input === nothing
        joinpath(@__DIR__, "rate_sve/load_rate_input_Nmax_$(Nmax)$(input_suffix).txt")
    else
        isabspath(rate_input) ? rate_input : joinpath(@__DIR__, "rate_sve", rate_input)
    end
    rate_list = readdlm(rate_file)
    if rate_input !== nothing
        for st in unique(string.(vec(rate_list[:, 1:3])))
            get_state_idx(st, Nmax) == -1 && error("rate_input $(rate_input): state $(st) is not among the Nmax = $(Nmax) modes; increase Nmax")
        end
    end
    # (SR_rates is no longer used: the former pruning of channels involving non-superradiant levels never
    #  removed anything, and a static t = 0 pruning would drop physical channels once spin/alpha evolve.
    #  The argument is kept so existing callers do not change.)

    Drate = Dict()
    
    if non_rel && self_grav && !SG_MISSING_WARNED[]
        println("load_rate_coeffs: self-gravity rates are only tabulated for non_rel=false (Leaver tables); ignoring self_grav.")
        SG_MISSING_WARNED[] = true
    end

    if non_rel
        include_m1 = true
        include_m2 = true
        include_m3 = true
        if include_m1
            Drate["211_211^322^BH"] = 4.2e-7 .* alph^11 .* faFac * rP
            Drate["322_322^211^Inf"] = 1.1e-8 * alph^8 .* faFac
            Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
            
            Drate["211_311^322^BH"] = 3.1e-10 .* alph^7 .* faFac * rP
            Drate["311_311^211^Inf"] = 5.1e-8 .* alph^8 * faFac
            Drate["311_322^211^Inf"] = 1.2e-8 .* alph^8 * faFac
            Drate["311_311^322^BH"] = 1.62e-10 .* alph^7 .* faFac * rP
            
        end
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
        n_sg_missing = 0
        
        # Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        
        
        # rates computed for fixed rP(a=0.9)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.95^2))
        
        dirN = joinpath(@__DIR__, "rate_sve/")
        ftag = "_LvrHc_"
        
        for i in 1:length(rate_list[:,1])
            nm_tag = string(rate_list[i, 1]) * "_" * string(rate_list[i, 2]) * "^" * string(rate_list[i, 3]) * "^" * string(rate_list[i, 4])
            fileT = dirN * string(rate_list[i, 1]) * "_" * string(rate_list[i, 2]) * "_" * string(rate_list[i, 3]) * "_" * string(rate_list[i, 4]) * ftag * ".dat"
            data = get!(LVR_FILE_CACHE, fileT) do
                isfile(fileT) || return nothing
                d = open(readdlm, fileT)
                d[d[:, 2] .!= 0.0, :]
            end
            sg = nothing
            if self_grav
                fileSG = dirN * string(rate_list[i, 1]) * "_" * string(rate_list[i, 2]) * "_" * string(rate_list[i, 3]) * "_" * string(rate_list[i, 4]) * ftag * "SG.dat"
                sg = get!(LVR_FILE_CACHE, fileSG) do
                    isfile(fileSG) || return nothing
                    d = open(readdlm, fileSG)        # alpha, Gamma_SI, Gamma_SG, Gamma_x, kappa
                    d[d[:, 3] .> 0.0, :]
                end
                sg === nothing && (n_sg_missing += 1)
            end
            if data !== nothing || sg !== nothing

                rate_out = 0.0
                g_si = 0.0
                if data !== nothing && alph .<= maximum(data[:,1])
                    itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                    g_si = 10 .^itp(log10.(alph))
                    rate_out = g_si .* faFac
                end
                if sg !== nothing && size(sg, 1) >= 2 && alph .<= maximum(sg[:, 1])
                    itp_sg = LinearInterpolation(log10.(sg[:, 1]), log10.(sg[:, 3]), extrapolation_bc=Line())
                    itp_ka = LinearInterpolation(log10.(sg[:, 1]), sg[:, 5], extrapolation_bc=Flat())
                    g_sg = 10 .^itp_sg(log10.(alph))
                    rate_out += si_sign * 2 * itp_ka(log10.(alph)) * sqrt(g_si * g_sg) * sqrt(faFac) + g_sg
                end
                if string(rate_list[i, 4]) == "BH"
                    rate_out *= rP_ratio
                end
                Drate[nm_tag] = rate_out
            end
        end
        if self_grav && n_sg_missing > 0 && !SG_MISSING_WARNED[]
            println("load_rate_coeffs: self_grav=true but $(n_sg_missing) channel(s) have no *_LvrHc_SG.dat table; ",
                    "those use the self-interaction rate only (warned once per session).")
            SG_MISSING_WARNED[] = true
        end
    end

    return Drate
end

function key_to_indx(keyN, Nmax)
    # Parse rate dictionary keys like "211_322^GW" or "211_211^322^BH"
    # Format is: state1_state2^state3_or_GW^TYPE
    # where each state can be "nlm" (old) or "n_l_m" (new)
    #
    # A degenerate N-body self-decay (e.g. "211_211_211^Inf": three identical
    # quanta annihilating to infinity, no intermediate bound state) is written
    # with all initial states underscore-joined before the single "^TYPE": all
    # of state1 and the underscore-joined "middle" states are annihilated.

    # Split by first underscore to separate state1 and remainder
    parts = split(keyN, "_", limit=2)
    if length(parts) != 2
        error("Invalid key format: $keyN")
    end

    state1 = parts[1]
    remainder = parts[2]

    # Split remainder by "^" to get state2 and rest
    caret_parts = split(remainder, "^")

    if length(caret_parts) == 2
        # Format: "state1_middle^GW/Inf/BH" or "state1_middle^state3"
        # `middle` is normally a single state (2-body process), but may itself
        # be several underscore-joined states for an N-body self-decay.
        middle_states = split(caret_parts[1], "_")
        state3_or_type = caret_parts[2]  # This is GW, Inf, BH, or a state
        n_initial = 1 + length(middle_states)
        totN = n_initial + 1

        outPix = zeros(Int, totN)
        sgn = zeros(totN)

        sgn[1] = -1.0
        outPix[1] = get_state_idx(state1, Nmax)
        for (i, st) in enumerate(middle_states)
            sgn[1 + i] = -1.0
            outPix[1 + i] = get_state_idx(String(st), Nmax)
        end

        if state3_or_type == "BH"
            outPix[totN] = -1
        elseif state3_or_type == "Inf" || state3_or_type == "GW"
            outPix[totN] = 0
        else
            # It's actually a state
            outPix[totN] = get_state_idx(state3_or_type, Nmax)
        end
        sgn[totN] = 1.0

        return outPix, sgn
    elseif length(caret_parts) == 3
        # Format: "state1_state2^state3^TYPE" (4 total components)
        totN = 4
        state2 = caret_parts[1]
        state3_or_type = caret_parts[2]  # This is state3
        type_str = caret_parts[3]  # This is BH, Inf, etc.

        outPix = zeros(Int, totN)
        sgn = zeros(totN)

        # Parse state 1
        sgn[1] = -1.0
        outPix[1] = get_state_idx(state1, Nmax)

        # Parse state 2
        sgn[2] = -1.0
        outPix[2] = get_state_idx(state2, Nmax)

        # state3_or_type is a quantum state
        sgn[3] = 1.0
        outPix[3] = get_state_idx(state3_or_type, Nmax)

        # Parse type (BH, Inf, etc.)
        sgn[4] = 1.0
        if type_str == "BH"
            outPix[4] = -1
        elseif type_str == "Inf" || type_str == "GW"
            outPix[4] = 0
        else
            error("Unknown type in key: $type_str")
        end

        return outPix, sgn
    else
        error("Invalid key format: $keyN (unexpected number of ^ separators)")
    end
end

function get_state_idx(str_nlm, Nmax)
    cnt = 1
    out_idx = -1

    for nn in 1:Nmax, l in 1:(nn - 1),  m in 1:l

        # Support both old format "211" and new format "2_1_1"
        state_str_new = format_state_string(nn, l, m)
        state_str_old = (nn < 10 && l < 10 && m < 10) ? format_state_string_legacy(nn, l, m) : ""

        if str_nlm == state_str_new || str_nlm == state_str_old
            out_idx = cnt
            found = true
            break
        end
        cnt += 1
    end


    if out_idx == -1
        seen_truncation_modes = Set()
        for nn in 1:Nmax, l in 1:(nn - 1), m in 1:l
            n_end = 2 * l + 1
            l_end = 2 * l
            m_end = 2 * m
            trunc_key = (n_end, l_end, m_end)
            if n_end > Nmax && !(trunc_key in seen_truncation_modes)
                state_str_new = format_state_string(n_end, l_end, m_end)
                state_str_old = (n_end < 10 && l_end < 10 && m_end < 10) ? format_state_string_legacy(n_end, l_end, m_end) : ""
                if str_nlm == state_str_new || str_nlm == state_str_old
                    out_idx = cnt
                    break
                end
                push!(seen_truncation_modes, trunc_key)
                cnt += 1
            end
        end
    end
    return out_idx
end
