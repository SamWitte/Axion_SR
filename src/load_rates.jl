using Glob
include("state_utils.jl")
include(joinpath(@__DIR__, "Numerics", "gw_rates.jl"))
include(joinpath(@__DIR__, "spin_table.jl"))
@isdefined(table_read) || include(joinpath(@__DIR__, "Core", "table_io.jl"))

# rate tables, read once per session
if !@isdefined(LVR_FILE_CACHE)
    const LVR_FILE_CACHE = Dict{String, Any}()
end

# channels with self_grav = true but no SG table (warned once per session)
if !@isdefined(SG_MISSING_WARNED)
    const SG_MISSING_WARNED = Ref(false)
end

# rate_sve/regime_split.txt: "S1 S2 S3 alpha_c [LOW]"; the fourth quantum is LOW (BH/Inf, default BH) below alpha_c.
if !@isdefined(REGIME_SPLIT_CACHE)
    const REGIME_SPLIT_CACHE = Dict{String, Dict{String, Tuple{Float64, String}}}()
end
function regime_split_table(dirN)
    get!(REGIME_SPLIT_CACHE, dirN) do
        tab = Dict{String, Tuple{Float64, String}}()
        split_file = joinpath(dirN, "regime_split.txt")
        if table_exists(split_file)
            for ln in table_read(readlines, split_file)
                p = split(ln)
                if length(p) in (4, 5) && !startswith(p[1], "#")
                    tab[p[1] * "_" * p[2] * "_" * p[3]] =
                        (parse(Float64, p[4]), length(p) == 5 ? String(p[5]) : "BH")
                end
            end
        end
        tab
    end
end

# Scattering (BH / Inf) rate coefficients. GW annihilations and transitions are
# not in this dictionary: they come from Numerics/gw_rates.jl (gw_build_cache /
# gw_rhs!) inside solve_system.
#
# self_grav = true: Gamma(f_a) = Gamma_SI c^2 + 2 kappa c sqrt(Gamma_SI Gamma_SG) + Gamma_SG with
# c = si_sign (M_pl/f_a)^2 + 12 pi G M_pl^2 (s-channel graviton term; NL_S_CHANNEL=0 drops it), from the SG tables
# (rate_sve/*_LvrHc_SG.dat, non_rel = false). NL_SI_FROM_SG=1 takes Gamma_SI from the SG table.
# rate_input: channel list in rate_sve/ or an absolute path (format of load_rate_input_Nmax_X.txt).
# regime_split: channels in rate_sve/regime_split.txt get BH and Inf rows, each used on its side of alpha_c
# (regime_mask = false: no cut). parts: filled with (SI, cross, SG) terms per channel (spin tables).
function load_rate_coeffs(mu, M, a, f_a, Nmax, SR_rates; non_rel=true, lm_only=false, self_grav=false, si_sign=1.0, rate_input=nothing,
                          regime_split=true, regime_mask=true, parts=nothing)
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4
    si_from_sg = get(ENV, "NL_SI_FROM_SG", "0") == "1"
    sg_overlay = get(ENV, "NL_SG_OVERLAY", "")
  # s-channel graviton term
    s_chan = self_grav && !non_rel && get(ENV, "NL_S_CHANNEL", "1") != "0"
    c_contact = si_sign * sqrt(faFac) + 12pi * (GNew / M_to_eV) * (M_pl * 1e9)^2
    si_fac = s_chan ? c_contact^2 : faFac          # Gamma_SI is multiplied by this

    rate_file = rate_input_path(Nmax; lm_only=lm_only, rate_input=rate_input)
    rate_list = table_read(readdlm, rate_file)
    if rate_input !== nothing
        for st in unique(string.(vec(rate_list[:, 1:3])))
            get_state_idx(st, Nmax) == -1 && error("rate_input $(rate_input): state $(st) is not among the Nmax = $(Nmax) modes; increase Nmax")
        end
    end
    # (SR_rates unused, kept for the callers)

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
        
        rows = [Tuple(string.(rate_list[i, 1:4])) for i in 1:size(rate_list, 1)]
        split_tab = regime_split ? regime_split_table(dirN) : Dict{String, Tuple{Float64, String}}()
        if !isempty(split_tab)
            have = Set(rows)
            for (s1, s2, s3, s4) in copy(rows)
                other = (s1, s2, s3, s4 == "BH" ? "Inf" : "BH")
                if s4 in ("BH", "Inf") && haskey(split_tab, s1 * "_" * s2 * "_" * s3) && !(other in have)
                    push!(rows, other)
                    push!(have, other)
                end
            end
        end

        for (s1, s2, s3, s4) in rows
            nm_tag = s1 * "_" * s2 * "^" * s3 * "^" * s4
            fileT = dirN * s1 * "_" * s2 * "_" * s3 * "_" * s4 * ftag * ".dat"
            # without the split, a split channel's unsplit original (kept by split_regime.py)
            orig = dirN * "regime_split_orig/" * basename(fileT)
            !regime_split && table_exists(orig) && (fileT = orig)
            data = get!(LVR_FILE_CACHE, fileT) do
                table_exists(fileT) || return nothing
                d = table_read(p -> open(readdlm, p), fileT)
                d[d[:, 2] .!= 0.0, :]
            end
            split_ch = get(split_tab, s1 * "_" * s2 * "_" * s3, nothing)
            sg = nothing
            if self_grav || si_from_sg  # NL_SI_FROM_SG=1 without self_grav: SG table read for Gamma_SI only
                fileSG = dirN * s1 * "_" * s2 * "_" * s3 * "_" * s4 * ftag * "SG.dat"
                # NL_SG_OVERLAY=<dir>: tables there replace those in rate_sve
                if !isempty(sg_overlay)
                    alt = joinpath(sg_overlay, basename(fileSG))
                    table_exists(alt) && (fileSG = alt)
                end
                sg = get!(LVR_FILE_CACHE, fileSG) do
                    table_exists(fileSG) || return nothing
                    d = table_read(p -> open(readdlm, p), fileSG)        # alpha, Gamma_SI, Gamma_SG, Gamma_x, kappa
                    d[d[:, 3] .> 0.0, :]
                end
                self_grav && sg === nothing && (n_sg_missing += 1)
            end
            if data !== nothing || sg !== nothing

                rate_out = 0.0
                g_si = 0.0
                p_x = 0.0; p_sg = 0.0       # cross and SG terms, for `parts`
                if data !== nothing && split_ch !== nothing
                    # LOW half: alpha < alpha_c. Other half: alpha_c <= alpha <= top of its grid.
                    a_c, low = split_ch
                    in_regime = (s4 == low) ? (alph < a_c) : (alph >= a_c)
                    if (in_regime || !regime_mask) && (s4 == low || alph <= maximum(data[:, 1]))
                        if size(data, 1) == 1
                            g_si = data[1, 2]
                        else
                            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                            g_si = 10 .^itp(log10.(alph))
                        end
                        rate_out = g_si .* si_fac
                    end
                elseif data !== nothing && alph .<= maximum(data[:,1])
                    itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                    g_si = 10 .^itp(log10.(alph))
                    rate_out = g_si .* si_fac
                end
                # a split channel has the same SG table under both labels: same mask as SI
                sg_in_regime = split_ch === nothing || !regime_mask ||
                               ((s4 == split_ch[2]) ? (alph < split_ch[1]) : (alph >= split_ch[1]))
                if sg !== nothing && size(sg, 1) >= 2 && alph .<= maximum(sg[:, 1]) && sg_in_regime
                    itp_sg = LinearInterpolation(log10.(sg[:, 1]), log10.(sg[:, 3]), extrapolation_bc=Line())
                    itp_ka = LinearInterpolation(log10.(sg[:, 1]), sg[:, 5], extrapolation_bc=Flat())
                    g_sg = 10 .^itp_sg(log10.(alph))
                    # NL_SI_FROM_SG=1: Gamma_SI from the SG table inside its alpha range.
                    if si_from_sg && alph >= minimum(sg[:, 1]) * (1 - 1e-9)
                        pos = sg[:, 2] .> 0.0
                        if count(pos) >= 2
                            itp_si2 = LinearInterpolation(log10.(sg[pos, 1]), log10.(sg[pos, 2]), extrapolation_bc=Line())
                            g_si = 10 .^itp_si2(log10.(alph))
                            rate_out = g_si .* si_fac
                        end
                    end
                    if s_chan
                        rate_out += 2 * itp_ka(log10.(alph)) * sqrt(g_si * g_sg) * c_contact + g_sg
                    else
                        self_grav && (rate_out += si_sign * 2 * itp_ka(log10.(alph)) * sqrt(g_si * g_sg) * sqrt(faFac) + g_sg)
                    end
                    if parts !== nothing && (s_chan || self_grav)
                        p_x = 2 * itp_ka(log10.(alph)) * sqrt(g_si * g_sg) * (s_chan ? c_contact : si_sign * sqrt(faFac))
                        p_sg = g_sg
                    end
                end
                if s4 == "BH"
                    rate_out *= rP_ratio
                end
                Drate[nm_tag] = rate_out
                if parts !== nothing
                    pf = s4 == "BH" ? rP_ratio : 1.0
                    parts[nm_tag] = (Float64(g_si * si_fac * pf), Float64(p_x * pf), Float64(p_sg * pf))
                end
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

# (n, l, m) of every bound state named in a rate key ("a_b^c^BH", "a_b^c^Inf", "a_b_c^Inf", ...).
rate_key_levels(k) = [parse_state_string(String(s)) for s in split(k, ['_', '^']) if !(s in ("BH", "Inf", "GW"))]

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

# level name ("n-l-m" or legacy "nlm") -> index in the Nmax basis (regular modes, then truncation modes); -1 if absent
if !@isdefined(STATE_IDX_CACHE)
    const STATE_IDX_CACHE = Dict{Int, Dict{String, Int}}()
end
function state_idx_table(Nmax)
    get!(STATE_IDX_CACHE, Nmax) do
        tab = Dict{String, Int}()
        add(n, l, m, i) = (get!(tab, format_state_string(n, l, m), i);
                           (n < 10 && l < 10 && m < 10) && get!(tab, format_state_string_legacy(n, l, m), i))
        cnt = 1
        for nn in 1:Nmax, l in 1:(nn - 1), m in 1:l
            add(nn, l, m, cnt); cnt += 1
        end
        seen = Set{NTuple{3, Int}}()
        for nn in 1:Nmax, l in 1:(nn - 1), m in 1:l
            key = (2 * l + 1, 2 * l, 2 * m)
            if key[1] > Nmax && !(key in seen)
                add(key..., cnt); push!(seen, key); cnt += 1
            end
        end
        tab
    end
end
get_state_idx(str_nlm, Nmax) = get(state_idx_table(Nmax), str_nlm, -1)
