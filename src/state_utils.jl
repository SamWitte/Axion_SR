# Utility functions for quantum state string parsing and formatting
# Supports states labeled by (n, l, m) where values can be multi-digit

"""
    parse_state_string(str_nlm::String) -> (n, l, m)

Parse "2-1-1" (hyphen-delimited, any number of digits) or legacy "211" (single digits) into (n, l, m).
"""
function parse_state_string(str_nlm::String)
    # Check if string contains delimiter
    if occursin("-", str_nlm)
        # New format: "2-1-1"
        parts = split(str_nlm, "-")
        if length(parts) != 3
            error("Invalid state string format: $str_nlm. Expected format: n-l-m")
        end
        n = parse(Int, parts[1])
        l = parse(Int, parts[2])
        m = parse(Int, parts[3])
    else
        # Legacy format: "211" (only works for single digits)
        if length(str_nlm) != 3
            error("Invalid state string format: $str_nlm. For multi-digit numbers, use delimiter format: n-l-m")
        end
        # Parse each character as a single digit
        n = parse(Int, str_nlm[1:1])
        l = parse(Int, str_nlm[2:2])
        m = parse(Int, str_nlm[3:3])
    end
    return (n, l, m)
end


"""
    format_state_string(n::Int, l::Int, m::Int) -> String

"n-l-m", e.g. "10-5-3" (hyphen, since underscores separate states in rate keys).
"""
function format_state_string(n::Int, l::Int, m::Int)
    return "$(n)-$(l)-$(m)"
end


"""
    format_state_string_legacy(n::Int, l::Int, m::Int) -> String

Undelimited "211"; errors if a quantum number is >= 10.
"""
function format_state_string_legacy(n::Int, l::Int, m::Int)
    if n >= 10 || l >= 10 || m >= 10
        error("Legacy format only supports single-digit quantum numbers. Use format_state_string() instead.")
    end
    return "$(n)$(l)$(m)"
end


"""
    parse_state_from_args(S::String) -> (n, l, m)

parse_state_string of a command-line argument.
"""
function parse_state_from_args(S::String)
    return parse_state_string(S)
end


"""
    get_m_from_state_string(S::String) -> Int

The m of a state string ("211" or "10-5-3").
"""
function get_m_from_state_string(S::String)
    n, l, m = parse_state_string(S)
    return m
end


"""
    format_state_for_filename(n::Int, l::Int, m::Int) -> String

Legacy "211" if all quantum numbers are < 10 (older file names), else "n-l-m".
"""
function format_state_for_filename(n::Int, l::Int, m::Int)
    if n < 10 && l < 10 && m < 10
        # Use legacy format for backward compatibility
        return format_state_string_legacy(n, l, m)
    else
        # Use new format for multi-digit numbers
        return format_state_string(n, l, m)
    end
end


"""
    rate_input_path(Nmax; lm_only=false, rate_input=nothing) -> String
Channel list: `rate_input` (file in rate_sve/ or absolute path), else ENV NL_DEFAULT_RATE_INPUT (output names
unchanged, e.g. to resume on a frozen list), else rate_sve/load_rate_input_Nmax_<Nmax>[_lm].txt.
"""
function rate_input_path(Nmax; lm_only=false, rate_input=nothing)
    if rate_input === nothing && !isempty(get(ENV, "NL_DEFAULT_RATE_INPUT", ""))
        rate_input = ENV["NL_DEFAULT_RATE_INPUT"]
    end
    rate_input === nothing && return joinpath(@__DIR__, "rate_sve", "load_rate_input_Nmax_$(Nmax)$(lm_only ? "_lm" : "").txt")
    return isabspath(rate_input) ? rate_input : joinpath(@__DIR__, "rate_sve", rate_input)
end


"""
    scattering_levels(Nmax; lm_only=false, rate_input=nothing) -> Set{NTuple{3,Int}}

Levels (n, l, m) in that channel list (first three columns); GW channels are kept only between these.
"""
function scattering_levels(Nmax; lm_only=false, rate_input=nothing)
    lv = Set{NTuple{3,Int}}()
    for line in eachline(rate_input_path(Nmax; lm_only=lm_only, rate_input=rate_input))
        f = split(line)
        (length(f) >= 3 && !startswith(f[1], "#")) || continue
        for s in f[1:3]
            push!(lv, parse_state_string(String(s)))
        end
    end
    return lv
end
