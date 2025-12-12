# Utility functions for quantum state string parsing and formatting
# Supports states labeled by (n, l, m) where values can be multi-digit

"""
    parse_state_string(str_nlm::String) -> (n, l, m)

Parse a state string into quantum numbers (n, l, m).
Supports both legacy format "211" (single digits) and new format "2-1-1" (delimited).
Note: Uses hyphen "-" as delimiter to avoid ambiguity with underscore in rate keys.

# Examples
```julia
parse_state_string("211")      # Returns (2, 1, 1)
parse_state_string("2-1-1")    # Returns (2, 1, 1)
parse_state_string("10-5-3")   # Returns (10, 5, 3)
```
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

Format quantum numbers (n, l, m) into a state string using hyphen delimiter.
Note: Uses hyphen "-" as delimiter to avoid ambiguity with underscore in rate keys.

# Examples
```julia
format_state_string(2, 1, 1)    # Returns "2-1-1"
format_state_string(10, 5, 3)   # Returns "10-5-3"
```
"""
function format_state_string(n::Int, l::Int, m::Int)
    return "$(n)-$(l)-$(m)"
end


"""
    format_state_string_legacy(n::Int, l::Int, m::Int) -> String

Format quantum numbers (n, l, m) into legacy state string without delimiter.
Only works for single-digit numbers (n, l, m < 10).

# Examples
```julia
format_state_string_legacy(2, 1, 1)  # Returns "211"
```
"""
function format_state_string_legacy(n::Int, l::Int, m::Int)
    if n >= 10 || l >= 10 || m >= 10
        error("Legacy format only supports single-digit quantum numbers. Use format_state_string() instead.")
    end
    return "$(n)$(l)$(m)"
end


"""
    parse_state_from_args(S::String) -> (n, l, m)

Parse a state string from command-line arguments.
Returns tuple of integers (n, l, m).

# Examples
```julia
parse_state_from_args("211")      # Returns (2, 1, 1)
parse_state_from_args("2-1-1")    # Returns (2, 1, 1)
parse_state_from_args("10-5-3")   # Returns (10, 5, 3)
```
"""
function parse_state_from_args(S::String)
    return parse_state_string(S)
end


"""
    get_m_from_state_string(S::String) -> Int

Extract the m quantum number from a state string.

# Examples
```julia
get_m_from_state_string("211")      # Returns 1
get_m_from_state_string("2-1-1")    # Returns 1
get_m_from_state_string("10-5-3")   # Returns 3
```
"""
function get_m_from_state_string(S::String)
    n, l, m = parse_state_string(S)
    return m
end


"""
    format_state_for_filename(n::Int, l::Int, m::Int) -> String

Format quantum numbers for use in filenames.
Uses legacy format (no delimiter) for single-digit numbers for backward compatibility.
Uses new format (with delimiter) for multi-digit numbers.

# Examples
```julia
format_state_for_filename(2, 1, 1)     # Returns "211" (legacy format)
format_state_for_filename(10, 5, 3)    # Returns "10-5-3" (new format)
format_state_for_filename(2, 10, 1)    # Returns "2-10-1" (new format)
```
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
