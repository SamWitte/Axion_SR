"""
Rate coefficients of the axion-BH scattering channels: RateCoefficient (one channel) and RateDatabase
(all channels for an Nmax).
"""

module RateCoefficients

using LinearAlgebra

"""
    RateCoefficient

Channel rate (key e.g. "211_211^322^BH", channel "BH" or "Inf") = coefficient * alpha^power_index [* fa_factor] [* rp].
"""
struct RateCoefficient
    key::String
    level1::String
    level2::String
    channel::String
    coefficient::Float64
    power_index::Int
    requires_rp::Bool
    requires_fa::Bool
end

"""
    evaluate_rate(rc::RateCoefficient, alpha::Float64, fa_factor::Float64=1.0, rp::Float64=1.0)

coefficient * alpha^power_index, times fa_factor = (M_pl/f_a)^4 and rp where the coefficient requires them.
"""
function evaluate_rate(rc::RateCoefficient, alpha::Float64, fa_factor::Float64=1.0, rp::Float64=1.0)
    rate = rc.coefficient * (alpha ^ rc.power_index)

    if rc.requires_fa
        rate *= fa_factor
    end

    if rc.requires_rp
        rate *= rp
    end

    return rate
end

"""
    RateDatabase

Rate coefficients by key for one Nmax (and non_relativistic flag).
"""
mutable struct RateDatabase
    coefficients::Dict{String, RateCoefficient}
    nmax::Int
    non_relativistic::Bool

    RateDatabase(nmax::Int; non_relativistic::Bool=true) = new(Dict(), nmax, non_relativistic)
end

"""
    add_rate!(db::RateDatabase, rc::RateCoefficient)

Add a rate coefficient to the database.
"""
function add_rate!(db::RateDatabase, rc::RateCoefficient)
    db.coefficients[rc.key] = rc
end

"""
    get_rate(db::RateDatabase, key::String)

The RateCoefficient stored under `key`, or nothing.
"""
function get_rate(db::RateDatabase, key::String)
    return get(db.coefficients, key, nothing)
end

"""
    list_rates(db::RateDatabase; channel::Union{String, Nothing}=nothing)

Keys of all rate coefficients, or of those of `channel` ("BH" or "Inf").
"""
function list_rates(db::RateDatabase; channel::Union{String, Nothing}=nothing)
    if channel === nothing
        return collect(keys(db.coefficients))
    else
        return [key for key in keys(db.coefficients)
                if db.coefficients[key].channel == channel]
    end
end

"""
    count_rates(db::RateDatabase; channel::Union{String, Nothing}=nothing)

Number of rate coefficients (optionally of one channel).
"""
function count_rates(db::RateDatabase; channel::Union{String, Nothing}=nothing)
    return length(list_rates(db; channel=channel))
end

"""
    build_default_rates(Nmax::Int; non_relativistic::Bool=true)

RateDatabase with the tabulated non-relativistic coefficients (Nmax 3-5 channels; empty if !non_relativistic).
"""
function build_default_rates(Nmax::Int; non_relativistic::Bool=true)
    db = RateDatabase(Nmax; non_relativistic=non_relativistic)

    if non_relativistic
        # Nmax >= 3 rates
        add_rate!(db, RateCoefficient("211_211^322^BH", "211", "211", "BH", 4.2e-7, 11, true, true))
        add_rate!(db, RateCoefficient("322_322^211^Inf", "322", "322", "Inf", 1.1e-8, 8, false, true))
        add_rate!(db, RateCoefficient("211_211_211^Inf", "211", "211", "Inf", 1.5e-8, 21, false, true))

        add_rate!(db, RateCoefficient("211_311^322^BH", "211", "311", "BH", 3.1e-10, 7, true, true))
        add_rate!(db, RateCoefficient("311_311^211^Inf", "311", "311", "Inf", 5.1e-8, 8, false, true))
        add_rate!(db, RateCoefficient("311_322^211^Inf", "311", "322", "Inf", 1.2e-8, 8, false, true))
        add_rate!(db, RateCoefficient("311_311^322^BH", "311", "311", "BH", 1.62e-10, 7, true, true))


        if Nmax >= 4
            # Nmax >= 4 rates
            add_rate!(db, RateCoefficient("211_411^322^BH", "211", "411", "BH", 2.5e-8, 11, true, true))
            add_rate!(db, RateCoefficient("322_411^211^Inf", "322", "411", "Inf", 3.7e-9, 8, false, true))
            add_rate!(db, RateCoefficient("211_211^422^BH", "211", "211", "BH", 1.5e-7, 11, true, true))
            add_rate!(db, RateCoefficient("411_422^211^Inf", "411", "422", "Inf", 2.2e-9, 8, false, true))
            add_rate!(db, RateCoefficient("411_411^322^BH", "411", "411", "BH", 1.7e-11, 11, true, true))
            add_rate!(db, RateCoefficient("411_411^422^BH", "411", "411", "BH", 2.2e-11, 7, true, true))
            add_rate!(db, RateCoefficient("211_422^433^BH", "211", "422", "BH", 7.83e-11, 7, true, true))

            # Additional Nmax >= 4 rates
            add_rate!(db, RateCoefficient("411_411^211^Inf", "411", "411", "Inf", 1.7e-9, 8, false, true))
            add_rate!(db, RateCoefficient("411_433^211^Inf", "411", "433", "Inf", 1.1e-10, 8, false, true))
            add_rate!(db, RateCoefficient("422_422^211^Inf", "422", "422", "Inf", 1.6e-9, 8, false, true))
            add_rate!(db, RateCoefficient("422_433^211^Inf", "422", "433", "Inf", 6.1e-10, 8, false, true))
            add_rate!(db, RateCoefficient("211_411^422^BH", "211", "411", "BH", 3.2e-11, 7, true, true))
            add_rate!(db, RateCoefficient("411_422^433^BH", "411", "422", "BH", 2.3e-11, 7, true, true))

            add_rate!(db, RateCoefficient("211_311^422^BH", "211", "311", "BH", 2.7e-7, 11, true, true))
            add_rate!(db, RateCoefficient("311_311^422^BH", "311", "311", "BH", 1.7e-11, 11, true, true))
            add_rate!(db, RateCoefficient("311_322^433^BH", "311", "322", "BH", 7.0e-8, 11, true, true))
            add_rate!(db, RateCoefficient("311_411^211^Inf", "311", "411", "Inf", 1.9e-8, 8, false, true))
            add_rate!(db, RateCoefficient("311_422^211^Inf", "311", "422", "Inf", 7.0e-9, 8, false, true))
            add_rate!(db, RateCoefficient("311_433^211^Inf", "311", "433", "Inf", 2.2e-10, 8, false, true))
            add_rate!(db, RateCoefficient("311_411^322^BH", "311", "411", "BH", 1.9e-10, 7, true, true))
            add_rate!(db, RateCoefficient("311_411^422^BH", "311", "411", "BH", 3.8e-13, 7, true, true))
            add_rate!(db, RateCoefficient("311_422^433^BH", "311", "422", "BH", 7.7e-12, 7, true, true))

            add_rate!(db, RateCoefficient("422_322^211^Inf", "422", "322", "Inf", 1.6e-8, 8, false, true))
            add_rate!(db, RateCoefficient("433_433^211^Inf", "433", "433", "Inf", 9.2e-11, 8, false, true))
            add_rate!(db, RateCoefficient("322_433^211^Inf", "322", "433", "Inf", 2.6e-9, 8, false, true))
            add_rate!(db, RateCoefficient("211_322^433^BH", "211", "322", "BH", 9.1e-8, 11, true, true))
            add_rate!(db, RateCoefficient("322_411^433^BH", "322", "411", "BH", 3.8e-11, 7, true, true))
        end

        if Nmax >= 5
            # Nmax >= 5 rates
            add_rate!(db, RateCoefficient("211_211^522^BH", "211", "211", "BH", 7.5e-8, 11, true, true))
            add_rate!(db, RateCoefficient("322_411^533^BH", "322", "411", "BH", 2.0e-8, 11, true, true))
            add_rate!(db, RateCoefficient("411_411^522^BH", "411", "411", "BH", 9.0e-11, 11, true, true))

            # Additional Nmax >= 5 rates
            add_rate!(db, RateCoefficient("211_311^522^BH", "211", "311", "BH", 1.0e-7, 11, true, true))
            add_rate!(db, RateCoefficient("211_411^522^BH", "211", "411", "BH", 9.9e-8, 11, true, true))
            add_rate!(db, RateCoefficient("211_322^533^BH", "211", "322", "BH", 3.1e-8, 11, true, true))
            add_rate!(db, RateCoefficient("211_422^533^BH", "211", "422", "BH", 1.1e-7, 11, true, true))
            add_rate!(db, RateCoefficient("211_433^544^BH", "211", "433", "BH", 1.1e-9, 11, true, true))
            add_rate!(db, RateCoefficient("211_511^322^BH", "211", "511", "BH", 2.9e-8, 11, true, true))
            add_rate!(db, RateCoefficient("211_511^422^BH", "211", "511", "BH", 2.1e-10, 11, true, true))
            add_rate!(db, RateCoefficient("211_522^433^BH", "211", "522", "BH", 6.5e-8, 11, true, true))
            add_rate!(db, RateCoefficient("211_511^522^BH", "211", "511", "BH", 6.6e-12, 7, true, true))
            add_rate!(db, RateCoefficient("211_522^533^BH", "211", "522", "BH", 2.6e-11, 7, true, true))
            add_rate!(db, RateCoefficient("211_533^544^BH", "211", "533", "BH", 4.6e-13, 7, true, true))
        end
    end

    return db
end

"""
    export_to_dict(db::RateDatabase, alpha::Float64, fa_factor::Float64=1.0, rp::Float64=1.0)

Dict of key => evaluate_rate(rc, alpha, fa_factor, rp) for every coefficient.
"""
function export_to_dict(db::RateDatabase, alpha::Float64, fa_factor::Float64=1.0, rp::Float64=1.0)
    result = Dict{String, Float64}()

    for (key, rc) in db.coefficients
        result[key] = evaluate_rate(rc, alpha, fa_factor, rp)
    end

    return result
end

# Export public API
export RateCoefficient, RateDatabase, evaluate_rate, add_rate!, get_rate, list_rates
export count_rates, build_default_rates, export_to_dict

end  # module
