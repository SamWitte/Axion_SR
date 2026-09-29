# Build the non-relativistic GW annihilation/transition table for a given Nmax.
#   julia scripts/build_gw_tables.jl 18
# Writes rate_sve/gw_nr_rates_Nmax_<Nmax>.txt (read by load_rate_coeffs with gw_model=:nonrel).
include(joinpath(@__DIR__, "..", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "state_utils.jl"))
include(joinpath(@__DIR__, "..", "Numerics", "gw_rates.jl"))

for arg in (isempty(ARGS) ? ["18"] : ARGS)
    gw_build_nr_table(parse(Int, arg))
end
