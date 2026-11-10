# Float64 vs BigFloat Teukolsky annihilation rates at the lowest tabulated alpha (0.05): the
# annihilation integrals cancel to ~alpha^(4l+10) of their bulk, so this decides whether Float64
# tables are trustworthy there for the higher-l pairs.
#   julia scripts/diag_gw_kerr_bigfloat.jl
using Suppressor
include(joinpath(@__DIR__, "..", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "state_utils.jl"))
@suppress include(joinpath(@__DIR__, "..", "Numerics", "gw_rates.jl"))
include(joinpath(@__DIR__, "..", "Numerics", "gw_kerr_build.jl"))

for (sa, α, a) in (((2, 1, 1), 0.05, 0.5), ((3, 2, 2), 0.05, 0.5), ((4, 3, 3), 0.05, 0.5), ((5, 4, 4), 0.05, 0.9), ((6, 5, 5), 0.05, 0.9),
                   ((3, 2, 2), 0.1, 0.9), ((5, 4, 4), 0.1, 0.9))
    ma, _ = kerr_mode(sa..., α, a)
    t0 = time(); Γ64, P64, m64 = gw_kerr_rate(ma, ma, :ann; lmax=sa[3] * 2 + 10)
    t1 = time(); ΓB, PB, mB = gw_kerr_rate(ma, ma, :ann; lmax=sa[3] * 2 + 10, T=BigFloat)
    t2 = time()
    C, p = gw_nr_annihilation_coeff(sa, sa)
    println("$(sa) α=$α a=$a: Float64 Γ=$(Γ64) ($(round(t1 - t0, digits=1)) s), BigFloat Γ=$(ΓB) ($(round(t2 - t1, digits=1)) s), " *
            "rel diff $(abs(Γ64 / ΓB - 1)), Γ/Γ_NR = $(ΓB / (C * α^p)), nℓ=$(length(m64))")
    flush(stdout)
end
