# Ehrenfest check for the transition 211 <-> 322 with the Kerr modes:  GM <y/y^3>  vs  w31^2 <y>  (radial parts)
include("exchange.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "solve_sr_rates.jl"))
using Printf
for (a, alph) in ((0.0, 0.035), (0.0, 0.1), (0.95, 0.1), (0.0, 0.2))
    mu = alph / GNew
    g1 = find_im_part(mu, 1.0, max(a, 0.5), 3, 2, 2; Ntot_force=5000, for_s_rates=true, return_both=true)
    g3 = find_im_part(mu, 1.0, max(a, 0.5), 2, 1, 1; Ntot_force=5000, for_s_rates=true, return_both=true)
    md1, _ = kerr_mode(3, 2, 2, alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
    md3, _ = kerr_mode(2, 1, 1, alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
    r = collect(10 .^ range(log10(rplus(a) * 1.0001), log10(200 / alph^2), 40000))
    R1 = [R_BL(md1, x) for x in r]; R3 = [R_BL(md3, x) for x in r]
    rho = R1 .* conj.(R3)
    trap(f) = sum(0.5 .* (f[1:end-1] .+ f[2:end]) .* diff(r))
    D1 = trap(rho)                 # ∫ rho r^2 * r / r^3 dr
    D2 = trap(rho .* r .^ 3)       # ∫ rho r^2 * r dr
    w31 = real(md1.w) - real(md3.w)
    # hydrogenic for comparison
    Rh1 = radial_bound_NR(3, 2, 2, mu, 1.0, r); Rh3 = radial_bound_NR(2, 1, 1, mu, 1.0, r)
    rh = Float64.(Rh1 .* Rh3)
    w31h = alph * (1 / 8 - 1 / 18) * alph^2
    @printf("a=%.2f alpha=%.3f:  Kerr modes: GM<y/y^3> / (w31^2 <y>) = %.5f    hydrogenic: %.5f   (w31 Kerr/hyd = %.5f)\n",
            a, alph, real(D1 / (w31^2 * D2)), trap(rh) / (w31h^2 * trap(rh .* r .^ 3)), w31 / w31h)
end
