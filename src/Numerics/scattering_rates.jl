"""
    scattering_rates.jl

Documented API stubs of the scattering-rate and frequency-shift functions; implementations in solve_sr_rates.jl.
"""

module ScatteringRates

export sr_rates, s_rate_bnd, s_rate_inf, freq_shifts

"""
    sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true) -> rate

SR rate of (n, l, m) for boson mass massB around a BH of mass MBH and spin aBH; impose_low_cut: minimum coupling.
"""
function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true)
    error("This function must be called from solve_sr_rates context")
end

"""
    s_rate_bnd(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; kwargs) -> rate

Bound-bound scattering rate (n1,l1,m1) -> (n2,l2,m2) with mode (n3,l3,m3), from wave-function overlaps.
"""
function s_rate_bnd(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3;
                    kpts=10, rpts=2000, rmaxT=100, inf_nr=true, Nang=100000,
                    Npts_Bnd=1000, debug=false, include_cont=true, Ntot_safe=5000,
                    sve_for_test=false, bnd_thresh=1e-3, use_analytic=false,
                    eps_r=1e-10, NON_REL=true)
    error("This function must be called from solve_sr_rates context")
end

"""
    s_rate_inf(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min; kwargs) -> rate

Rate of scattering from bound states into the continuum (escape to infinity), final l from lF_min.
"""
function s_rate_inf(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min;
                    rpts=4000, rmaxT=90, sve_for_test=false, inf_nr=false,
                    Npts_Bnd=1000, Nang=300000, debug=false, Ntot_safe=2000,
                    xtol=1e-2, ftol=1e-2, iter=20)
    error("This function must be called from solve_sr_rates context")
end

"""
    freq_shifts(mu, M, a, n1, l1, m1, n2, l2, m2; kwargs) -> shifts

Frequency shifts of the bound levels (n1,l1,m1) and (n2,l2,m2) from their interaction.
"""
function freq_shifts(mu, M, a, n1, l1, m1, n2, l2, m2;
                     rpts=500, rmaxT=100, Nang=100000, epsil_2=1.0,
                     Npts_Bnd=500)
    error("This function must be called from solve_sr_rates context")
end

end # module ScatteringRates
