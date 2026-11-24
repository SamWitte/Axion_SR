"""
    eigenvalue_computation.jl

Documented API stubs of the eigenvalue solvers; the implementations are in solve_sr_rates.jl.
"""

module EigenvalueComputation

export solve_radial, find_im_part, find_im_zero, compute_gridded, radial_inf, spheroidals

# Note: Functions are defined in solve_sr_rates.jl which this module complements
# This is a structural organization module for cleaner code navigation

"""
    solve_radial(mu, M, a, n, l, m; kwargs) -> (wf, wf_inf, erg)

Radial wave functions and, if return_erg, the eigenvalue ω (Im ω > 0: superradiant growth) of level (n, l, m).
"""
function solve_radial(mu, M, a, n, l, m; rpts=1000, rmaxT=50, debug=false,
                      iter=500, xtol=1e-50, ftol=1e-90, sve=false,
                      fnm="test_store/WF_", return_erg=false, Ntot_safe=5000,
                      eps_r=1e-10, QNM=false, QNM_ergs=nothing,
                      pre_compute_erg=nothing, prec=300, use_heunc=true)
    # This function is defined in solve_sr_rates.jl
    # This docstring provides documentation for the API
    error("This function must be called from solve_sr_rates context")
end

"""
    find_im_part(mu, M, a, n, l, m; kwargs) -> im_part::Float64

Imaginary part of the eigenvalue of (n, l, m) (> 0: superradiant growth, < 0: damping).
"""
function find_im_part(mu, M, a, n, l, m; debug=false, Ntot_force=5000,
                      iter=500, xtol=1e-50, ftol=1e-90, return_both=false,
                      for_s_rates=true, QNM=false, QNM_E=1.0,
                      erg_Guess=nothing, max_n_qnm=10)
    error("This function must be called from solve_sr_rates context")
end

"""
    find_im_zero(mu, M, n, l, m; kwargs) -> zero_spin::Float64

BH spin at which Im ω of (n, l, m) crosses zero (onset of superradiance).
"""
function find_im_zero(mu, M, n, l, m; debug=false, Ntot=3000,
                      iter=1000, xtol=1e-6, ftol=1e-20)
    error("This function must be called from solve_sr_rates context")
end

"""
    compute_gridded(mu, M, a, n, l, m; kwargs) -> eigenvalues

Eigenvalues of (n, l, m) on a grid of npts spins from amin (also negative spins if compute_neg).
"""
function compute_gridded(mu, M, a, n, l, m; Ntot=2000, iter=50,
                        xtol=1e-7, npts=30, amin=0.0, compute_neg=false)
    error("This function must be called from solve_sr_rates context")
end

"""
    radial_inf(erg, mu, M, a, l, m; kwargs) -> wf_inf

Radial wave function at large r for eigenvalue erg (asymptotic expansion and matching).
"""
function radial_inf(erg, mu, M, a, l, m; rpts=1000, rmax_val=1e4,
                    debug=false, iter=50, xtol=1e-120, ftol=1e-120,
                    sve_for_test=false, fnm="test_store/test_radial",
                    is_infin=true)
    error("This function must be called from solve_sr_rates context")
end

"""
    spheroidals(l, m, a, erg, alph) -> coefficients

Angular eigenvalue and spheroidal harmonic coefficients; erg = ω G M, alph = μ G M, c = a sqrt(erg^2 - alph^2).
"""
function spheroidals(l, m, a, erg, alph)
    error("This function must be called from solve_sr_rates context")
end

end # module EigenvalueComputation
