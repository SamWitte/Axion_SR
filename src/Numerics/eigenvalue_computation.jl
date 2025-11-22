"""
    eigenvalue_computation.jl

Eigenvalue computation module for black hole superradiance analysis.

This module provides functions for solving eigenvalue problems related to
scalar field perturbations around rotating black holes, including:
- Direct eigenvalue computation via radial equation solving
- Imaginary part extraction for stability analysis
- Gridded eigenvalue computations for parameter space exploration

**Main Functions**:
- `solve_radial()` - Primary eigenvalue solver using Heun equation recursion
- `find_im_part()` - Extract imaginary eigenvalue component
- `find_im_zero()` - Find zero of imaginary part
- `compute_gridded()` - Compute eigenvalues on spin grid
- `radial_inf()` - Solve radial equation at infinity
- `spheroidals()` - Compute spheroidal harmonic coefficients

**Note**: This module provides high-level eigenvalue computation functions.
It complements the rate_computation module which handles pre-computed
data interpolation and the broader numerical infrastructure in solve_sr_rates.jl.

**Usage Example**:
```julia
using .EigenvalueComputation

# Compute eigenvalue for given parameters
_, _, erg = solve_radial(mu, M_BH, a, n, l, m;
                         rpts=1000,
                         rmaxT=50,
                         use_heunc=true,
                         return_erg=true)

# Extract just imaginary part
im_part = find_im_part(mu, M, a, n, l, m)
```
"""

module EigenvalueComputation

export solve_radial, find_im_part, find_im_zero, compute_gridded, radial_inf, spheroidals

# Note: Functions are defined in solve_sr_rates.jl which this module complements
# This is a structural organization module for cleaner code navigation

"""
    solve_radial(mu, M, a, n, l, m; kwargs) -> (wf, wf_inf, erg)

Solve radial eigenvalue equation for scalar perturbations around Kerr black hole.

Computes the radial wave function and eigenvalue (energy) for a scalar field
with quantum numbers (n,l,m) in the background of a Kerr black hole with mass M
and spin parameter a, at boson mass mu.

**Arguments**:
- `mu::Float64` - Boson mass in natural units
- `M::Float64` - Black hole mass
- `a::Float64` - Black hole spin parameter (dimensionless)
- `n::Int` - Principal quantum number
- `l::Int` - Orbital angular momentum quantum number
- `m::Int` - Azimuthal quantum number

**Keyword Arguments**:
- `rpts::Int` - Number of radial grid points (default: 1000)
- `rmaxT::Float64` - Maximum radius for computation (default: 50)
- `use_heunc::Bool` - Use Heun equation recursion relation (default: true)
- `return_erg::Bool` - Whether to return eigenvalue (default: false)
- `Ntot_safe::Int` - Safety limit for iterations (default: 5000)
- `debug::Bool` - Verbose output (default: false)
- `iter::Int` - Maximum iterations (default: 500)
- `xtol::Float64` - Position tolerance (default: 1e-50)
- `ftol::Float64` - Function tolerance (default: 1e-90)

**Returns**:
- `wf::Vector` - Radial wave function at boundary
- `wf_inf::Vector` - Radial wave function at infinity
- `erg::Complex` - Eigenvalue (if return_erg=true)

**Physical Interpretation**:
The eigenvalue ω = ωR + i·ωI has:
- Real part ωR: Oscillation frequency
- Imaginary part ωI: Decay rate (negative) or growth rate (positive)
- Im(ω) > 0 indicates superradiant instability

**Implementation Notes**:
- Uses Heun equation recursion relation by default for efficiency
- Solves boundary value problem with asymptotic matching
- Provides both bound state and quasi-normal mode solutions
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

Find imaginary part of eigenvalue for stability analysis.

Solves for the eigenvalue and returns just the imaginary component,
which indicates the instability/decay timescale.

**Arguments**:
- `mu::Float64` - Boson mass
- `M::Float64` - Black hole mass
- `a::Float64` - Black hole spin parameter
- `n, l, m::Int` - Quantum numbers

**Returns**:
- `im_part::Float64` - Imaginary part of eigenvalue

**Physical Meaning**:
- Im(ω) > 0: Superradiant instability (exponential growth)
- Im(ω) ≈ 0: Marginally stable
- Im(ω) < 0: Stable/damped mode
"""
function find_im_part(mu, M, a, n, l, m; debug=false, Ntot_force=5000,
                      iter=500, xtol=1e-50, ftol=1e-90, return_both=false,
                      for_s_rates=true, QNM=false, QNM_E=1.0,
                      erg_Guess=nothing, max_n_qnm=10)
    error("This function must be called from solve_sr_rates context")
end

"""
    find_im_zero(mu, M, n, l, m; kwargs) -> zero_spin::Float64

Find spin at which imaginary eigenvalue component crosses zero.

For a given boson mass and quantum numbers, finds the black hole spin
at which the eigenvalue transitions from stable to unstable (or vice versa).

**Arguments**:
- `mu::Float64` - Boson mass
- `M::Float64` - Black hole mass
- `n, l, m::Int` - Quantum numbers

**Returns**:
- `zero_spin::Float64` - Spin at which Im(ω) = 0

**Physical Significance**:
This spin value represents the critical point for superradiance onset.
Below this spin, modes are stable; above, they grow superradiantly.
"""
function find_im_zero(mu, M, n, l, m; debug=false, Ntot=3000,
                      iter=1000, xtol=1e-6, ftol=1e-20)
    error("This function must be called from solve_sr_rates context")
end

"""
    compute_gridded(mu, M, a, n, l, m; kwargs) -> eigenvalues

Compute eigenvalues on a grid of spin parameters.

Efficiently computes eigenvalues across a range of black hole spins
for parameter space exploration and mapping.

**Arguments**:
- `mu::Float64` - Boson mass
- `M::Float64` - Black hole mass
- `a::Float64` - Reference spin parameter
- `n, l, m::Int` - Quantum numbers

**Keyword Arguments**:
- `Ntot::Int` - Total grid points (default: 2000)
- `iter::Int` - Iterations per point (default: 50)
- `xtol::Float64` - Tolerance (default: 1e-7)
- `npts::Int` - Number of spin points (default: 30)
- `amin::Float64` - Minimum spin (default: 0.0)
- `compute_neg::Bool` - Compute negative spin (default: false)

**Returns**:
- `eigenvalues::Vector` - Eigenvalues at each spin point
"""
function compute_gridded(mu, M, a, n, l, m; Ntot=2000, iter=50,
                        xtol=1e-7, npts=30, amin=0.0, compute_neg=false)
    error("This function must be called from solve_sr_rates context")
end

"""
    radial_inf(erg, mu, M, a, l, m; kwargs) -> wf_inf

Solve radial equation at spatial infinity.

Computes the radial wave function behavior at large radii using
asymptotic expansions and matching conditions.

**Arguments**:
- `erg::Complex` - Eigenvalue
- `mu::Float64` - Boson mass
- `M::Float64` - Black hole mass
- `a::Float64` - Black hole spin
- `l::Int` - Orbital angular momentum
- `m::Int` - Azimuthal quantum number

**Returns**:
- `wf_inf::Vector` - Wave function at infinity
"""
function radial_inf(erg, mu, M, a, l, m; rpts=1000, rmax_val=1e4,
                    debug=false, iter=50, xtol=1e-120, ftol=1e-120,
                    sve_for_test=false, fnm="test_store/test_radial",
                    is_infin=true)
    error("This function must be called from solve_sr_rates context")
end

"""
    spheroidals(l, m, a, erg) -> coefficients

Compute spheroidal harmonic expansion coefficients.

Calculates the angular eigenvalue and expansion coefficients for
spheroidal harmonics in Kerr geometry.

**Arguments**:
- `l::Int` - Orbital angular momentum
- `m::Int` - Azimuthal quantum number
- `a::Float64` - Black hole spin
- `erg::Complex` - Eigenvalue (frequency)

**Returns**:
- `coefficients` - Spheroidal harmonic coefficients
"""
function spheroidals(l, m, a, erg)
    error("This function must be called from solve_sr_rates context")
end

end # module EigenvalueComputation
