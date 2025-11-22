"""
    scattering_rates.jl

Scattering rate computation module for black hole-axion interactions.

This module provides functions for computing scattering rates between different
quantum states of axions around rotating black holes, including:
- Bound state to bound state scattering
- Bound state to continuum scattering
- Frequency shift computations

**Main Functions**:
- `sr_rates()` - Primary scattering rate calculator
- `s_rate_bnd()` - Bound state scattering rates
- `s_rate_inf()` - Continuum scattering rates
- `freq_shifts()` - Frequency shift calculations

**Dependencies**:
- Eigenvalue functions from solve_sr_rates.jl
- Integration routines for wave function overlaps
- Special functions (Bessel, hypergeometric, etc.)

**Physical Background**:
Scattering rates quantify the rate at which an axion transitions from one
quantum state to another. Key processes include:
- Superradiant transitions (extracting energy from black hole)
- Decay transitions (axion loses energy)
- Bound state oscillations

**Rate Formulas**:
The transition rates between states are computed from:
- Wave function overlaps at different radii
- Integration over scattering phase space
- Special function evaluations for angular momentum couplings
"""

module ScatteringRates

export sr_rates, s_rate_bnd, s_rate_inf, freq_shifts

"""
    sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true) -> rate

Compute superradiance scattering rate between specified quantum states.

Primary interface for computing the rate at which axions scatter between
different bound states around the black hole, accounting for superradiance.

**Arguments**:
- `n::Int` - Principal quantum number
- `l::Int` - Orbital angular momentum
- `m::Int` - Azimuthal quantum number
- `massB::Float64` - Boson mass (axion mass)
- `MBH::Float64` - Black hole mass
- `aBH::Float64` - Black hole spin parameter (dimensionless, 0 ≤ a < 1)

**Keyword Arguments**:
- `impose_low_cut::Float64` - Minimum allowed coupling strength (default: 0.001)
- `solve_322::Bool` - Special handling for (3,2,2) state (default: true)

**Returns**:
- `rate::Float64` - Transition rate in appropriate units

**Physical Significance**:
The returned rate determines how quickly the axion cloud evolves:
- Larger rate: Faster transition between states
- Rate determines characteristic timescales for superradiance

**Special Cases**:
- The (3,2,2) state requires special numerical treatment due to near-zero modes
- Rates vanish above certain couplings (alpha_threshold)
"""
function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true)
    error("This function must be called from solve_sr_rates context")
end

"""
    s_rate_bnd(n1, l1, m1, n2, l2, m2, n3, l3, m3; kwargs) -> rate

Compute scattering rate for bound-bound transitions.

Calculates the transition rate from state (n1,l1,m1) to state (n2,l2,m2)
with photon/graviton emission/absorption in mode (n3,l3,m3).

**Arguments**:
- `n1, l1, m1::Int` - Initial quantum state
- `n2, l2, m2::Int` - Final quantum state
- `n3, l3, m3::Int` - Emitted/absorbed quantum
- Plus: mu, M, a parameters

**Keyword Arguments**:
- `kpts::Int` - Momentum grid points (default: 10)
- `rpts::Int` - Radial grid points (default: 2000)
- `rmaxT::Float64` - Maximum radius (default: 100)
- `inf_nr::Bool` - Include infinite radius contributions (default: true)
- `Nang::Int` - Angular integration points (default: 100000)
- `Npts_Bnd::Int` - Boundary points (default: 1000)
- `debug::Bool` - Verbose output (default: false)
- `include_cont::Bool` - Include continuum (default: true)
- `Ntot_safe::Int` - Safety limit (default: 5000)
- `sve_for_test::Bool` - Save for testing (default: false)
- `bnd_thresh::Float64` - Boundary threshold (default: 1e-3)
- `use_analytic::Bool` - Use analytic formulas where possible (default: false)
- `eps_r::Float64` - Radial epsilon (default: 1e-10)
- `NON_REL::Bool` - Non-relativistic approximation (default: true)

**Returns**:
- `rate::Float64` - Bound state scattering rate

**Computation Method**:
- Computes wave function overlaps
- Integrates over intermediate states
- Includes QED/field theory corrections
"""
function s_rate_bnd(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3;
                    kpts=10, rpts=2000, rmaxT=100, inf_nr=true, Nang=100000,
                    Npts_Bnd=1000, debug=false, include_cont=true, Ntot_safe=5000,
                    sve_for_test=false, bnd_thresh=1e-3, use_analytic=false,
                    eps_r=1e-10, NON_REL=true)
    error("This function must be called from solve_sr_rates context")
end

"""
    s_rate_inf(n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min; kwargs) -> rate

Compute scattering rate to continuum states.

Calculates the rate at which axions scatter from a bound state into the
continuum (plane waves at infinity).

**Arguments**:
- `n1, l1, m1::Int` - Initial bound state quantum numbers
- `n2, l2, m2::Int` - Final state quantum numbers (may be continuum)
- `n3, l3, m3::Int` - Intermediate state quantum numbers
- `lF_min::Int` - Minimum final orbital angular momentum
- Plus: mu, M, a parameters

**Keyword Arguments**:
- `rpts::Int` - Radial grid resolution (default: 4000)
- `rmaxT::Float64` - Maximum integration radius (default: 90)
- `sve_for_test::Bool` - Save intermediate results (default: false)
- `inf_nr::Bool` - Infinite radius non-relativistic (default: false)
- `Npts_Bnd::Int` - Boundary integration points (default: 1000)
- `Nang::Int` - Angular integration points (default: 300000)
- `debug::Bool` - Verbose output (default: false)
- `Ntot_safe::Int` - Safety limit for iteration (default: 2000)
- `xtol::Float64` - Position tolerance (default: 1e-2)
- `ftol::Float64` - Function tolerance (default: 1e-2)
- `iter::Int` - Maximum iterations (default: 20)

**Returns**:
- `rate::Float64` - Continuum scattering rate

**Physical Context**:
These rates describe ionization-like processes where the axion escapes
the black hole's potential well into the continuum.
"""
function s_rate_inf(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min;
                    rpts=4000, rmaxT=90, sve_for_test=false, inf_nr=false,
                    Npts_Bnd=1000, Nang=300000, debug=false, Ntot_safe=2000,
                    xtol=1e-2, ftol=1e-2, iter=20)
    error("This function must be called from solve_sr_rates context")
end

"""
    freq_shifts(mu, M, a, n1, l1, m1, n2, l2, m2; kwargs) -> shifts

Compute frequency shifts due to interactions.

Calculates energy shifts and frequency corrections to bound state energies
due to interactions with other quantum modes.

**Arguments**:
- `mu::Float64` - Boson mass (axion mass)
- `M::Float64` - Black hole mass
- `a::Float64` - Black hole spin parameter
- `n1, l1, m1::Int` - State 1 quantum numbers
- `n2, l2, m2::Int` - State 2 quantum numbers

**Keyword Arguments**:
- `rpts::Int` - Radial grid points (default: 500)
- `rmaxT::Float64` - Maximum radius (default: 100)
- `Nang::Int` - Angular integration points (default: 100000)
- `Npts_Bnd::Int` - Boundary points (default: 500)
- `epsil_2::Float64` - Energy epsilon (default: 1.0)

**Returns**:
- `shifts` - Frequency shift corrections

**Applications**:
- Self-energy calculations
- Many-body effects in axion clouds
- Precision spectroscopy predictions
"""
function freq_shifts(mu, M, a, n1, l1, m1, n2, l2, m2;
                     rpts=500, rmaxT=100, Nang=100000, epsil_2=1.0,
                     Npts_Bnd=500)
    error("This function must be called from solve_sr_rates context")
end

end # module ScatteringRates
