"""
    Core parameter structures for axion superradiance simulations.

Defines structured types replacing tuple-based parameter passing throughout the codebase.
"""

"""
    QMState

Represents a quantum mechanical state with principal, orbital, and azimuthal quantum numbers.

# Fields
- `n::Int`: Principal quantum number (n ≥ 1)
- `l::Int`: Orbital angular momentum quantum number (0 ≤ l < n)
- `m::Int`: Azimuthal quantum number (-l ≤ m ≤ l)
"""
struct QMState
    n::Int
    l::Int
    m::Int

    function QMState(n::Int, l::Int, m::Int)
        @assert n ≥ 1 "Principal quantum number n must be ≥ 1"
        @assert 0 ≤ l < n "Orbital quantum number l must satisfy 0 ≤ l < n"
        @assert abs(m) ≤ l "Azimuthal quantum number |m| must be ≤ l"
        new(n, l, m)
    end
end

"""
    EvolutionParameters

Parameters controlling the ODE system evolution for axion-black hole dynamics.

# Fields
- `mu::Float64`: Axion mass in eV
- `f_a::Float64`: Axion decay constant (PQ scale)
- `e_max_2::Float64`: Second energy level maximum
- `spin_initial::Float64`: Initial black hole spin parameter (0 ≤ a ≤ 1)
- `mass_initial::Float64`: Initial black hole mass (M_⊙)
- `impose_low_cut::Float64`: Low-coupling cutoff for alpha parameter
"""
struct EvolutionParameters
    mu::Float64
    f_a::Float64
    e_max_2::Float64
    spin_initial::Float64
    mass_initial::Float64
    impose_low_cut::Float64

    function EvolutionParameters(mu, f_a, e_max_2, spin_i, mass_i, impose_low_cut)
        @assert mu > 0 "Axion mass must be positive"
        @assert f_a > 0 "Decay constant must be positive"
        @assert 0 ≤ spin_i ≤ 1 "Spin must be in [0, 1]"
        @assert mass_i > 0 "Black hole mass must be positive"
        new(mu, f_a, e_max_2, spin_i, mass_i, impose_low_cut)
    end
end

"""
    SolverSettings

Configuration options for ODE solver and numerical algorithms.

# Fields
- `tau_max::Float64`: Maximum evolution time (years)
- `alpha_max_cut::Float64`: Maximum coupling cutoff
- `impose_low_cut::Float64`: Low-coupling cutoff
- `stop_on_a::Float64`: Stop evolution when spin drops below threshold
- `eq_threshold::Float64`: Energy equilibration threshold
- `abstol::Float64`: Absolute tolerance for ODE solver
- `reltol::Float64`: Relative tolerance for ODE solver
- `use_chebyshev::Bool`: Use Chebyshev polynomials for interpolation
- `include_spinone::Bool`: Include spin-1 modes
- `include_nonrelativistic::Bool`: Include non-relativistic corrections
- `nmax::Int`: Maximum principal quantum number to include
- `n_interp_pts::Int`: Number of interpolation points for spin
- `n_interp_pts_l::Int`: Number of interpolation points for l-levels
- `debug::Bool`: Enable debug output
"""
struct SolverSettings
    tau_max::Float64
    alpha_max_cut::Float64
    impose_low_cut::Float64
    stop_on_a::Float64
    eq_threshold::Float64
    abstol::Float64
    reltol::Float64
    use_chebyshev::Bool
    include_spinone::Bool
    include_nonrelativistic::Bool
    nmax::Int
    n_interp_pts::Int
    n_interp_pts_l::Int
    debug::Bool

    function SolverSettings(;
        tau_max = 1e4,
        alpha_max_cut = 100.0,
        impose_low_cut = 0.01,
        stop_on_a = 0.0,
        eq_threshold = 1e-100,
        abstol = 1e-30,
        reltol = 1e-5,
        use_chebyshev = true,
        include_spinone = false,
        include_nonrelativistic = true,
        nmax = 3,
        n_interp_pts = 100,
        n_interp_pts_l = 100,
        debug = false
    )
        @assert tau_max > 0 "Maximum time must be positive"
        @assert abstol > 0 && reltol > 0 "Tolerances must be positive"
        @assert nmax ≥ 1 "Nmax must be at least 1"
        new(tau_max, alpha_max_cut, impose_low_cut, stop_on_a, eq_threshold,
            abstol, reltol, use_chebyshev, include_spinone, include_nonrelativistic,
            nmax, n_interp_pts, n_interp_pts_l, debug)
    end
end

"""
    BlackHoleProperties

Represents a black hole's observable and dynamical properties.

# Fields
- `mass::Float64`: Mass in solar masses
- `spin::Float64`: Dimensionless spin parameter (0 ≤ a ≤ 1)
- `age::Float64`: Age in years (optional, for timing constraints)
- `name::String`: Identifier/name for the object
"""
struct BlackHoleProperties
    mass::Float64
    spin::Float64
    age::Union{Float64, Nothing}
    name::String

    function BlackHoleProperties(mass, spin; age=nothing, name="")
        @assert mass > 0 "Black hole mass must be positive"
        @assert 0 ≤ spin ≤ 1 "Spin must be in [0, 1]"
        if age !== nothing
            @assert age > 0 "Age must be positive"
        end
        new(mass, spin, age, name)
    end
end

"""
    TransitionRate

Represents a single quantum transition rate (spontaneous or stimulated).

# Fields
- `initial_state::QMState`: Initial quantum state
- `final_state::QMState`: Final quantum state
- `rate::Float64`: Transition rate (1/time)
- `destination::Symbol`: Either :BH (stays near BH) or :Inf (escapes to infinity)
"""
struct TransitionRate
    initial_state::QMState
    final_state::QMState
    rate::Float64
    destination::Symbol  # :BH or :Inf

    function TransitionRate(initial::QMState, final::QMState, rate::Float64, dest::Symbol)
        @assert rate ≥ 0 "Transition rate must be non-negative"
        @assert dest ∈ (:BH, :Inf) "Destination must be :BH or :Inf"
        new(initial, final, rate, dest)
    end
end
