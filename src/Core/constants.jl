"""
    Core physical constants and numerical tolerances for axion superradiance calculations.
"""

using ForwardDiff: Dual
using DelimitedFiles
using Interpolations

# Physical Constants (SI-compatible units)
const c_km = 2.99792e5              # Speed of light (km/s)
const c_cm = 2.99792e10             # Speed of light (cm/s)
const hbar = 6.582119e-16           # Reduced Planck constant (eV·s)
const GNew = 7484169213.942707      # Gravitational constant (1/(M_⊙·eV))
const M_to_eV = 1.1219176708765e+66 # Solar mass to eV conversion
const M_pl = 1.22e19                # Planck mass (GeV)
const m_elec = 5.11e5               # Electron mass (eV)
const G_to_eV2 = 1.95e-2            # G in eV² units
const e_charge = 0.3                # Elementary charge
const a_fs = 1.0 / 137.0            # Fine structure constant
const θR = 0.0                      # Strong CP phase

# Black hole spin bounds
const maxSpin = 0.998               # Maximum allowed black hole spin
const minSpin = 0.01                # Minimum allowed black hole spin

# Numerical Solver Tolerances
"""
Tolerance parameters for ODE solvers and numerical algorithms.
These values control convergence, boundary condition enforcement, and energy tracking.
"""
const SOLVER_TOLERANCES = (
    bosenova_threshold = 1e-2,      # Bosenova detection: log-space tolerance for binding energy
    energy_floor = 1e-75,           # Minimum energy threshold to avoid underflow
    ode_absolute_tight = 1e-30,     # Tight absolute tolerance for sensitive regions
    ode_absolute_default = 1e-10,   # Default absolute tolerance
    ode_relative_default = 1e-5,    # Default relative tolerance
    ode_relative_high = 1e-3,       # High-speed integration tolerance
    spin_mismatch = 1e-2,           # Spin boundary condition tolerance
)

# Constants for physical processes
const YEAR_IN_SECONDS = 3.15e7      # Conversion factor: 1 year = 3.15e7 seconds
const err_CS = 1e-5                 # Coulomb sum error threshold

"""
    seed(x::Matrix)

Create vector with Dual numbers for automatic differentiation of 3D vectors.
Constructs partials for gradient computation.
"""
const seed = x -> [
    map(y -> Dual(y, (1., 0., 0.)), x[:, 1])
    map(y -> Dual(y, (0., 1., 0.)), x[:, 2])
    map(y -> Dual(y, (0., 0., 1.)), x[:, 3])
]

"""
    grad(x::Vector{Dual})

Extract gradient vector from Dual numbers (used for AD).
"""
const grad = x -> [
    map(x -> x.partials[1], x)
    map(x -> x.partials[2], x)
    map(x -> x.partials[3], x)
]
