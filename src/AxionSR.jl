"""
    AxionSR

Main module for Axion Superradiance simulations and analysis.

This module provides a unified interface to the complete Axion Superradiance
research toolkit, including:
- Black hole superradiance dynamics
- Axion cloud evolution and timescale calculations
- Eigenvalue computation and analysis
- Rate coefficient calculations
- Visualization and analysis tools

**Key Components**:

- `super_rad_check()` - Main interface for superradiance simulations
- `solve_system()` - Time evolution of black hole and axion cloud
- Rate coefficients and interpolation functions
- Eigenvalue solvers for stability analysis
- Statistical analysis tools

**Quick Start**:
```julia
using AxionSR

# Run a superradiance simulation
final_spin, final_mass = super_rad_check(
    M_BH=1.0,           # Solar masses
    aBH=0.8,            # Black hole spin
    massB=1e-20,        # Axion mass
    f_a=1e16,           # Axion decay constant
    tau_max=100.0       # Evolution time in Gyrs
)
```

**Modules**:

- `Constants` - Physical and numerical constants
- `Core.EvolutionHelpers` - Boundary condition and evolution utilities
- `Core.RateCoefficients` - Type-safe rate coefficient management
- `Core.LoadRatesStructured` - Backward compatibility for rate loading
- `Numerics.RateComputation` - Rate interpolation and computation
- `Numerics.EigenvalueComputation` - Eigenvalue solver APIs
- `Numerics.ScatteringRates` - Scattering dynamics

**Main Functions**:

The public API includes:
- `super_rad_check()` - Primary simulation interface
- `solve_system()` - System evolution solver
- `solve_radial()` - Eigenvalue computation
- Rate loading and interpolation functions
- Statistical analysis functions

**Documentation**:

For detailed information on specific functions, use Julia's help system:
```julia
help(super_rad_check)
help(solve_system)
```

See also:
- `README_PROJECT_STATUS.md` - Project overview and phase completion
- `EIGENVALUE_PLOTS_README.md` - Visualization tools guide
- Module documentation in respective source files

**Dependencies**:

- OrdinaryDiffEq.jl - Differential equation solving
- Interpolations.jl - Interpolation functions
- Statistics.jl - Statistical computations
- DelimitedFiles.jl - File I/O
- NPZ.jl - NumPy-compatible file format
- Plots.jl - Visualization (for plotting functions)
- LaTeXStrings.jl - LaTeX string formatting

**Version**: 1.0
**Last Updated**: November 2025
**Status**: Production Ready
"""

module AxionSR

# Standard library imports
using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Dates
using Interpolations
using Printf

# Core modules
include("Core/constants.jl")
include("Core/evolution_helpers.jl")
include("Core/rate_coefficients.jl")
include("Core/load_rates_structured.jl")

# Numerics modules
include("Numerics/rate_computation.jl")
include("Numerics/eigenvalue_computation.jl")
include("Numerics/scattering_rates.jl")

# Main computation modules
include("solve_sr_rates.jl")
include("load_rates.jl")
include("solve_system_unified.jl")
include("super_rad.jl")
include("stat_analysis.jl")

# Public exports - Main interfaces
export super_rad_check
export solve_system

# Public exports - Rate computation
export RateComputation
export RateCoefficients

# Public exports - Eigenvalue computation
export EigenvalueComputation
export ScatteringRates

# Public exports - Helper functions
export get_clamped_spin
export enforce_bosenova_boundary!
export enforce_all_boundaries!
export check_energy_floor
export is_near_bosenova

# Public exports - Constants
export GNew, c_light, hbar, M_sun, year_to_sec, au_to_m
export minSpin, maxSpin, SOLVER_TOLERANCES

end # module AxionSR
