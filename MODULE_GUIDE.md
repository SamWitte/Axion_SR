# AxionSR Module Usage Guide

**Version**: 1.0
**Status**: Production Ready
**Last Updated**: November 2025

---

## Overview

AxionSR is a comprehensive Julia package for studying black hole superradiance of ultralight bosons (axions). This guide explains how to use the module system and access the various computational tools.

## Installation & Import

### Basic Import

```julia
using AxionSR
```

This imports the main AxionSR module and all its public functions.

### Accessing Sub-modules

```julia
using AxionSR.RateComputation
using AxionSR.RateCoefficients
using AxionSR.EigenvalueComputation
```

## Core Functions

### Main Simulation Interface

#### `super_rad_check()`

The primary interface for running superradiance simulations.

```julia
final_spin, final_mass = super_rad_check(
    M_BH::Float64,              # Black hole mass (solar masses)
    aBH::Float64,               # Black hole spin (0 ≤ a < 1)
    massB::Float64,             # Boson mass (natural units)
    f_a::Float64,               # Axion decay constant
    ;
    tau_max=1e4,                # Evolution time (Gyrs)
    Nmax=3,                     # Quantum level truncation
    spinone=false,              # Spin-1 field mode
    cheby=true,                 # Chebyshev decomposition
    # ... additional parameters
)
```

**Example**:
```julia
using AxionSR

# Simulate 1 solar mass black hole at 80% spin
# with ultra-light boson
final_a, final_M = super_rad_check(
    M_BH=1.0,
    aBH=0.8,
    massB=1e-20,
    f_a=1e16,
    tau_max=100.0
)

println("Final spin: $final_a")
println("Final mass: $final_M")
```

### System Evolution

#### `solve_system()`

Direct access to time-evolution solver with more control.

```julia
final_spin, final_mass = solve_system(
    massB, f_a, aBH, M_BH, tau_max;
    spinone=false,
    non_rel=true,
    # ... other parameters
)
```

## Working with Rate Coefficients

### Structured Rate System (Phase 3.1)

```julia
using AxionSR.RateCoefficients

# Build rate database for Nmax=4
db = build_default_rates(4; non_relativistic=true)

# Get rate for specific transition
rate = get_rate(db, "211_211^GW")

# List all available rates
rates = list_rates(db)

# Evaluate rate with parameters
value = evaluate_rate(rate, alpha_parameter)
```

### Legacy Dictionary System

```julia
using AxionSR.Core.LoadRatesStructured

# Load rates in dictionary format (backward compatible)
rate_dict = load_rate_coeffs_structured()
```

## Eigenvalue Computation

### Direct Eigenvalue Solving

```julia
using AxionSR.EigenvalueComputation

# Compute eigenvalue for quantum state (n,l,m)
_, _, eigenvalue = solve_radial(
    mu,        # Boson mass
    M_BH,      # Black hole mass
    a,         # Black hole spin
    n, l, m,   # Quantum numbers
    ;
    return_erg=true
)

# Extract just imaginary part (stability indicator)
im_part = find_im_part(mu, M_BH, a, n, l, m)

# Find threshold for superradiance
critical_spin = find_im_zero(mu, M_BH, n, l, m)
```

### Rate Computation

```julia
using AxionSR.RateComputation

# Compute rates with interpolation
quantum_configs = [(2,1,1,0.01), (3,2,1,0.02), (3,2,2,0.03)]
rates, funcs, dict = compute_sr_rates(
    quantum_configs,
    M_BH, aBH, alpha
)

# Use interpolation function for custom spin
func = dict[(n, l, m)]
rate_value = func(new_spin)
```

## Constants & Physical Parameters

```julia
using AxionSR

# Physical constants
GNew          # Gravitational constant
c_light       # Speed of light
hbar          # Reduced Planck constant
M_sun         # Solar mass
year_to_sec   # Year to seconds conversion
au_to_m       # AU to meters conversion

# Numerical constraints
minSpin       # Minimum allowed spin
maxSpin       # Maximum allowed spin (< 1)

# Solver tolerances
SOLVER_TOLERANCES  # Dict of numerical tolerances
```

## Working with Evolution Helpers

```julia
using AxionSR

# Clamp spin to allowed range
clamped_a = get_clamped_spin(a)

# Apply boundary conditions during evolution
enforce_bosenova_boundary!(mass, spin)
enforce_all_boundaries!(M, a, E)

# Check special conditions
if is_near_bosenova(M, a)
    println("Warning: Near bosenova")
end

if check_energy_floor(E)
    println("Warning: Energy too low")
end
```

## Test Suites

### Running All Tests

```bash
cd /path/to/AxionSR
julia test/runtests.jl
```

### Running Individual Test Suites

```bash
# Unit tests
julia test/unit_tests.jl

# Integration tests
julia test/integration_tests.jl

# Rate coefficient tests
julia test/coefficient_tests.jl

# Computation tests
julia test/computation_tests.jl
```

## Visualization & Analysis Scripts

### Eigenvalue Visualization

```bash
# Quick test (10-20 minutes)
julia scripts/plot_eigenvalue_quick.jl

# Full analysis (2-4 hours)
julia scripts/plot_eigenvalue_imaginary.jl
```

See `scripts/README.md` for details.

## Module Architecture

### Directory Structure

```
src/
├── AxionSR.jl              # Main module
├── Constants.jl            # Physical constants
├── Core/
│   ├── evolution_helpers.jl
│   ├── rate_coefficients.jl
│   └── load_rates_structured.jl
├── Numerics/
│   ├── rate_computation.jl
│   ├── eigenvalue_computation.jl
│   └── scattering_rates.jl
├── super_rad.jl            # Main simulation logic
├── solve_sr_rates.jl       # Rate calculations
├── solve_system_unified.jl # Evolution solver
├── load_rates.jl           # Rate loading
└── stat_analysis.jl        # Statistical tools
```

## Common Workflows

### Workflow 1: Quick Superradiance Check

```julia
using AxionSR

# Simulate with default parameters
final_spin, final_mass = super_rad_check(
    M_BH=1.0,
    aBH=0.8,
    massB=1e-20,
    f_a=1e16
)
```

### Workflow 2: Rate Coefficient Analysis

```julia
using AxionSR.RateCoefficients

# Build database for different Nmax values
for nmax in 3:5
    db = build_default_rates(nmax)
    n_rates = count_rates(db)
    println("Nmax=$nmax: $n_rates rates")
end
```

### Workflow 3: Eigenvalue Analysis

```julia
using AxionSR.EigenvalueComputation

# Find superradiance thresholds
for (n, l, m) in [(2,1,1), (3,2,1), (3,2,2)]
    threshold = find_im_zero(1e-20, 1.0, n, l, m)
    println("($n,$l,$m): Critical spin = $threshold")
end
```

## Documentation & Resources

- `README_PROJECT_STATUS.md` - Project overview
- `EIGENVALUE_PLOTS_README.md` - Visualization guide
- Module docstrings available via Julia `help()`

---

**Module Version**: 1.0
**Production Status**: ✅ Ready for use
**Last Updated**: November 2025
