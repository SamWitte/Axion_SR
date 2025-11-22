# Code Cleanup & Refactoring Guide

This document provides a detailed, prioritized roadmap for cleaning up the Axion_SR codebase.

## New Directory Structure (Created)

```
src/
├── Core/
│   ├── constants.jl         # Physical constants and tolerances ✓
│   ├── parameters.jl        # Parameter structures (QMState, EvolutionParameters, etc) ✓
│   ├── physics.jl           # TODO: Core physics calculations
│   └── evolution.jl         # TODO: Merged ODE solving code
├── Numerics/
│   ├── ode_solvers.jl       # TODO: ODE callbacks and setup
│   ├── heun_equations.jl    # TODO: Move heunc.jl here (already clean)
│   └── interpolation.jl     # TODO: Interpolation utilities
├── DataIO/
│   ├── loaders.jl           # TODO: Rate loading functions
│   ├── rate_coefficients.jl # TODO: Structured rate data
│   └── serialization.jl     # TODO: File I/O utilities
├── Statistics/
│   ├── mcmc.jl              # TODO: MCMC sampler
│   ├── likelihoods.jl       # TODO: Likelihood functions
│   └── priors.jl            # TODO: Prior definitions
├── Models/
│   ├── tde.jl               # TODO: Tidal disruption events
│   └── bh_models.jl         # TODO: Black hole property models
└── AxionSR.jl              # TODO: Main module file with include statements
scripts/                     # TODO: Executable scripts
├── run_mcmc.jl
├── compute_rates.jl
├── quick_limits.jl
└── ...
tests/                       # TODO: Proper test suite
├── test_physics.jl
├── test_evolution.jl
└── ...
```

## Phase 1: Critical Cleanups (High Impact, Lower Effort)

### 1.1 Remove Debug Print Statements (1-2 hours)

**Files affected**: `super_rad.jl`, `stat_analysis.jl`, others

**Before**:
```julia
if debug
    print("Alpha \t", alph, "\n")
    print("Solving system... \n")
end
```

**After** (create `src/Utils/logging.jl`):
```julia
using Logging

function debug_log(logger::AbstractLogger, msg::String, debug_flag::Bool)
    if debug_flag
        @debug msg
    end
end
```

**Action items**:
- [ ] Create `src/Utils/logging.jl` with proper logging setup
- [ ] Replace all `if debug; print(...)` blocks with `debug_log()` calls
- [ ] Use Julia's `Logging` module configuration instead of print statements
- [ ] Search for patterns: `print(`, `println(`, `@time` in production code

**Files to update**:
- `super_rad.jl` (26 print statements)
- `stat_analysis.jl` (4 print statements)
- `solve_sr_rates.jl` (check for any)
- Any others with debug output

### 1.2 Remove All Commented-Out Code (30 minutes - 1 hour)

**Files affected**: `solve_sr_rates.jl`, `sr_matter_ints.jl`, `super_rad.jl`

**Action items**:
- [ ] Search and remove all commented code blocks (50+ lines)
- [ ] In `solve_sr_rates.jl`: Remove commented lines 30-31, 41, 62, 183-184, 199-200, 211-212
- [ ] In `sr_matter_ints.jl`: Remove 50+ commented-out lines
- [ ] In `super_rad.jl`: Remove commented test code at line 558-560
- [ ] Verify with git that removed code can be recovered if needed

**Commands**:
```bash
# Find commented lines (rough check)
grep -n "^\s*#" src/solve_sr_rates.jl | head -20
grep -n "^\s*#" src/sr_matter_ints.jl | head -30
```

### 1.3 Fix Obvious Bugs (15-30 minutes)

**Bug 1**: `solve_sr_rates.jl:262` - Floating point equality
```julia
# BEFORE (wrong):
if erg_4G != erg_ind  # Floating point comparison!

# AFTER:
const FP_TOL = 1e-10
if abs(erg_4G - erg_ind) > FP_TOL
```

**Bug 2**: Unexplained magic number at `super_rad.jl:298`
```julia
# BEFORE:
u_fake = u_real * 1.1 # was 2!

# AFTER:
# Growth estimate with 10% safety margin (previously used 2x, too aggressive)
GROWTH_SAFETY_MARGIN = 1.1
u_fake = u_real * GROWTH_SAFETY_MARGIN
```

## Phase 2: Structure & Duplication (Medium Impact, Medium Effort)

### 2.1 Merge `solve_system` and `solve_system_spinone` (3-4 hours)

**Issue**: ~200 lines of duplicate code across two functions

**Strategy**:
1. Create unified `solve_system()` function with parameter `spinone::Bool`
2. Extract duplicate code into helpers:
   - `initialize_quantum_levels()` - Shared setup
   - `enforce_boundary_conditions!()` - Repeated boundary checking
   - `enforce_spin_bounds()` - Repeated spin clamping
   - `setup_ode_problem()` - Shared ODE setup

**Before** (solve_system):
```julia
function solve_system(mu, fa, aBH, M_BH, t_max; ..., spinone=false)
    # ... 400 lines ...
end

function solve_system_spinone(mu, aBH, M_BH, tau_max; ...)
    # ... 200 lines of similar code ...
end
```

**After** (create `src/Core/evolution.jl`):
```julia
function solve_system(
    params::EvolutionParameters,
    settings::SolverSettings;
    spinone::Bool = false
)
    # Single unified implementation
    if spinone
        setup = setup_spinone_system(params)
    else
        setup = setup_standard_system(params)
    end
    # ... shared evolution code ...
end
```

**Action items**:
- [ ] Extract `enforce_boundary_conditions!()` helper (handles both cases)
- [ ] Extract `enforce_spin_bounds()` helper
- [ ] Extract `initialize_quantum_levels()` helper
- [ ] Merge the two functions
- [ ] Test with both spinone=true and spinone=false
- [ ] Update references in `super_rad.jl`

### 2.2 Extract Helper Functions from `super_rad.jl` (2-3 hours)

**Pattern 1**: Repeated boundary enforcement (appears 4+ times)
```julia
# Extract to src/Core/evolution.jl
function enforce_bounds!(u, u_real, bounds_info::BoundaryInfo)
    for i in 1:length(u)
        clamp_to_bounds!(u, u_real, i, bounds_info)
    end
end
```

**Pattern 2**: Repeated spin bound checking (appears 3+ times)
```julia
# Extract to src/Core/evolution.jl
function get_clamped_spin(spin_val::Float64)::Tuple{Float64, Float64}
    """Returns (clamped_spin, rP_factor)"""
    if spin_val > maxSpin
        return (maxSpin, 1.0 + sqrt(1 - maxSpin^2))
    elseif spin_val < 0.0
        return (0.0, 2.0)
    else
        return (spin_val, 1.0 + sqrt(1 - spin_val^2))
    end
end
```

**Pattern 3**: Parameter unpacking (repeated in RHS_ax! functions)
```julia
# Create src/Core/parameters.jl helper
function unpack_evolution_params(packed::EvolutionParameters)
    return packed.mu, packed.f_a, packed.e_max_2, packed.spin_initial, packed.mass_initial
end
```

**Action items**:
- [ ] Create `src/Core/evolution_helpers.jl` or add to `evolution.jl`
- [ ] Extract all 5-7 repeated patterns
- [ ] Replace calls in `super_rad.jl` with helper function calls
- [ ] Test that behavior is unchanged

### 2.3 Consolidate Tolerance Constants (1 hour)

**Before**: Magic numbers scattered throughout
```julia
1e-2   # Line 205, 261, 301, 356, 384
1e-75  # Line 311
1e-30  # ODE tolerance
1e-100 # eq_threshold
```

**After** (Already started in `src/Core/constants.jl`):
```julia
const SOLVER_TOLERANCES = (
    bosenova_threshold = 1e-2,      # Binding energy boundary detection
    energy_floor = 1e-75,            # Underflow prevention
    ode_absolute_tight = 1e-30,      # Sensitive integration
    ode_absolute_default = 1e-10,
    ode_relative_default = 1e-5,
    ode_relative_high = 1e-3,
    spin_mismatch = 1e-2,            # Spin boundary tolerance
)

# Usage in super_rad.jl:
if abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold
    # Boundary condition enforcement
end
```

**Action items**:
- [ ] Move all magic number constants to `src/Core/constants.jl`
- [ ] Update all references to use named constants
- [ ] Add documentation comments explaining physical/numerical meaning
- [ ] Files to update: `super_rad.jl`, `solve_sr_rates.jl`, `stat_analysis.jl`

## Phase 3: Architecture & Organization (High Impact, Higher Effort)

### 3.1 Extract Rate Coefficients (2-3 hours)

**Current problem** (`load_rates.jl` lines 27-156):
```julia
# 130 lines of hardcoded coefficients with no explanation
Drate["211_211^322^BH"] = 4.2e-7 .* alph^11 .* faFac * rP
Drate["322_322^211^Inf"] = 1.1e-8 * alph^8 .* faFac
# ... 128 more lines ...
```

**Solution**: Create structured data file `src/DataIO/rate_coefficients.jl`

```julia
"""
Transition rate coefficients for axion superradiance.
Format: Dict mapping (initial_state, final_state, destination) -> (coefficient, power_of_alpha)
"""
const RATE_COEFFICIENTS = Dict(
    ("211", "322", "BH") => (coeff = 4.2e-7, alpha_power = 11, has_spin_factor = true),
    ("322", "211", "Inf") => (coeff = 1.1e-8, alpha_power = 8, has_spin_factor = false),
    # ... etc ...
)

function compute_rate(initial::String, final::String, dest::String, alpha::Float64, f_a::Float64, spin::Float64)
    """Compute transition rate from stored coefficients"""
    key = (initial, final, dest)
    haskey(RATE_COEFFICIENTS, key) || error("Unknown transition: $key")

    entry = RATE_COEFFICIENTS[key]
    rate = entry.coeff * alpha^entry.alpha_power

    if entry.has_spin_factor
        rate *= (1.0 + sqrt(1.0 - spin^2))
    end

    return rate
end
```

**Action items**:
- [ ] Create `src/DataIO/rate_coefficients.jl` with structured rate data
- [ ] Extract all 130+ hardcoded rates from `load_rates.jl` into structured format
- [ ] Create `compute_rate()` helper function
- [ ] Rewrite `load_rates.jl` to use structured data
- [ ] Document the physical meaning of each coefficient
- [ ] Consider alternative: Load from external YAML/JSON file for better maintainability

### 3.2 Move Scripts Out of `src/` (1 hour)

**Create `scripts/` directory with**:
- [ ] `scripts/run_mcmc.jl` - Moved from `src/MCMC.jl`
- [ ] `scripts/compute_rates.jl` - Moved from `src/TEST_compute_rates.jl`
- [ ] `scripts/quick_limits.jl` - Moved from `src/estimate_lim.jl`
- [ ] `scripts/analyze_output.jl` - For post-processing MCMC results

**Update each script**:
```julia
# scripts/run_mcmc.jl (instead of src/MCMC.jl)
import Pkg; Pkg.activate(".")

# Include modules from src/ instead of execute inline
include("../src/AxionSR.jl")
using .AxionSR

# Main script logic
function main()
    # Parse arguments, run MCMC, save results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
```

### 3.3 Separate `TEST.jl` into Proper Test Suite (2-3 hours)

**Current**: 2,463 line file mixing tests and production code

**Action items**:
- [ ] Create `tests/runtests.jl` with standard Julia test structure
- [ ] Separate into:
  - [ ] `tests/test_physics.jl` - Physics calculations
  - [ ] `tests/test_evolution.jl` - ODE solving
  - [ ] `tests/test_rates.jl` - Rate loading and computation
  - [ ] `tests/test_mcmc.jl` - MCMC functionality
- [ ] Use Julia's `Test` package with proper `@test` assertions
- [ ] Document test coverage and expected behavior

**Template**:
```julia
# tests/test_physics.jl
using Test
using AxionSR

@testset "Physics Calculations" begin
    @testset "Energy calculations" begin
        # Tests here
    end

    @testset "Rate computations" begin
        # Tests here
    end
end
```

## Phase 4: Documentation (Lower Impact, Medium Effort)

### 4.1 Add Docstrings to Major Functions (2-3 hours)

**Template**:
```julia
"""
    solve_system(params::EvolutionParameters, settings::SolverSettings; spinone::Bool=false)

Solve the coupled ODE system for axion-black hole superradiance evolution.

# Arguments
- `params::EvolutionParameters`: Evolution parameters (mass, spin, coupling, etc.)
- `settings::SolverSettings`: ODE solver configuration and tolerances
- `spinone::Bool`: Include spin-1 modes if true (default: false)

# Returns
- `Tuple{Float64, Float64}`: Final black hole spin and mass after evolution

# Physical Details
Solves the coupled ODEs for N quantum levels of the axion cloud, tracking:
- Binding energy evolution for each level
- Black hole spin and mass change
- Gravitational wave extraction

The system includes:
- Superradiant growth (exponential amplification)
- Nonlinear saturation via level mixing
- Angular momentum extraction (spin-down)
"""
function solve_system(params::EvolutionParameters, settings::SolverSettings; spinone::Bool=false)
    # Implementation
end
```

**Files to document**:
- [ ] `super_rad_check()` - Main entry point
- [ ] `solve_system()` - ODE solver
- [ ] `compute_sr_rates()` - Rate calculation
- [ ] `log_probability()` - MCMC likelihood
- [ ] All other public functions

### 4.2 Document Physical Constants & Units (1 hour)

Update `src/Core/constants.jl` to clarify all units:

```julia
"""
    const GNew = 7484169213.942707

Gravitational constant in natural units where ℏ = c = 1.
Units: [1/(M_⊙·eV)]

This is G * (M_⊙ to eV conversion) to allow direct computation of dimensionless
coupling strength α = G * M_BH * m_axion in natural units.
"""
const GNew = 7484169213.942707
```

**Action items**:
- [ ] Add docstring to each constant explaining physical meaning
- [ ] Clarify unit system being used
- [ ] Add example calculations

## Phase 5: Validation & Testing (Final Step)

- [ ] Run full test suite after each major phase
- [ ] Compare numerical results before/after refactoring
- [ ] Test on reference black hole cases (GRS 1915+105, Cygnus X-1, GW231123)
- [ ] Verify MCMC chains produce same posteriors
- [ ] Check git diff to ensure no logic changes

## Recommended Execution Order

1. **Day 1**: Phase 1 (Quick wins - Phases 1.1-1.3) - 3-4 hours
   - Remove debug prints
   - Remove commented code
   - Fix obvious bugs

2. **Day 2-3**: Phase 2 (Structure - Phases 2.1-2.3) - 8-10 hours
   - Merge duplicate functions
   - Extract helpers
   - Consolidate constants

3. **Day 4**: Phase 3 Part 1 (Architecture - Phases 3.1-3.2) - 3-4 hours
   - Extract rate coefficients
   - Move scripts

4. **Day 5**: Phase 3 Part 2 + Phase 4 (Organization & Docs - Phases 3.3-4.2) - 6-8 hours
   - Separate tests
   - Add docstrings
   - Document constants

5. **Day 6**: Phase 5 (Validation) - 4-6 hours
   - Run tests
   - Verify numerical results
   - Final cleanup

**Total estimated time**: 25-35 hours of focused development

## Notes

- Use git branches: `git checkout -b cleanup/remove-debug-prints`, etc.
- Commit frequently with clear messages: "Remove debug prints from super_rad.jl"
- Test after each major change
- Keep a parallel version working until refactoring is complete
- Consider pair programming or code review for critical sections
