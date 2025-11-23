# Axion Superradiance: Black Hole Evolution & Constraints

A computational framework for modeling axion superradiance in rotating black holes and constraining axion physics using multi-messenger astrophysical observations.

### Key Physics

**Superradiance** is an instability where bosonic fields (like axions) extract rotational energy from spinning black holes, causing the black hole to spin down while the axion cloud grows exponentially. This project:

- Computes the eigenvalues of quasi-bound states, and the rate of energy transfer induced by self-interactions 
- Simulates black hole mass and spin evolution
- Compares predictions against observed black hole populations (stellar binaries, AGN, gravitational wave sources)
- Derives constraints on axion mass and coupling strength

## Project Structure

```
src/
├── Core Physics Computation
│   ├── solve_sr_rates.jl          # Primary rate calculation engine
│   ├── heunc.jl                   # Heun equation solver (used in radial wave function solutions)
│   ├── super_rad.jl               # Black hole evolution dynamics
│   └── Constants.jl               # Physical constants & utilities
│
├── Rate Computation & Storage
│   ├── Compute_all_rates.jl       # CLI for computing transition rates
│   ├── load_rates.jl              # Load pre-computed rate tables for evaluation during evolution
│   └── rate_sve/                  # Pre-computed rate files [must be pre-computed!]
│
├── Statistical Analysis & Inference
│   ├── MCMC.jl                    # Affine-invariant MCMC sampler
│   ├── stat_analysis.jl           # Likelihood & statistical functions
│   ├── stat_1dL.jl                # 1D likelihood analysis
│   └── output_mcmc/               # MCMC chains and posterior output
│
├── Gravitational Wave Analysis
│   ├── Ligo_limit.jl              # GW event constraint derivation
│   ├── Ligo_1dL.jl                # 1D GW likelihood analysis
│   └── Ligo_1dL_cosma.jl          # HPC-optimized version
│
├── Auxiliary Computation
│   ├── grid_sr_growth.jl          # Superradiance growth grids
│   ├── spin_down_plot_files.jl    # Spin evolution visualization data
│   ├── Make_Regge_Plt_data.jl     # Regge diagram generation
│   ├── binary_mixing_test.jl      # Quantum level mixing validation
│   └── tde_input.jl               # Tidal disruption event modeling
│
├── Utilities & Testing
│   ├── install_pkgs.jl            # Dependency installation
│   ├── test_sr_rates.jl           # Unit tests
│   └── estimate_lim.jl            # Quick limit estimation
│
├── Data & Results
│   ├── BH_data/                   # Astrophysical black hole parameters
│   ├── stored_limits/             # Compiled constraint results
│   └── plts/                      # Jupyter notebooks & visualizations
│
└── Python Utilities
    ├── Print_all_levels.py        # Generate valid quantum transitions [useful for identifying relevant scattering permutations]
    ├── latex_line_gen.py          # LaTeX output generation [for writing papers]
    ├── Run_Ligo.py                # GW data processing
    └── Combine_data.py            # Data combination utilities
```

## Main Components

### 1. **Quantum Superradiance Physics** (`solve_sr_rates.jl`, `heunc.jl`)

Computes superradiant transition rates for axion-black hole systems:

- **Quantum states**: |nlm⟩ representation (principal, orbital, azimuthal quantum numbers)
- **Transition modes**: Spontaneous emission to lower quantum levels
- **Three-body processes**: Axion-axion scattering enabling higher-level excitations
- **Numerical solvers**: Leaver eigenvalue method + Heun ODE solver for solving Teukolsky equations in Kerr spacetime
- **Supported regime**: Up to Nmax=8 principal quantum levels, axion coupling α from 10⁻⁴ to 10⁻¹

### 2. **Black Hole Evolution Dynamics** (`super_rad.jl`)

Models coupled evolution of black hole properties under axion extraction:

- **Coupled ODE system**: Tracks M(t), a(t), and axion cloud populations simultaneously
- **Physical processes**:
  - Superradiant growth (exponential amplification phase)
  - Nonlinear saturation via level mixing
  - Angular momentum extraction (spin-down)
  - Mass loss via gravitational wave emission
- **Time integration**: Adaptive Runge-Kutta with automatic differentiation support
- **Evolution scales**: up to 10⁸ years of black hole lifetime

### 3. **Statistical Inference** (`MCMC.jl`, `stat_analysis.jl`)

Bayesian parameter estimation to constrain axion properties:

- **Observable** **constraints**:
  - LIGO gravitational wave spin measurements
  - Stellar binary black hole mass/spin (Cygnus X-1, GRS 1915+105, etc.)
  - AGN supermassive black hole spins
  - Tidal disruption event timescales
- **MCMC method**: Affine-invariant ensemble sampler with 10+ parallel chains
- **Parameters**: Axion mass (m_a), axion decay constant (f_a)
- **Output**: Posterior distributions, 95% confidence intervals, exclusion regions

### 4. **Pre-Computed Rate Tables** (`rate_sve/`)

Extensive database of pre-calculated transition rates:

- **740+ rate files** spanning parameter space (axion coupling, black hole spin)
- **File format**: `n1l1m1_n2l2m2_n3l3m3_Destination.dat` (nlm indices + GW destination)
- **Companion data**: Imaginary eigenfrequency data (`.npz` format) for growth rate calculations
- **Destinations**: "BH" (energy stays near black hole) or "Inf" (radiates to infinity)

### 5. **Astrophysical Data** (`BH_data/`)

Real black hole observations:

- **Stellar-mass**: Cygnus X-1, GRS 1915+105, M33 X-7, LMC X-3
- **Gravitational waves**: GW231123, LIGO test events
- **Tidal disruption events**: Star-black hole encounters
- **Parameters**: Mass, spin, age, chirp mass, orbital period
 
## Core Workflow

```
Input: Black hole parameters (M, a, age)
   ↓
Compute transition rates (or load from pre-computed tables)
   ↓
Simulate BH + axion cloud evolution over time
   ↓
Compare final predicted spin vs. observed spin
   ↓
MCMC: Sample (m_a, f_a) parameter space
   ↓
Derive posterior distributions & exclusion limits
```

## Usage

## Dependencies

- **Julia** (v1.7+)
  - DifferentialEquations.jl (ODE solving)
  - SpecialFunctions.jl (Hypergeometric, spheroidal functions)
  - Distributions.jl (Statistical distributions)
  - StatsBase.jl (MCMC diagnostics)
  - ForwardDiff.jl (Automatic differentiation)
  - HDF5.jl, NPZ.jl (Data I/O)

- **Python** (v3.7+)
  - numpy, scipy (numerical computation)
  - matplotlib, seaborn (visualization)

Install Julia dependencies:
```bash
julia install_pkgs.jl
```

## Output & Results

### MCMC Results (`output_mcmc/`)
- MCMC chains (samples from posterior distribution)
- KDE probability density estimates
- Marginalized 1D/2D posterior distributions
- Constraint contours for publication

### Computed Rates (`rate_sve/`)
- Precomputed transition rate databases
- Eigenfrequency tables
- Interpolation data for efficient simulation

### Constraints (`stored_limits/`)
- Compiled axion mass limits
- 95% confidence intervals
- Multi-messenger constraint compilation

**Last Updated**: November 2025
