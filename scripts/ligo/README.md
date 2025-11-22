# LIGO Statistical Analysis Scripts

This directory contains scripts for running MCMC-based statistical analysis on LIGO gravitational wave detections to derive constraints on axion properties.

## Files

### Main Analysis Scripts

- **Ligo_1dL.jl** - Single run LIGO limit computation with custom parameter override
  - Performs 1D likelihood profile analysis on binary black hole data
  - Output: MCMC samples saved to `../../output_mcmc/`

- **Ligo_1dL_cosma.jl** - Batch processing variant for HPC/COSMA environments
  - Supports parameter sweeps over axion mass ranges
  - Designed for parallel job submission
  - Output: Multiple MCMC chains saved to `../../output_mcmc/`

- **Ligo_limit.jl** - Post-processing script to derive exclusion limits
  - Computes limits from MCMC chains using KDE
  - Supports multiple statistical methods (KDE, GA)
  - Output: Compiled limits saved to `../../output_mcmc/Lim_*.dat`

### Supporting Files

- **stat_1dL.jl** - Core statistical functions for likelihood profiling
  - Implements MCMC sampler with Turing.jl
  - Handles MCMC diagnostic checking and convergence
  - Called by all main analysis scripts

## Usage

### Basic Usage

```bash
# Run single analysis
cd /path/to/scripts/ligo
julia Ligo_1dL.jl --dataname LIGO_test --ax_mass 1e-13 --Nmax 3

# Compute limits from existing chains
julia Ligo_limit.jl --dataname LIGO_test
```

### Command Line Arguments

**Ligo_1dL.jl**:
- `--dataname`: Input data file name (default: LIGO_test)
- `--ax_mass`: Axion mass in GeV (default: 1e-13)
- `--fa_min`: Min axion decay constant (default: 1e11)
- `--fa_max`: Max axion decay constant (default: 1e20)
- `--tau_max_override`: Evolution time in years (default: 4.5e7)
- `--Nmax`: Max principal quantum number (default: 3)
- `--Nsamples`: MCMC samples per walker (default: 1)
- `--non_rel`: Use non-relativistic approximation (default: false)
- `--use_kde`: Use KDE for limit computation (default: true)
- `--one_BH`: Analyze single BH population (default: false)
- `--cut_high_spin`: Apply spin cutoff at a=0.89 (default: false)

**Ligo_1dL_cosma.jl**:
- Supports parameter sweep: `--ax_mass_min`, `--ax_mass_max`, `--ax_mass_N`, `--mass_idx`
- Other parameters same as Ligo_1dL.jl

**Ligo_limit.jl**:
- `--dataname`: Matches MCMC output prefix
- `--tau_max`: Evolution time (must match MCMC runs)
- `--Nmax`: Quantum number range (must match MCMC)
- `--kde_test`: Enable diagnostic output (default: false)

## Input Data

LIGO data files should be in `../../BH_data/` directory with format:
```
M_1    M_2    chi_1    chi_2
[data rows]
```

Where:
- M_1, M_2: Component masses in solar masses
- chi_1, chi_2: Dimensionless spins (0-1)

## Output

All output files are saved to `../../output_mcmc/`:

- `LIGO_<name>_..._mcmc.dat` - MCMC samples from Turing.jl sampler
- `Lim_LIGO_<name>_....dat` - Derived exclusion limits (M, f_a)
- `KDE_*.dat` - Diagnostic KDE data (if --kde_test enabled)

## Performance Notes

- Single analysis typically requires 30-60 min for Nmax=3
- MCMC convergence checked via MCMCDiagnosticTools
- For large parameter sweeps, use HPC submission with Ligo_1dL_cosma.jl
- Output MCMC chains can be appended (see over_run setting in stat_1dL.jl)

## Path Handling

All paths are relative to `scripts/ligo/`:
- Input data: `../../BH_data/`
- Output: `../../output_mcmc/`
- Code access: stat_1dL.jl (same directory)

These scripts assume they are run from the `scripts/ligo/` directory or that the working directory is set correctly.

## References

- Primary likelihood in super_rad.jl (via stat_1dL.jl)
- MCMC diagnostics: MCMCDiagnosticTools.jl
- Turing.jl documentation: https://turinglang.org/
