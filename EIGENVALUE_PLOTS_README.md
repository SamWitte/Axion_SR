# Eigenvalue Imaginary Component Plotting Suite

**Purpose**: Generate comprehensive visualizations of eigenvalue (energy level) imaginary components as functions of boson mass for various black hole spins.

**Status**: Ready for use
**Date Created**: November 2025

---

## Overview

This suite provides tools to study the superradiance instability through eigenvalue analysis. For a 1 solar mass black hole at various spins, the scripts compute and plot the imaginary parts of all eigenvalues across a range of boson masses.

### Physical Interpretation

- **Im(ω) > 0**: Unstable mode (superradiant growth)
- **Im(ω) ≈ 0**: Marginally stable mode
- **Im(ω) < 0**: Stable/decaying mode

The imaginary component indicates the growth/decay timescale of quantum axion modes around the black hole.

---

## Scripts Provided

### 1. `plot_eigenvalue_imaginary.jl` (Full Comprehensive Version)

**Purpose**: Generate complete eigenvalue analysis with all quantum levels and spin values.

**Features**:
- 6 spin parameters: a ∈ {0.0, 0.5, 0.7, 0.9, 0.95, 0.99}
- 50 boson mass points (log-spaced from 1e-22 to 1e-20)
- 7 quantum levels: (2,1,1), (3,1,1), (3,2,1), (3,2,2), (4,1,1), (4,3,1), (4,3,3)
- Uses Heun equation recursion relation approach
- Total computations: ~2,100

**Usage**:
```bash
julia plot_eigenvalue_imaginary.jl
```

**Output**:
- PDF plots in `plts/` directory
- One plot per quantum level: `Imag_erg_nlm.pdf`
- Summary statistics printed to console
- Computation time: ~2-4 hours (depending on system)

**Example Output Files**:
```
plts/Imag_erg_211.pdf    # (n=2, l=1, m=1)
plts/Imag_erg_312.pdf    # (n=3, l=1, m=2)
plts/Imag_erg_322.pdf    # (n=3, l=2, m=2)
... (one per quantum level)
```

---

### 2. `plot_eigenvalue_quick.jl` (Quick Test Version)

**Purpose**: Fast validation of the plotting pipeline.

**Features**:
- 3 spin parameters: a ∈ {0.5, 0.9, 0.95}
- 20 boson mass points (reduced from 50)
- 3 quantum levels (subset): (2,1,1), (3,2,1), (3,2,2)
- Uses Heun equation approach
- Total computations: 180 (vs 2,100 for full)

**Usage**:
```bash
julia plot_eigenvalue_quick.jl
```

**Output**:
- Quick test plots in `plts/` directory with `_quick` suffix
- Computation time: ~10-20 minutes
- Useful for verifying setup before running full analysis

---

## Technical Details

### Computational Method

Both scripts use the **Heun equation recursion relation** approach:

1. **Eigenvalue Computation**: `solve_radial()` function
   - Solves radial ODE with Heun equation
   - Uses `use_heunc=true` parameter
   - Computes complex eigenvalue ω = ωR + i·ωI

2. **Parameter Space**:
   - Black hole mass: M = 1 M☉ (fixed)
   - Boson mass: μ ∈ [1e-22, 1e-20] (log-spaced)
   - Black hole spin: a ∈ [0, 0.99] (various values)
   - Quantum numbers: (n, l, m) with n ≤ 4

3. **Data Organization**:
   - Results stored by (n,l,m,a) tuple
   - Sorted by boson mass for plotting
   - Complex eigenvalues extracted for imaginary component

### Plot Features

**For Each Quantum Level**:
- X-axis: log₁₀(μ) - logarithmic boson mass
- Y-axis: Im(ω) - imaginary eigenvalue component
- Multiple curves: One per spin value
- Different colors/markers for different spins
- Legend showing spin values

**Layout**:
- Left panel: Heun equation recursion relation results
- Right panel: Prepared for Chebyshev comparison (future enhancement)
- Comprehensive title with quantum numbers

---

## Customization Guide

### Modify Quantum Levels

Edit the `quantum_levels` list in either script:

```julia
quantum_levels = [
    (2, 1, 1),  # Add or remove levels as needed
    (3, 2, 2),
    (4, 3, 3),
]
```

### Adjust Boson Mass Range

Modify the `mu_values` computation:

```julia
mu_min = 1e-23  # Change lower bound
mu_max = 1e-19  # Change upper bound
n_mu_points = 100  # Change number of points
mu_values = 10 .^ (range(log10(mu_min), log10(mu_max), n_mu_points))
```

### Change Spin Parameters

Modify the `spins` array:

```julia
spins = [0.0, 0.3, 0.6, 0.9, 0.99]  # Custom spin list
```

### Adjust Numerical Precision

Fine-tune ODE solver parameters in `solve_radial()` call:

```julia
_, _, erg = solve_radial(mu, M_BH, a, n, l, m;
                         rpts=2000,         # Grid resolution
                         rmaxT=100,         # Max radius
                         use_heunc=true,    # Method
                         return_erg=true,   # Return eigenvalue
                         Ntot_safe=10000,   # Safety limit
                         debug=false)       # Verbosity
```

---

## Output Interpretation Guide

### Reading the Plots

1. **Boson Mass Dependence**:
   - X-axis shows where eigenvalues become complex
   - Transition from real to complex indicates instability onset

2. **Spin Dependence**:
   - Different colors show how spin affects Im(ω)
   - Higher spins typically show larger imaginary components
   - Indicates faster growth/instability rates

3. **Quantum Level Dependence**:
   - Each plot shows one (n,l,m) combination
   - Compare across levels to understand hierarchy

### Statistics Output

Console output includes summary table:

```
Level (n,l,m) | Spin a  | Im(ω) range            | Mean Im(ω)
(2,1,1)       | 0.50    | [1.23e-30, 4.56e-25]  | 2.34e-25
(2,1,1)       | 0.90    | [3.45e-28, 1.23e-22]  | 5.67e-23
...
```

---

## Advanced Features

### Future Enhancements

The scripts are designed to accommodate:

1. **Chebyshev Decomposition Comparison**:
   - Right panel prepared for second method
   - Would show if different approaches agree
   - Add via `use_heunc=false` and appropriate Chebyshev solver

2. **More Quantum Levels**:
   - Can extend to n ≤ 6 (Nmax=6)
   - Requires more computation time
   - 52+ rates available per Nmax

3. **Extended Spin Range**:
   - Can go to a → 1 (near-extremal)
   - Requires more careful numerics
   - Physical interest for near-extremal black holes

4. **Statistical Analysis**:
   - Compute growth rates directly
   - Analyze bifurcation points
   - Study level crossing phenomena

---

## Practical Usage Guide

### Quick Start (Testing)

1. Verify code setup:
   ```bash
   julia plot_eigenvalue_quick.jl
   ```

2. Check output in `plts/` directory for `_quick` files

3. Examine statistics and plots for correctness

### Full Analysis

1. Run comprehensive version when ready:
   ```bash
   julia plot_eigenvalue_imaginary.jl
   ```

2. Monitor progress (printed every 10% of computations)

3. Collect all `Imag_erg_*.pdf` files from `plts/`

4. Use plots for:
   - Publication figures
   - Presentation materials
   - Parameter space exploration
   - Physics interpretation

### Batch Processing

To run multiple parameter sets:

```bash
# Create wrapper script
for mass in 0.5 1.0 2.0 5.0 10.0
    # Modify M_BH in script
    sed -i "s/M_BH = .*/M_BH = $mass/" plot_eigenvalue_imaginary.jl
    julia plot_eigenvalue_imaginary.jl
    # Save outputs with mass label
    mv plts/ "plts_M${mass}/"
done
```

---

## Requirements

### Julia Packages
- `Plots.jl` - Plotting
- `LaTeXStrings.jl` - LaTeX formatting
- `Statistics.jl` - Basic statistics
- `Printf.jl` - Formatting

### Internal Modules
- `solve_sr_rates.jl` - Eigenvalue computation
- `heunc.jl` - Heun equation solver
- `Constants.jl` - Physical constants

### Computational
- ~2-4 hours for full version
- ~10-20 minutes for quick version
- ~1-2 GB RAM
- Multi-core system recommended

---

## Troubleshooting

### Common Issues

**Issue**: "Package X not found"
```
Solution: Standard packages should be available
If missing, might indicate Julia installation issue
```

**Issue**: Computation takes too long
```
Solution: Use plot_eigenvalue_quick.jl instead
Or reduce n_mu_points or number of quantum levels
```

**Issue**: NaN or Inf in results
```
Solution: May occur at extreme parameter values
Try adjusting rpts or Ntot_safe parameters
Or reduce spin value if a > 0.99
```

**Issue**: Plot won't display/save
```
Solution: Ensure plts/ directory exists and is writable
Check PDF viewer can handle generated files
Verify Plots.jl is properly configured
```

---

## File Structure

```
Axion_SR/
├── plot_eigenvalue_imaginary.jl    # Full version (2-4 hours)
├── plot_eigenvalue_quick.jl        # Quick test (10-20 min)
├── EIGENVALUE_PLOTS_README.md      # This file
└── plts/
    ├── Imag_erg_211.pdf            # Output files
    ├── Imag_erg_312.pdf
    └── ...
```

---

## References

### Physical Background
- Superradiance instability for rotating black holes
- Scalar field perturbations and quasi-normal modes
- Heun equation formalism for radial equation

### Computational Methods
- `solve_radial()`: Radial ODE solver with eigenvalue extraction
- `heunc.jl`: Heun equation recursion relation implementation
- `Plots.jl`: PDF output generation

### Related Scripts
- `solve_sr_rates.jl`: Scattering rate calculations
- `super_rad.jl`: Superradiance timescale estimates
- `grid_sr_growth.jl`: Evolution simulations

---

## Author Notes

**Creation Date**: November 2025
**Purpose**: Enable efficient visualization and analysis of eigenvalue behavior across parameter space
**Maintenance**: Update quantum_levels or parameter ranges as needed for specific studies

**Key Design Decisions**:
- Separated into quick and full versions for flexibility
- Heun equation method provides reliable eigenvalues
- Log-spaced boson mass for better resolution in unstable region
- Comprehensive output for publication use

---

## Citation

If using these scripts for research, please cite:
- Original black hole superradiance work
- Axion dark matter research
- This visualization suite with version and date

---

**Status**: Ready for immediate use
**Last Updated**: November 2025
**Tested**: ✅ Yes (quick version verified)
