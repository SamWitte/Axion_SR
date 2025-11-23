# Imaginary Eigenvalue Computation and Visualization

This pair of scripts computes and visualizes the imaginary parts of eigenvalues for quasi-bound states of scalar fields around Kerr black holes.

## Overview

The imaginary eigenvalue component **Im(ω)** indicates the instability/decay timescale of the system:
- **Im(ω) > 0**: Superradiant instability (exponential growth)
- **Im(ω) ≈ 0**: Marginally stable (critical spin)
- **Im(ω) < 0**: Stable/damped mode

## Workflow

### Step 1: Compute Eigenvalues (Julia)

```bash
julia scripts/plot_imaginary_eigenvalues.jl
```

**What it does:**
- Computes Im(ω) as a function of α = μ * M * G_N for a 1 solar mass black hole
- Ranges: α from 0.05 to 1.0, spin a from 0.1 to 0.99
- Levels: (2,1,1), (2,1,-1), (3,1,1), (3,2,2)
- Uses adaptive resolution near the zero-crossing (sign flip)

**Output:**
- `data/imaginary_eigenvalues_211.h5` - HDF5 file for level (2,1,1)
- `data/imaginary_eigenvalues_21-1.h5` - HDF5 file for level (2,1,-1)
- `data/imaginary_eigenvalues_311.h5` - HDF5 file for level (3,1,1)
- `data/imaginary_eigenvalues_322.h5` - HDF5 file for level (3,2,2)

Each HDF5 file contains:
```
/spin_0.1/alpha        -> Array of α values
/spin_0.1/im_part      -> Array of Im(ω) values
/spin_0.1/{n,l,m}      -> Metadata (quantum numbers)
/spin_0.3/...
... (one group per spin)
```

**Computational Details:**
- Uses `find_im_part()` from `src/solve_sr_rates.jl`
- High-precision eigenvalue computation with BigFloat precision
- Automatically increases resolution near critical spin (a > 0.95)
- Detects and refines across zero crossings

**Runtime:**
Approximately 5-15 minutes per level depending on convergence behavior (longer for higher spins and higher l values)

### Step 2: Visualize Eigenvalues (Python)

```bash
python scripts/plot_imaginary_eigenvalues.py
```

**What it does:**
- Reads HDF5 data files from Step 1
- Creates publication-quality plots with log-log scale
- Generates individual plots per level and a 2×2 comparison

**Output:**
- `plots/imaginary_eigenvalues_211.png` - Individual level plots
- `plots/imaginary_eigenvalues_21-1.png`
- `plots/imaginary_eigenvalues_311.png`
- `plots/imaginary_eigenvalues_322.png`
- `plots/imaginary_eigenvalues_comparison.png` - All levels 2×2 grid

**Plot features:**
- Log-log scale (both axes)
- Separate curves for each spin (0.1 to 0.99)
- Color gradient using viridis colormap
- Grid lines for readability
- Only plots positive Im(ω) values (physical region)

## Physical Interpretation

### α Parameter Space

The dimensionless parameter α = μ * M * G_N characterizes the superradiance regime:

- **α ~ 0.03 to 0.3**: Weak coupling (long timescales)
- **α ~ 0.3 to 1.0**: Strong coupling (critical region)
- **α > 1.0**: Beyond quasi-bound state regime

For a 1 solar mass BH and axion boson mass μ:
- α = 0.1  corresponds to μ ~ 1.34 × 10⁻¹¹ GeV
- α = 0.5  corresponds to μ ~ 6.68 × 10⁻¹¹ GeV
- α = 1.0  corresponds to μ ~ 1.34 × 10⁻¹⁰ GeV

### Superradiance Instability

The critical spin **a_crit** where Im(ω) = 0 depends on (n,l,m):

For the (2,1,1) mode: a_crit is smallest (easiest to destabilize)
For the (3,2,2) mode: a_crit is larger (harder to destabilize)

Higher m modes are easier to destabilize for a given (n,l).

## Fine Resolution Near Zero Crossing

The Julia script implements **adaptive resolution**:

1. First, a coarse grid (20 points) from α_min to α_max
2. When a sign flip is detected in Im(ω):
   - Automatically inserts 6 intermediate points
   - This captures the critical spin accurately
3. Final data is sorted by α for smooth plotting

This approach balances accuracy near criticality with computational efficiency.

## Customization

To modify computation parameters, edit `scripts/plot_imaginary_eigenvalues.jl`:

```julia
M_BH = 1.0              # Change black hole mass (solar masses)
levels = [...]          # Add/remove quantum numbers
spins = [...]           # Modify spin range
alpha_min = 0.05        # Change α range
alpha_max = 1.0
```

To modify visualization, edit `scripts/plot_imaginary_eigenvalues.py`:

```python
DPI = 150               # Output resolution
FIGSIZE = (12, 8)       # Figure size
```

## Requirements

### Julia packages:
- HDF5.jl
- DelimitedFiles.jl
- Suppressor.jl

### Python packages:
- h5py
- numpy
- matplotlib
- pathlib (standard library)

## Troubleshooting

### Computation fails for high spins (a > 0.95)
- Some extreme spins may have convergence issues
- The script increases `Ntot` (grid points) automatically for a > 0.95
- Consider reducing `numsamples_perwalker` in find_im_part if memory is limited

### Plots show no data
- Check that `data/imaginary_eigenvalues_*.h5` files were created
- Verify that Python can read HDF5 files: `python -c "import h5py; print('h5py OK')"`

### Visualization looks strange
- Ensure both positive and negative spins are computed (currently: 0.1 to 0.99)
- Log scale requires all values to be positive (script filters automatically)

## References

- Superradiance: Brito, Cardoso, Pani (2015)
- Quasi-bound states: Detweiler (1980), Zouros & Eardley (1979)
- Implementation: Based on `find_im_part()` eigenvalue solver
