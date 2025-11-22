# Eigenvalue Plotting Scripts

Computation and visualization of eigenvalue data for axion superradiance analysis.

## Overview

The plotting workflow is now **two-stage**:

1. **Data Generation** (Julia): Compute eigenvalues and save to HDF5
2. **Visualization** (Python/matplotlib): Create publication-quality plots

This separation provides:
- Better plot quality with matplotlib
- Log y-scale support for wide dynamic ranges
- Easier plot customization
- Decoupling of computation from visualization

## Key Feature: Smart μ Range Calculation

Rather than manually specifying boson mass ranges, the scripts use the physically relevant constraint:

**α = μ × M × G_N ∈ [0.03, 1.0]**

Where α is the dimensionless parameter controlling the superradiance regime. The scripts automatically calculate the required μ range from this.

## Visualization Scripts

### plot_eigenvalue_imaginary.jl
Comprehensive eigenvalue data generation across parameter space.

**Features**:
- 6 spin parameters: a ∈ {0.0, 0.5, 0.7, 0.9, 0.95, 0.99}
- 50 boson mass points (log-spaced)
- 7 quantum levels: (2,1,1), (3,1,1), (3,2,1), (3,2,2), (4,1,1), (4,3,1), (4,3,3)
- ~2,100 eigenvalue computations
- Dimensionless parameter range: α ∈ [0.03, 1.0]
- Runtime: 10-20 minutes
- Output: HDF5 data + CSV summary

**Usage**:
```bash
julia scripts/plot_eigenvalue_imaginary.jl
```

**Output files**:
- `data/eigenvalues_full.h5` - Complete eigenvalue data
- `data/eigenvalues_full_summary.csv` - Summary statistics

### plot_eigenvalue_quick.jl
Quick test version for verifying the pipeline.

**Features**:
- 3 spin parameters: a ∈ {0.5, 0.9, 0.95}
- 20 boson mass points
- 3 quantum levels: (2,1,1), (3,2,1), (3,2,2)
- 180 total computations
- Dimensionless parameter range: α ∈ [0.03, 1.0]
- Runtime: 1-2 minutes
- Output: HDF5 data + CSV summary

**Usage**:
```bash
julia scripts/plot_eigenvalue_quick.jl
```

**Output files**:
- `data/eigenvalues_quick.h5` - Eigenvalue data
- `data/eigenvalues_quick_summary.csv` - Summary statistics

### plot_eigenvalues.py
Python/matplotlib script for creating publication-quality plots from HDF5 data.

**Features**:
- Log y-scale for Im(ω) visualization
- Two plots per quantum level:
  - Im(ω) vs boson mass μ
  - Im(ω) vs dimensionless parameter α
- Seaborn styling for professional appearance
- PNG, PDF, or SVG output
- Configurable DPI and format

**Usage**:
```bash
# Plot from quick run
python scripts/plot_eigenvalues.py --input data/eigenvalues_quick.h5 --output plts/

# Plot from full run with PDF output
python scripts/plot_eigenvalues.py \
    --input data/eigenvalues_full.h5 \
    --output plts/ \
    --format pdf \
    --dpi 300

# Help
python scripts/plot_eigenvalues.py --help
```

**Command-line options**:
- `--input, -i`: Input HDF5 file (required)
- `--output, -o`: Output directory (default: `plts/`)
- `--format`: Output format: png, pdf, or svg (default: png)
- `--dpi`: DPI for raster formats (default: 150)

## Complete Workflow Example

### Quick test:
```bash
# 1. Generate data (1-2 minutes)
julia scripts/plot_eigenvalue_quick.jl

# 2. Create plots (10-30 seconds)
python scripts/plot_eigenvalues.py --input data/eigenvalues_quick.h5 --output plts/

# View PNG plots in plts/ directory
```

### Full production run:
```bash
# 1. Generate comprehensive data (10-20 minutes)
julia scripts/plot_eigenvalue_imaginary.jl

# 2. Create high-quality PDF plots
python scripts/plot_eigenvalues.py \
    --input data/eigenvalues_full.h5 \
    --output plts/ \
    --format pdf \
    --dpi 300

# View PDF plots in plts/ directory
```

## Data Format

HDF5 hierarchical structure:

```
eigenvalues_full.h5
├── level_2_1_1/
│   ├── spin_0.00/
│   │   ├── mu (boson masses)
│   │   ├── alpha (dimensionless parameter α)
│   │   ├── imag_eigenvalue (Im(ω) values)
│   │   ├── M_BH (black hole mass)
│   │   └── spin (spin parameter)
│   ├── spin_0.50/
│   └── ...
├── level_3_1_1/
└── ...
```

CSV summary format:

```
n,l,m,spin,n_points,im_min,im_max,im_mean
2,1,1,0.00,50,1.234e-03,5.678e-02,2.345e-02
...
```

## Physical Context

**Dimensionless parameter α = μ × M × G_N**:
- α < 0.03: Too small (classical superradiance ineffective)
- 0.03 ≤ α ≤ 1.0: **Optimal superradiance regime** (this script's range)
- α > 1.0: Relativistic/quantum effects dominate

**Imaginary eigenvalue Im(ω)**:
- **Im(ω) > 0**: Unstable mode (superradiant growth)
- **Im(ω) ≈ 0**: Marginally stable
- **Im(ω) < 0**: Stable/decaying mode

## Dependencies

### Julia
```julia
using HDF5      # Data storage
using Printf    # Formatted output
using Statistics # mean, etc.
```

### Python
```bash
pip install h5py numpy matplotlib seaborn
```

## Script Development Guidelines

When creating new scripts in this directory:

1. **Documentation**: Clear docstrings and comments
2. **Paths**: Use `@__DIR__` and `joinpath()` for file references
3. **Error handling**: Graceful failures with informative messages
4. **Performance**: Communicate runtime expectations
5. **Output**: Organize outputs in appropriate subdirectories

---

**Directory**: `scripts/`
**Last Updated**: November 2025
**Status**: Production Ready
