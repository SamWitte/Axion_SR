#!/usr/bin/env python3
"""
Plot eigenvalue data from HDF5 files using matplotlib.

Reads computed eigenvalue data from HDF5 files and creates publication-quality
plots with log y-scale using matplotlib and seaborn styling.

Usage:
    python scripts/plot_eigenvalues.py --input data/eigenvalues_full.h5 --output plts/
    python scripts/plot_eigenvalues.py --input data/eigenvalues_quick.h5 --output plts/ --quick
"""

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from itertools import cycle

# Set up matplotlib style
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 7)
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Plot eigenvalue data from HDF5 files with matplotlib'
    )
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Path to input HDF5 file'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='plts/',
        help='Output directory for plots (default: plts/)'
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Quick mode (skip some plots for speed)'
    )
    parser.add_argument(
        '--format',
        type=str,
        default='png',
        choices=['png', 'pdf', 'svg'],
        help='Output image format (default: png)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='DPI for raster formats (default: 150)'
    )
    return parser.parse_args()


def load_data(h5_file):
    """Load all data from HDF5 file and organize by quantum level.

    Returns:
        Dict mapping (n,l,m) -> {spin -> (mu, alpha, imag_eigenvalue)}
    """
    data_by_level = {}

    with h5py.File(h5_file, 'r') as f:
        for group_key in f.keys():
            if group_key.startswith('level_'):
                # Parse group name: level_n_l_m
                parts = group_key.split('_')
                n, l, m = int(parts[1]), int(parts[2]), int(parts[3])
                quantum_level = (n, l, m)

                if quantum_level not in data_by_level:
                    data_by_level[quantum_level] = {}

                # Load each spin subgroup
                level_group = f[group_key]
                for spin_key in level_group.keys():
                    if spin_key.startswith('spin_'):
                        spin_val = float(spin_key.split('_')[1])
                        spin_group = level_group[spin_key]

                        mu = np.array(spin_group['mu'])
                        alpha = np.array(spin_group['alpha'])
                        imag_erg = np.array(spin_group['imag_eigenvalue'])

                        data_by_level[quantum_level][spin_val] = {
                            'mu': mu,
                            'alpha': alpha,
                            'imag_eigenvalue': imag_erg
                        }

    return data_by_level


def create_plots(data_by_level, output_dir, output_format='png', dpi=150):
    """Create plots for each quantum level.

    Args:
        data_by_level: Dictionary of data organized by quantum level
        output_dir: Directory to save plots
        output_format: Image format (png, pdf, svg)
        dpi: DPI for raster formats
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Color palette for different spins
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    for (n, l, m), spin_data in sorted(data_by_level.items()):
        print(f"\nCreating plot for (n,l,m) = ({n},{l},{m})...")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

        # Plot 1: vs mu
        for i, (spin, data) in enumerate(sorted(spin_data.items())):
            mu = data['mu']
            imag_erg = (data['imag_eigenvalue'])  # Use absolute value for log scale

            ax1.semilogy(
                mu, imag_erg,
                marker='o',
                label=f'a = {spin:.2f}',
                color=colors[i % len(colors)],
                alpha=0.8
            )

        ax1.set_xlabel(r'Boson mass $\mu$ (GeV)', fontsize=12)
        ax1.set_ylabel(r'$|\mathrm{Im}(\omega)|$ (log scale)', fontsize=12)
        ax1.set_title(f'Eigenvalue Im(ω) vs Boson Mass\n(n,l,m) = ({n},{l},{m})', fontsize=13)
        ax1.grid(True, which='both', alpha=0.3)
        ax1.legend(loc='best', framealpha=0.95)
        ax1.set_xscale('log')

        # Plot 2: vs alpha
        for i, (spin, data) in enumerate(sorted(spin_data.items())):
            alpha = data['alpha']
            imag_erg = (data['imag_eigenvalue'])

            ax2.semilogy(
                alpha, imag_erg,
                marker='s',
                label=f'a = {spin:.2f}',
                color=colors[i % len(colors)],
                alpha=0.8
            )

        ax2.set_xlabel(r'Dimensionless parameter $\alpha = \mu M G_N$', fontsize=12)
        ax2.set_ylabel(r'$|\mathrm{Im}(\omega)|$ (log scale)', fontsize=12)
        ax2.set_title(f'Eigenvalue Im(ω) vs α\n(n,l,m) = ({n},{l},{m})', fontsize=13)
        ax2.grid(True, which='both', alpha=0.3)
        ax2.legend(loc='best', framealpha=0.95)
        ax2.set_xscale('log')

        plt.tight_layout()

        # Save plot
        filename = f"eigenvalue_{n:d}{l:d}{m:d}.{output_format}"
        filepath = output_path / filename

        if output_format in ['png', 'jpg', 'jpeg']:
            plt.savefig(filepath, dpi=dpi, bbox_inches='tight')
        else:
            plt.savefig(filepath, bbox_inches='tight')

        print(f"  ✓ Saved: {filename}")
        plt.close(fig)


def print_summary(data_by_level):
    """Print summary statistics."""
    print("\n" + "="*80)
    print("EIGENVALUE DATA SUMMARY")
    print("="*80)
    print(f"{'Level':<12} {'Spin':<8} {'N pts':<8} {'Im(ω) range':<25} {'Mean Im(ω)':<15}")
    print("-"*80)

    for (n, l, m), spin_data in sorted(data_by_level.items()):
        for spin, data in sorted(spin_data.items()):
            imag_erg = data['imag_eigenvalue']
            n_points = len(imag_erg)
            im_min = np.min(imag_erg)
            im_max = np.max(imag_erg)
            im_mean = np.mean(imag_erg)

            level_str = f"({n},{l},{m})"
            range_str = f"[{im_min:.3e}, {im_max:.3e}]"
            print(f"{level_str:<12} {spin:<8.2f} {n_points:<8} {range_str:<25} {im_mean:<15.3e}")


def main():
    """Main function."""
    args = parse_arguments()

    h5_file = Path(args.input)
    if not h5_file.exists():
        print(f"Error: Input file not found: {h5_file}")
        return 1

    print("="*80)
    print("EIGENVALUE PLOTTING SCRIPT")
    print("="*80)
    print(f"Input file: {h5_file}")
    print(f"Output directory: {args.output}")
    print(f"Output format: {args.format}")

    print("\n[LOADING] Reading HDF5 data...")
    data_by_level = load_data(str(h5_file))
    print(f"  ✓ Loaded {len(data_by_level)} quantum levels")

    # Print summary
    print_summary(data_by_level)

    print(f"\n[PLOTTING] Creating {len(data_by_level)} plots...")
    create_plots(
        data_by_level,
        args.output,
        output_format=args.format,
        dpi=args.dpi
    )

    print("\n" + "="*80)
    print("PLOTTING COMPLETE")
    print("="*80)
    print(f"Plots saved to: {args.output}")

    return 0


if __name__ == '__main__':
    exit(main())
