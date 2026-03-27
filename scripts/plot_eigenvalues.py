#!/usr/bin/env python3
"""
Plot eigenvalue data from HDF5 files using matplotlib.

Reads computed eigenvalue data from HDF5 files and creates publication-quality
plots. States can be plotted individually (default) or grouped onto a single
set of axes with --combine.

Usage:
    # One plot per state (default):
    python scripts/plot_eigenvalues.py --input data/eigenvalues_full.h5 --output plts/

    # All states on one plot:
    python scripts/plot_eigenvalues.py --input data/eigenvalues_full.h5 --output plts/ --combine

    # Only selected states on one plot:
    python scripts/plot_eigenvalues.py --input data/eigenvalues_full.h5 --output plts/ \
        --combine --states 211 322 321
"""

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
HBAR_EV_S      = 6.582119e-16   # ħ in eV·s
AGE_UNIVERSE_S = 13.8e9 * 3.15e7  # 13.8 Gyr in seconds
# Superradiance rate (in eV) corresponding to one Hubble time
GAMMA_UNIVERSE = HBAR_EV_S / AGE_UNIVERSE_S         # ≈ 1.51e-33 eV  (1/τ_U)
GAMMA_YMIN     = GAMMA_UNIVERSE / 10.0              # 10× age of universe (y-axis floor)
GAMMA_1YR      = HBAR_EV_S / 3.15e7                 # 1/yr ≈ 2.09e-23 eV

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
plt.rcParams['text.usetex']       = True
plt.rcParams['font.family']       = 'serif'
plt.rcParams['font.serif']        = ['Computer Modern']
plt.rcParams['figure.figsize']    = (8, 6)
plt.rcParams['font.size']         = 12
plt.rcParams['axes.labelsize']    = 14
plt.rcParams['axes.titlesize']    = 16
plt.rcParams['xtick.labelsize']   = 12
plt.rcParams['ytick.labelsize']   = 12
plt.rcParams['legend.fontsize']   = 11
plt.rcParams['lines.linewidth']   = 1.0
plt.rcParams['axes.facecolor']    = 'white'
plt.rcParams['figure.facecolor']  = 'white'
plt.rcParams['axes.edgecolor']    = 'black'
plt.rcParams['axes.grid']         = False
plt.rcParams['xtick.direction']   = 'in'
plt.rcParams['ytick.direction']   = 'in'

# Publication-quality color cycle (colorblind-friendly)
STATE_COLORS = [
    '#1f77b4',  # blue
    '#d62728',  # red
    '#2ca02c',  # green
    '#9467bd',  # purple
    '#e377c2',  # pink
    '#8c564b',  # brown
    '#17becf',  # cyan
    '#ff7f0e',  # orange
    '#bcbd22',  # olive
    '#7f7f7f',  # grey
]
# Line styles cycle across spins within a state
SPIN_STYLES = ['-', '--', ':', '-.']


def _apply_tick_style(ax):
    """Major and minor ticks on all four sides, both pointing inward."""
    ax.minorticks_on()
    for which, length in [('major', 6), ('minor', 3)]:
        ax.tick_params(which=which, direction='in', length=length,
                       top=True, right=True, bottom=True, left=True)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Plot eigenvalue data from HDF5 files'
    )
    parser.add_argument('--input',   '-i', type=str, required=True,
                        help='Path to input HDF5 file')
    parser.add_argument('--output',  '-o', type=str, default='plts/',
                        help='Output directory for plots (default: plts/)')
    parser.add_argument('--format',        type=str, default='png',
                        choices=['png', 'pdf', 'svg'],
                        help='Output image format (default: png)')
    parser.add_argument('--dpi',           type=int, default=150,
                        help='DPI for raster formats (default: 150)')
    parser.add_argument('--combine',       action='store_true',
                        help='Plot all selected states on a single set of axes')
    parser.add_argument('--states', nargs='*', metavar='NLM',
                        help='States to include, e.g. --states 211 322 321. '
                             'Omit to include all states in the file.')
    return parser.parse_args()


def parse_state_filter(states_arg):
    """Convert ['211', '322'] -> {(2,1,1), (3,2,2)} or None for all."""
    if states_arg is None:
        return None
    result = set()
    for s in states_arg:
        s = s.strip()
        if len(s) == 3:
            result.add((int(s[0]), int(s[1]), int(s[2])))
        else:
            raise ValueError(f"Cannot parse state '{s}' — expected 3-digit string like '211'")
    return result


def load_data(h5_file, state_filter=None):
    """Load all data from HDF5 file, organised by quantum level.

    Returns:
        dict: (n,l,m) -> {spin_float -> {'mu', 'alpha', 'imag_eigenvalue'}}
    """
    data_by_level = {}
    with h5py.File(h5_file, 'r') as f:
        for group_key in f.keys():
            if not group_key.startswith('level_'):
                continue
            parts = group_key.split('_')
            n, l, m = int(parts[1]), int(parts[2]), int(parts[3])
            if state_filter is not None and (n, l, m) not in state_filter:
                continue
            if (n, l, m) not in data_by_level:
                data_by_level[(n, l, m)] = {}
            level_group = f[group_key]
            for spin_key in level_group.keys():
                if not spin_key.startswith('spin_'):
                    continue
                spin_val   = float(spin_key.split('_')[1])
                spin_group = level_group[spin_key]
                data_by_level[(n, l, m)][spin_val] = {
                    'mu':              np.array(spin_group['mu']),
                    'alpha':           np.array(spin_group['alpha']),
                    'imag_eigenvalue': np.array(spin_group['imag_eigenvalue']),
                }
    return data_by_level


ALPHA_XMIN = 0.06  # fixed x-axis lower cutoff


def _plot_state_on_ax(ax, spin_data, color, label_prefix='', alpha_lower=ALPHA_XMIN):
    """Draw all spin curves for one state onto *ax*."""
    for j, (spin, data) in enumerate(sorted(spin_data.items())):
        alpha    = data['alpha']
        imag_erg = data['imag_eigenvalue']
        mask     = (alpha >= alpha_lower) & (imag_erg >= GAMMA_YMIN)
        if not np.any(mask):
            continue
        spin_label = f'{label_prefix}$a={spin:.2f}$' if label_prefix else f'$a={spin:.2f}$'
        ax.semilogy(
            alpha[mask], imag_erg[mask],
            linestyle=SPIN_STYLES[j % len(SPIN_STYLES)],
            color=color,
            linewidth=1.0,
            label=spin_label,
        )


def _apply_axis_limits(ax):
    """Enforce fixed x/y limits and draw age-of-universe and 1-yr reference lines."""
    ax.set_xlim(left=ALPHA_XMIN)
    ax.set_ylim(bottom=GAMMA_YMIN)
    ax.axhline(GAMMA_UNIVERSE, color='k', linestyle='--', linewidth=0.8, alpha=0.4,
               label=r'$\tau_{\rm U}^{-1}$')
    ax.axhline(GAMMA_1YR, color='k', linestyle=':', linewidth=0.8, alpha=0.4,
               label=r'$1\,\mathrm{yr}^{-1}$')


def _savefig(fig, path, fmt, dpi):
    if fmt in ('png', 'jpg', 'jpeg'):
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
    else:
        fig.savefig(path, bbox_inches='tight')
    print(f"  Saved: {path.name}")


# ---------------------------------------------------------------------------
# Per-state plots (default mode)
# ---------------------------------------------------------------------------
def create_individual_plots(data_by_level, output_path, output_format, dpi):
    for i, ((n, l, m), spin_data) in enumerate(sorted(data_by_level.items())):
        print(f"  Plotting ({n},{l},{m})...")
        fig, ax = plt.subplots()

        _plot_state_on_ax(ax, spin_data, color=STATE_COLORS[i % len(STATE_COLORS)])

        ax.set_xlabel(r'$\alpha$')
        ax.set_ylabel(r'$\Gamma\ [\mathrm{eV}]$')
        ax.set_title(f'$({n},{l},{m})$')
        ax.set_xscale('log')
        _apply_axis_limits(ax)
        ax.legend(loc='best', framealpha=0.95)
        _apply_tick_style(ax)
        plt.tight_layout()

        _savefig(fig, output_path / f"eigenvalue_{n}{l}{m}.{output_format}",
                 output_format, dpi)
        plt.close(fig)


# ---------------------------------------------------------------------------
# Combined plot (--combine mode)
# ---------------------------------------------------------------------------
def create_combined_plot(data_by_level, output_path, output_format, dpi):
    fig, ax = plt.subplots()

    for i, ((n, l, m), spin_data) in enumerate(sorted(data_by_level.items())):
        color = STATE_COLORS[i % len(STATE_COLORS)]
        label_prefix = f'$({n},{l},{m}),\\ $'
        _plot_state_on_ax(ax, spin_data, color=color, label_prefix=label_prefix)

    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\Gamma\ [\mathrm{eV}]$')
    ax.set_xscale('log')
    _apply_axis_limits(ax)
    ax.legend(loc='best', framealpha=0.95, fontsize=9)
    _apply_tick_style(ax)
    plt.tight_layout()

    state_tag = '_'.join(f'{n}{l}{m}' for (n, l, m) in sorted(data_by_level))
    _savefig(fig, output_path / f"eigenvalue_{state_tag}.{output_format}",
             output_format, dpi)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def print_summary(data_by_level):
    print("\n" + "="*80)
    print("EIGENVALUE DATA SUMMARY")
    print("="*80)
    print(f"{'Level':<12} {'Spin':<8} {'N pts':<8} {'Im(ω) range':<25}")
    print("-"*80)
    for (n, l, m), spin_data in sorted(data_by_level.items()):
        for spin, data in sorted(spin_data.items()):
            erg = data['imag_eigenvalue']
            print(f"({n},{l},{m}){'':<7} {spin:<8.2f} {len(erg):<8} "
                  f"[{erg.min():.3e}, {erg.max():.3e}]")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_arguments()

    h5_file = Path(args.input)
    if not h5_file.exists():
        print(f"Error: input file not found: {h5_file}")
        return 1

    state_filter = parse_state_filter(args.states)

    print("="*80)
    print("EIGENVALUE PLOTTING SCRIPT")
    print("="*80)
    print(f"Input : {h5_file}")
    print(f"Output: {args.output}  format={args.format}")
    if state_filter:
        print(f"States: {sorted(state_filter)}")
    if args.combine:
        print("Mode  : combined (all selected states on one axes)")
    else:
        print("Mode  : individual (one plot per state)")

    print("\n[LOADING]")
    data_by_level = load_data(str(h5_file), state_filter)
    print(f"  Loaded {len(data_by_level)} quantum levels")

    print_summary(data_by_level)

    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    print("\n[PLOTTING]")
    if args.combine:
        create_combined_plot(data_by_level, output_path, args.format, args.dpi)
    else:
        create_individual_plots(data_by_level, output_path, args.format, args.dpi)

    print("\n" + "="*80)
    print(f"Done. Plots in: {args.output}")
    return 0


if __name__ == '__main__':
    exit(main())
