#!/usr/bin/env python3
"""Decent (not paper-final) spin + state evolution plots for N-levels runs.

Pure matplotlib (no LaTeX / pgfplots). One colour per Nmax on spin panels
(default Matplotlib cycle); one colour per mode on occupation panels with
legend labels ``|n l m>``.

Requirements
------------
Python 3 with ``numpy`` and ``matplotlib`` (``pandas`` optional, faster I/O).

Input layout under ``--output-root``
------------------------------------
::

    BH_<tag>/fa_<tag>/alpha_<tag>/
        Time_*Nmax_K.dat  Spin_*Nmax_K.dat  Modes_*Nmax_K.dat  States_*Nmax_K.dat
    and/or
        NL_output_Nmax_K.tar.gz   # unpacked on demand

``BH_tag`` is ``10``, ``1e4``, …; ``fa_tag`` is ``1e14``, …; ``alpha_tag`` is
``0.2``, ``0.4``, …

``--Nmax`` ranges use the campaign set ``{3,4,5,6,7,8,15,18}`` only, so
``6-15`` → ``6 7 8 15`` (no empty panels for 9…14).

Outputs (under ``--figure-dir``, default ``<output-root>/../figures/preview``)
-----------------------------------------------------------------------------
- ``spin_BH_*_fa_*.{png,pdf,jpg}`` — one panel per ``--alpha``, curves per Nmax
- ``states_BH_*_fa_*_alpha_*.{png,pdf,jpg}`` — one file per alpha; rows = Nmax

Examples
--------
::

    cd src/Nlevels_Runs

    # Spin grid + state stacks
    python plot_evolution_preview.py \\
      --output-root /path/to/output_or_archive \\
      --MassBH 10 --f_a 1e18 \\
      --alpha 0.2 0.4 0.6 0.8 \\
      --Nmax 5-8 \\
      --figure-dir ./figures/preview

    # States only
    python plot_evolution_preview.py \\
      --output-root ./output/run_nodrag \\
      --MassBH 1e4 --f_a 1e14 --alpha 0.6 \\
      --Nmax 6-15 --what states --fmt pdf

    python plot_evolution_preview.py -h   # all flags
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import tarfile
from typing import Iterable, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from matplotlib import colormaps as _mpl_cmaps

    def _get_cmap(name: str):
        return _mpl_cmaps[name]

except Exception:  # pragma: no cover
    from matplotlib.cm import get_cmap as _get_cmap  # type: ignore

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
# Canonical Nmax values produced by the N-levels campaign / archives.
NMAX_AVAILABLE = (3, 4, 5, 6, 7, 8, 15, 18)
N_POINTS = 800
PEAK_THRESHOLD = 5e-12
YMIN_STATES = PEAK_THRESHOLD
YMAX_STATES = 1e5

# Colourblind-friendly extras beyond tab20 for many modes
_EXTRA = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
    "#666666",
]


def _fast_loadtxt(path: str) -> np.ndarray:
    """Whitespace load; prefer pandas C engine when available."""
    try:
        import pandas as pd

        return pd.read_csv(path, sep=r"\s+", header=None, engine="c").values
    except ImportError:
        return np.loadtxt(path)


def bh_tag(mbh: float) -> str:
    mbh = float(mbh)
    if abs(mbh - 10.0) < 1e-6:
        return "10"
    return f"1e{int(round(math.log10(mbh)))}"


def fa_tag(fa: float) -> str:
    return f"1e{int(round(math.log10(float(fa))))}"


def alpha_tag(alpha: float) -> str:
    return "%g" % float(alpha)


def leaf_dir(output_root: str, mbh: float, fa: float, alpha: float) -> str:
    return os.path.join(
        output_root, f"BH_{bh_tag(mbh)}", f"fa_{fa_tag(fa)}", f"alpha_{alpha_tag(alpha)}"
    )


def parse_nmax_list(values: Sequence[str]) -> list[int]:
    """Accept '6-15' and/or discrete ints.

    Ranges are intersected with the canonical campaign set
    ``NMAX_AVAILABLE = (3,4,5,6,7,8,15,18)``, so ``6-15`` → ``[6,7,8,15]``
    (not every integer in between).
    """
    out: list[int] = []
    available = list(NMAX_AVAILABLE)
    for v in values:
        if "-" in v:
            a, b = v.split("-", 1)
            try:
                lo, hi = int(a), int(b)
            except ValueError as exc:
                raise SystemExit(f"Bad --Nmax range {v!r}") from exc
            if hi < lo:
                lo, hi = hi, lo
            chosen = [n for n in available if lo <= n <= hi]
            if not chosen:
                raise SystemExit(
                    f"--Nmax {v!r} matches none of the available values {available}"
                )
            out.extend(chosen)
        else:
            n = int(v)
            if n not in available:
                print(
                    f"  warning: Nmax={n} is not in the usual set {available}; "
                    "will still try to load it",
                    flush=True,
                )
            out.append(n)
    return sorted(dict.fromkeys(out))


def log_sample_indices(time: np.ndarray, n_points: int) -> np.ndarray:
    pos = np.where(time > 0)[0]
    if len(pos) == 0:
        return np.array([], dtype=int)
    if len(pos) <= n_points:
        return pos
    t0, t1 = float(time[pos[0]]), float(time[pos[-1]])
    if t1 <= t0:
        return pos
    targets = np.logspace(np.log10(t0), np.log10(t1), n_points)
    idx = np.searchsorted(time[pos], targets)
    idx = np.clip(idx, 0, len(pos) - 1)
    return pos[np.unique(idx)]


def _glob_one(directory: str, prefix: str, nmax: int) -> Optional[str]:
    hits = [
        p
        for p in glob.glob(os.path.join(directory, f"{prefix}_*Nmax_{nmax}.dat"))
        if not os.path.basename(p).startswith("._")
    ]
    return hits[0] if hits else None


def ensure_unpacked(leaf: str, nmax: int, prefixes: Iterable[str]) -> bool:
    """Ensure requested prefixes exist; unpack tar.gz if needed."""
    missing = [p for p in prefixes if _glob_one(leaf, p, nmax) is None]
    if not missing:
        return True
    tar = os.path.join(leaf, f"NL_output_Nmax_{nmax}.tar.gz")
    if not os.path.isfile(tar):
        print(f"  [skip] Nmax={nmax}: missing {missing} and no {os.path.basename(tar)}")
        return False
    print(f"  unpacking {os.path.basename(tar)} (need {missing}) ...", flush=True)
    # Prefer selective extract of needed prefixes when possible
    with tarfile.open(tar, "r:gz") as tf:
        members = []
        for m in tf:
            base = os.path.basename(m.name)
            if base.startswith("._") or not base.endswith(".dat"):
                continue
            if any(base.startswith(f"{p}_") and f"Nmax_{nmax}.dat" in base for p in prefixes):
                members.append(m)
        if not members:
            # fall back: full extract of .dat (can be huge for States)
            members = [
                m
                for m in tf.getmembers()
                if m.name.endswith(".dat") and not os.path.basename(m.name).startswith("._")
            ]
        tf.extractall(leaf, members=members)
    return all(_glob_one(leaf, p, nmax) is not None for p in prefixes)


def load_time_spin(leaf: str, nmax: int) -> Optional[tuple[np.ndarray, np.ndarray]]:
    if not ensure_unpacked(leaf, nmax, ("Time", "Spin")):
        return None
    tp, sp = _glob_one(leaf, "Time", nmax), _glob_one(leaf, "Spin", nmax)
    if not tp or not sp:
        return None
    return _fast_loadtxt(tp).ravel(), _fast_loadtxt(sp).ravel()


def load_states_panel(
    leaf: str, nmax: int, n_points: int = N_POINTS, peak_thr: float = PEAK_THRESHOLD
) -> Optional[tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]]:
    """Return (t_ds, modes, states_ds[n_active, n_t], labels) for active modes."""
    if not ensure_unpacked(leaf, nmax, ("Time", "Modes", "States")):
        return None
    tp = _glob_one(leaf, "Time", nmax)
    mp = _glob_one(leaf, "Modes", nmax)
    sp = _glob_one(leaf, "States", nmax)
    if not (tp and mp and sp):
        return None

    time = _fast_loadtxt(tp).ravel()
    modes = np.atleast_2d(_fast_loadtxt(mp))
    idx = log_sample_indices(time, n_points)
    if len(idx) == 0:
        return None
    t_ds = time[idx]

    # Stream States rows; keep only those above peak threshold on sampled cols
    active_rows: list[np.ndarray] = []
    active_labels: list[str] = []
    nmodes = modes.shape[0]
    with open(sp) as fh:
        for r, line in enumerate(fh):
            if r >= nmodes:
                break
            if not line.strip():
                continue
            toks = line.split()
            if len(toks) < int(idx[-1]) + 1:
                continue
            row = np.array([float(toks[i]) for i in idx], dtype=float)
            if np.nanmax(row) < peak_thr:
                continue
            n, l, m = int(modes[r, 0]), int(modes[r, 1]), int(modes[r, 2])
            active_rows.append(row)
            active_labels.append(rf"$|{n}\,{l}\,{m}\rangle$")

    if not active_rows:
        return t_ds, modes, np.zeros((0, len(idx))), []
    return t_ds, modes, np.vstack(active_rows), active_labels


def mode_colors(n: int) -> list:
    base = list(plt.rcParams["axes.prop_cycle"].by_key().get("color", []))
    cmap = _get_cmap("tab20")
    tab = [cmap(i / 20.0) for i in range(20)]
    palette = base + tab + _EXTRA
    return [palette[i % len(palette)] for i in range(n)]


def nmax_colors(nmax_list: Sequence[int]) -> dict[int, object]:
    """Default Matplotlib colour cycle (C0, C1, ...), one colour per Nmax."""
    cycle = list(plt.rcParams["axes.prop_cycle"].by_key().get("color", []))
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {n: cycle[i % len(cycle)] for i, n in enumerate(nmax_list)}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_spin_grid(
    output_root: str,
    mbh: float,
    fa: float,
    alphas: Sequence[float],
    nmax_list: Sequence[int],
    figure_path: str,
    ncols: int = 2,
) -> None:
    alphas = list(alphas)
    nmax_list = list(nmax_list)
    colors = nmax_colors(nmax_list)
    ncols = max(1, min(int(ncols), len(alphas)))
    nrows = int(math.ceil(len(alphas) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.2 * ncols, 3.4 * nrows),
        sharex=False,
        sharey=False,
        squeeze=False,
    )
    fig.suptitle(
        rf"Spin evolution — $M={bh_tag(mbh)}\,M_\odot$, $f_a={fa_tag(fa)}\,\mathrm{{eV}}$",
        fontsize=12,
        y=1.01,
    )

    for i, alpha in enumerate(alphas):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        leaf = leaf_dir(output_root, mbh, fa, alpha)
        ax.set_xscale("log")
        ax.set_ylim(0.0, 1.0)
        ax.set_title(rf"$\alpha={alpha_tag(alpha)}$", fontsize=11)
        ax.set_ylabel(r"$\tilde{a}$")
        if not os.path.isdir(leaf):
            ax.text(0.5, 0.5, "no leaf", transform=ax.transAxes, ha="center", color="gray")
            continue
        t_min, t_max = None, None
        for nmax in nmax_list:
            loaded = load_time_spin(leaf, nmax)
            if loaded is None:
                continue
            time, spin = loaded
            idx = log_sample_indices(time, N_POINTS)
            if len(idx) == 0:
                continue
            t, a = time[idx], spin[idx]
            ax.plot(t, a, color=colors[nmax], lw=1.6, label=rf"$N_{{\rm max}}={nmax}$")
            t_min = float(t.min()) if t_min is None else min(t_min, float(t.min()))
            t_max = float(t.max()) if t_max is None else max(t_max, float(t.max()))
        if t_min is not None and t_max is not None and t_max > t_min:
            ax.set_xlim(t_min, t_max)
        ax.legend(fontsize=8, loc="lower left", framealpha=0.75)
        ax.set_xlabel(r"$t$ [yr]")

    # Remove unused subplot slots entirely (avoids empty framed panel).
    for j in range(len(alphas), nrows * ncols):
        r, c = divmod(j, ncols)
        fig.delaxes(axes[r][c])

    os.makedirs(os.path.dirname(figure_path) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {figure_path}", flush=True)


def plot_states_stack(
    output_root: str,
    mbh: float,
    fa: float,
    alpha: float,
    nmax_list: Sequence[int],
    figure_path: str,
    peak_thr: float = PEAK_THRESHOLD,
) -> None:
    nmax_list = list(nmax_list)
    leaf = leaf_dir(output_root, mbh, fa, alpha)
    n = len(nmax_list)
    fig, axes = plt.subplots(
        n,
        1,
        figsize=(11, 2.5 * n),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    fig.suptitle(
        rf"States — $M={bh_tag(mbh)}\,M_\odot$, $f_a={fa_tag(fa)}\,\mathrm{{eV}}$, "
        rf"$\alpha={alpha_tag(alpha)}$",
        fontsize=12,
        y=1.01,
    )

    if not os.path.isdir(leaf):
        axes[0].text(0.5, 0.5, f"missing leaf\n{leaf}", transform=axes[0].transAxes, ha="center")
        for ax in axes[1:]:
            ax.axis("off")
    else:
        t_xlim = None
        for i, nmax in enumerate(nmax_list):
            ax = axes[i]
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim(YMIN_STATES, YMAX_STATES)
            ax.set_ylabel(rf"$N_{{\rm max}}={nmax}$" + "\n" + r"$M_c/M_i$")

            panel = load_states_panel(leaf, nmax, N_POINTS, peak_thr)
            if panel is None:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", color="gray")
                continue
            t_ds, _modes, states, labels = panel
            if states.shape[0] == 0:
                ax.text(
                    0.5,
                    0.5,
                    "no active modes",
                    transform=ax.transAxes,
                    ha="center",
                    color="gray",
                )
                continue
            cols = mode_colors(states.shape[0])
            for k in range(states.shape[0]):
                ax.plot(t_ds, np.clip(states[k], YMIN_STATES, None), color=cols[k], lw=1.0, label=labels[k])
            ncol = 2 if states.shape[0] <= 12 else 3
            if states.shape[0] > 24:
                ncol = 4
            ax.legend(
                fontsize=6,
                loc="upper right",
                ncol=ncol,
                framealpha=0.65,
                handlelength=1.1,
            )
            t_xlim = (float(t_ds.min()), float(t_ds.max())) if t_xlim is None else (
                min(t_xlim[0], float(t_ds.min())),
                max(t_xlim[1], float(t_ds.max())),
            )

        if t_xlim is not None and t_xlim[1] > t_xlim[0]:
            axes[-1].set_xlim(*t_xlim)
        axes[-1].set_xlabel(r"$t$ [yr]")

    os.makedirs(os.path.dirname(figure_path) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {figure_path}", flush=True)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--output-root",
        required=True,
        help="Root with BH_*/fa_*/alpha_*/ leaves (archive or campaign output)",
    )
    ap.add_argument("--MassBH", type=float, required=True)
    ap.add_argument("--f_a", type=float, required=True)
    ap.add_argument(
        "--alpha",
        type=float,
        nargs="+",
        required=True,
        help="One or more alpha values (each = one spin panel / one states figure)",
    )
    ap.add_argument(
        "--Nmax",
        nargs="+",
        default=["5-8"],
        help="Nmax list or range over {3,4,5,6,7,8,15,18}, e.g. 6-15 → 6 7 8 15 "
        "(default: 5-8)",
    )
    ap.add_argument(
        "--what",
        choices=("both", "spin", "states"),
        default="both",
        help="Which figures to write (default: both)",
    )
    ap.add_argument(
        "--figure-dir",
        default=None,
        help="Output directory (default: <output-root>/../figures/preview)",
    )
    ap.add_argument(
        "--peak-threshold",
        type=float,
        default=PEAK_THRESHOLD,
        help=f"Min peak occupation to draw a state (default {PEAK_THRESHOLD:g})",
    )
    ap.add_argument("--ncols", type=int, default=2, help="Columns in spin grid (default 2)")
    ap.add_argument(
        "--fmt",
        choices=("png", "pdf", "jpg"),
        default="png",
        help="Image format (default png)",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    nmax_list = parse_nmax_list(args.Nmax)
    alphas = list(args.alpha)
    root = os.path.abspath(args.output_root)
    fig_dir = args.figure_dir
    if fig_dir is None:
        fig_dir = os.path.join(os.path.dirname(root), "figures", "preview")
    fig_dir = os.path.abspath(fig_dir)
    os.makedirs(fig_dir, exist_ok=True)

    tag = f"BH_{bh_tag(args.MassBH)}_fa_{fa_tag(args.f_a)}"
    print(f"output-root: {root}")
    print(f"M={args.MassBH}  f_a={args.f_a}  alpha={alphas}  Nmax={nmax_list}")
    print(f"figure-dir: {fig_dir}")

    if args.what in ("both", "spin"):
        spin_path = os.path.join(fig_dir, f"spin_{tag}.{args.fmt}")
        plot_spin_grid(root, args.MassBH, args.f_a, alphas, nmax_list, spin_path, ncols=args.ncols)

    if args.what in ("both", "states"):
        for alpha in alphas:
            states_path = os.path.join(
                fig_dir, f"states_{tag}_alpha_{alpha_tag(alpha)}.{args.fmt}"
            )
            plot_states_stack(
                root,
                args.MassBH,
                args.f_a,
                alpha,
                nmax_list,
                states_path,
                peak_thr=args.peak_threshold,
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
