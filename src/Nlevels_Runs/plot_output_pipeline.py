#!/usr/bin/env python3
"""
Unpack tarballs -> plot_runs (all Nmax panels) -> repack for complete leaves only.

When running under SLURM, unpack/prune/plot use node-local disk ($SLURM_TMPDIR
or /tmp); only tarballs, .pruned markers, and figures are written back to lustre.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time as _time

_HERE = os.path.dirname(os.path.abspath(__file__))
# Default to this directory so collaborators can run without lustre paths.
# Override with NL_BASE_DIR or --nl-base.
NL_BASE = os.environ.get("NL_BASE_DIR", _HERE)
RUN_TAG = os.environ.get("RUN_TAG", "").strip()


def _output_root(nl_base=None, run_tag=None):
    base = nl_base if nl_base is not None else NL_BASE
    tag = run_tag if run_tag is not None else RUN_TAG
    if tag:
        return os.path.join(base, "output", tag)
    return os.path.join(base, "output")


def _figures_root(nl_base=None, run_tag=None):
    base = nl_base if nl_base is not None else NL_BASE
    tag = run_tag if run_tag is not None else RUN_TAG
    if tag:
        return os.path.join(base, "figures", tag)
    return os.path.join(base, "figures")


OUTPUT_ROOT = _output_root()
FIGURES_ROOT = _figures_root()
COMPRESS_SH = os.path.join(NL_BASE, "compress_nmax_output.sh")


def configure_paths(nl_base=None, run_tag=None):
    """Set module-level paths (call after parsing CLI / env)."""
    global NL_BASE, RUN_TAG, OUTPUT_ROOT, FIGURES_ROOT, COMPRESS_SH
    if nl_base is not None:
        NL_BASE = nl_base
    if run_tag is not None:
        RUN_TAG = run_tag.strip()
    OUTPUT_ROOT = _output_root()
    FIGURES_ROOT = _figures_root()
    COMPRESS_SH = os.path.join(NL_BASE, "compress_nmax_output.sh")

NMAX_LIST = [3, 4, 5, 6, 7, 8, 15, 18]

from prune_states import prune_states_modes


def _log(msg):
    print(msg, flush=True)


def _init_plot_runs():
    import plot_runs

    plot_runs.OUTPUT_ROOT = OUTPUT_ROOT
    plot_runs.FIGURES_ROOT = FIGURES_ROOT
    plot_runs._HERE = NL_BASE
    return plot_runs


def _scratch_root():
    """Node-local scratch for unpack/prune/plot (never lustre)."""
    for key in ("SLURM_TMPDIR", "TMPDIR"):
        val = os.environ.get(key)
        if val and not val.startswith("/lustre"):
            return val
    job = os.environ.get("SLURM_JOB_ID")
    if job:
        task = os.environ.get("SLURM_ARRAY_TASK_ID", "0")
        return os.path.join("/tmp", f"nl_scratch_{job}_{task}")
    return None


def use_local_work():
    """Use node scratch only when output lives on lustre (cluster jobs)."""
    if not OUTPUT_ROOT.startswith("/lustre"):
        return False
    return _scratch_root() is not None


def work_leaf_path(lustre_leaf):
    tmp = _scratch_root()
    rel = os.path.relpath(lustre_leaf, OUTPUT_ROOT)
    safe = rel.replace(os.sep, "__")
    return os.path.join(tmp, f"leaf_{safe}")


def ensure_work_leaf(lustre_leaf):
    work = work_leaf_path(lustre_leaf)
    os.makedirs(work, exist_ok=True)
    return work


def remove_work_leaf(lustre_leaf):
    work = work_leaf_path(lustre_leaf)
    if os.path.isdir(work):
        shutil.rmtree(work, ignore_errors=True)


def tarball_path(directory, nmax):
    return os.path.join(directory, f"NL_output_Nmax_{nmax}.tar.gz")


def sync_tarball_from_lustre(lustre_leaf, work_leaf, nmax):
    src = tarball_path(lustre_leaf, nmax)
    if not os.path.isfile(src):
        return False
    dst = tarball_path(work_leaf, nmax)
    if os.path.isfile(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
        return True
    _log(f"    copy tarball Nmax={nmax} lustre -> node scratch")
    shutil.copy2(src, dst)
    return True


def sync_tarball_to_lustre(work_leaf, lustre_leaf, nmax):
    src = tarball_path(work_leaf, nmax)
    dst = tarball_path(lustre_leaf, nmax)
    _log(f"    copy tarball Nmax={nmax} node scratch -> lustre")
    shutil.copy2(src, dst)


def has_tarball(directory, nmax, min_bytes=1024):
    path = tarball_path(directory, nmax)
    return os.path.isfile(path) and os.path.getsize(path) >= min_bytes


def is_valid_tarball(directory, nmax):
    """False if missing, tiny, or tar/gzip cannot read the archive."""
    if not has_tarball(directory, nmax):
        return False
    path = tarball_path(directory, nmax)
    if subprocess.run(["gzip", "-t", path], capture_output=True).returncode != 0:
        return False
    result = subprocess.run(
        ["tar", "-tzf", path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def has_loose_nmax(directory, nmax):
    return any(
        f.endswith(f"Nmax_{nmax}.dat")
        for f in os.listdir(directory)
        if f.startswith(("Spin_", "Time_", "States_", "Modes_", "MassBH_"))
    )


def is_complete_leaf(directory):
    return all(os.path.isfile(tarball_path(directory, n)) for n in NMAX_LIST)


def figure_path(directory):
    try:
        rel = os.path.relpath(
            os.path.realpath(directory), os.path.realpath(OUTPUT_ROOT)
        )
    except ValueError:
        rel = os.path.basename(directory)
    parts = rel.split(os.sep)
    bh, fa, alpha = (parts + ["unknown", "unknown", "plot"])[:3]
    return os.path.join(FIGURES_ROOT, bh, fa, alpha + ".jpg")


def find_complete_leaves(root, bh_mass=None):
    leaves = []
    for bh in sorted(os.listdir(root)):
        if bh_mass is not None and bh != bh_mass:
            continue
        bh_path = os.path.join(root, bh)
        if not os.path.isdir(bh_path):
            continue
        for fa in sorted(os.listdir(bh_path)):
            fa_path = os.path.join(bh_path, fa)
            if not os.path.isdir(fa_path):
                continue
            for alpha in sorted(os.listdir(fa_path)):
                leaf = os.path.join(fa_path, alpha)
                if os.path.isdir(leaf) and is_complete_leaf(leaf):
                    leaves.append(leaf)
    return leaves


def find_pruneable_leaves(root, bh_mass=None):
    """Leaves with at least one present tarball that is not yet pruned."""
    leaves = []
    for bh in sorted(os.listdir(root)):
        if bh_mass is not None and bh != bh_mass:
            continue
        bh_path = os.path.join(root, bh)
        if not os.path.isdir(bh_path):
            continue
        for fa in sorted(os.listdir(bh_path)):
            fa_path = os.path.join(bh_path, fa)
            if not os.path.isdir(fa_path):
                continue
            for alpha in sorted(os.listdir(fa_path)):
                leaf = os.path.join(fa_path, alpha)
                if not os.path.isdir(leaf):
                    continue
                if any(
                    has_tarball(leaf, n) and not is_nmax_pruned(leaf, n)
                    for n in NMAX_LIST
                ):
                    leaves.append(leaf)
    return leaves


def find_leaves_missing_figures(root, bh_mass=None, include_partial=False):
    """Leaves without a figure yet; partial leaves need at least one tarball."""
    leaves = []
    for bh in sorted(os.listdir(root)):
        if bh_mass is not None and bh != bh_mass:
            continue
        bh_path = os.path.join(root, bh)
        if not os.path.isdir(bh_path):
            continue
        for fa in sorted(os.listdir(bh_path)):
            fa_path = os.path.join(bh_path, fa)
            if not os.path.isdir(fa_path):
                continue
            for alpha in sorted(os.listdir(fa_path)):
                leaf = os.path.join(fa_path, alpha)
                if not os.path.isdir(leaf):
                    continue
                if os.path.isfile(figure_path(leaf)):
                    continue
                if include_partial:
                    if any(has_tarball(leaf, n) for n in NMAX_LIST):
                        leaves.append(leaf)
                elif is_complete_leaf(leaf):
                    leaves.append(leaf)
    return leaves


def nmax_pruned_marker(directory, nmax):
    return os.path.join(directory, f"NL_output_Nmax_{nmax}.pruned")


def is_nmax_pruned(directory, nmax):
    """True only if a prune marker exists and is at least as new as the tarball.

    Stale markers (e.g. Jul-11 prune + Jul-27 Nmax=18 rerun) are treated as
    not pruned so re-runs get re-pruned.
    """
    marker = nmax_pruned_marker(directory, nmax)
    if not os.path.isfile(marker):
        return False
    tar = tarball_path(directory, nmax)
    if not os.path.isfile(tar):
        return True
    # Allow a few seconds of clock skew between marker touch and tar write.
    return os.path.getmtime(marker) + 5.0 >= os.path.getmtime(tar)


def clear_nmax_prune_marker(directory, nmax):
    marker = nmax_pruned_marker(directory, nmax)
    if os.path.isfile(marker):
        os.remove(marker)


def mark_nmax_pruned(directory, nmax):
    open(nmax_pruned_marker(directory, nmax), "w").close()


def _tar_members(directory, nmax, prefixes):
    archive = tarball_path(directory, nmax)
    result = subprocess.run(
        ["tar", "-tzf", archive], capture_output=True, text=True, check=True
    )
    members = []
    for line in result.stdout.strip().splitlines():
        base = os.path.basename(line)
        for prefix in prefixes:
            if base.startswith(f"{prefix}_") and base.endswith(f"Nmax_{nmax}.dat"):
                members.append(line)
                break
    return members


def unpack_nmax(directory, nmax, prefixes=None):
    if prefixes is None:
        prefixes = ("Time", "Spin", "States", "Modes", "MassBH")
    archive = tarball_path(directory, nmax)
    if not os.path.isfile(archive):
        return False
    members = _tar_members(directory, nmax, prefixes)
    if not members:
        return False
    missing = [
        m
        for m in members
        if not os.path.isfile(os.path.join(directory, os.path.basename(m)))
    ]
    if not missing:
        return True
    _log(
        f"    unpack Nmax={nmax} ({len(missing)} file(s)) from "
        f"{os.path.basename(archive)}"
    )
    subprocess.run(
        ["tar", "-xzf", archive, "-C", directory] + missing,
        check=True,
    )
    return True


def cleanup_loose_nmax(directory, nmax):
    """Remove unpacked files for one Nmax (e.g. after a failed prune/repack)."""
    removed = 0
    token = f"Nmax_{nmax}.dat"
    for name in os.listdir(directory):
        if token not in name:
            continue
        path = os.path.join(directory, name)
        if os.path.isfile(path):
            os.remove(path)
            removed += 1
    if removed:
        _log(f"    cleaned {removed} loose file(s) for Nmax={nmax}")


def repack_nmax(directory, nmax, lustre_leaf=None):
    if not has_loose_nmax(directory, nmax):
        return
    _log(f"    repack Nmax={nmax}")
    env = os.environ.copy()
    env["FORCE"] = "1"
    subprocess.run(
        ["bash", COMPRESS_SH, directory, str(nmax)],
        check=True,
        env=env,
    )
    marker_leaf = lustre_leaf if lustre_leaf is not None else directory
    if lustre_leaf is not None:
        sync_tarball_to_lustre(directory, lustre_leaf, nmax)
    mark_nmax_pruned(marker_leaf, nmax)


def unpack_leaf(directory):
    for nmax in NMAX_LIST:
        if not unpack_nmax(directory, nmax):
            raise RuntimeError(f"missing tarball for Nmax={nmax}")


def repack_leaf(directory):
    for nmax in NMAX_LIST:
        repack_nmax(directory, nmax)


def is_fully_pruned(directory):
    """True if every *present* Nmax tarball has a fresh prune marker."""
    present = [n for n in NMAX_LIST if has_tarball(directory, n)]
    if not present:
        return False
    return all(is_nmax_pruned(directory, n) for n in present)


def find_plottable_leaves(root, bh_mass=None):
    return [
        leaf
        for leaf in find_complete_leaves(root, bh_mass=bh_mass)
        if is_fully_pruned(leaf) and not os.path.isfile(figure_path(leaf))
    ]


def prepare_nmax_unpack_only(directory, nmax):
    if not unpack_nmax(directory, nmax):
        raise RuntimeError(f"missing tarball for Nmax={nmax}")


def process_leaf_prune_only(lustre_leaf, force=False):
    rel = os.path.relpath(lustre_leaf, OUTPUT_ROOT)
    local = use_local_work()
    work = ensure_work_leaf(lustre_leaf) if local else lustre_leaf
    if local:
        _log(f"  prune-only ({_scratch_root()}): {rel}")
    else:
        _log(f"  prune-only: {rel}")

    pruned_any = False
    try:
        for nmax in NMAX_LIST:
            if not has_tarball(lustre_leaf, nmax):
                _log(f"    skip Nmax={nmax} (no tarball)")
                continue
            if force:
                clear_nmax_prune_marker(lustre_leaf, nmax)
            if is_nmax_pruned(lustre_leaf, nmax):
                _log(f"    skip Nmax={nmax} (already pruned, marker fresh)")
                continue
            if not is_valid_tarball(lustre_leaf, nmax):
                _log(f"    skip Nmax={nmax} (corrupt tarball — re-rsync from cluster)")
                continue
            # Drop stale marker so repack can rewrite it after prune.
            clear_nmax_prune_marker(lustre_leaf, nmax)
            try:
                if local:
                    if not sync_tarball_from_lustre(lustre_leaf, work, nmax):
                        raise RuntimeError(f"missing tarball for Nmax={nmax}")
                if not unpack_nmax(work, nmax, ("States", "Modes")):
                    raise RuntimeError(f"missing tarball for Nmax={nmax}")
                # Fresh unpack: drop any leftover States_*.pruned from prior runs.
                for name in os.listdir(work):
                    if name.endswith(f"Nmax_{nmax}.dat.pruned"):
                        os.remove(os.path.join(work, name))
                prune_states_modes(work, nmax)
                if not unpack_nmax(work, nmax, ("Time", "Spin", "MassBH")):
                    raise RuntimeError(f"missing small files for Nmax={nmax}")
                repack_nmax(
                    work,
                    nmax,
                    lustre_leaf=lustre_leaf if local else None,
                )
                cleanup_loose_nmax(work, nmax)
                pruned_any = True
            except Exception:
                cleanup_loose_nmax(work, nmax)
                raise
    finally:
        if local:
            remove_work_leaf(lustre_leaf)
    return "pruned" if pruned_any else "skipped"


def process_leaf(directory, skip_existing=True, prune_only=False, require_pruned=True, force_reprune=False):
    rel = os.path.relpath(directory, OUTPUT_ROOT)
    out = figure_path(directory)
    if skip_existing and os.path.isfile(out):
        _log(f"  SKIP (figure exists): {rel}")
        return "skipped"

    if prune_only:
        return process_leaf_prune_only(directory, force=force_reprune)

    if require_pruned and not is_fully_pruned(directory):
        _log(f"  SKIP (not fully pruned): {rel}")
        return "not_ready"

    local = use_local_work()
    work = ensure_work_leaf(directory) if local else directory
    if local:
        _log(f"  plot (node scratch, read-only tarballs on lustre): {rel}")
    else:
        _log(f"  plot (unpack/repack per Nmax): {rel}")

    def prepare_nmax(nmax_dir, nmax, prefixes=None):
        if not has_tarball(directory, nmax):
            return
        if local:
            if not sync_tarball_from_lustre(directory, work, nmax):
                return
        if not unpack_nmax(work, nmax, prefixes=prefixes):
            return

    def finalize_nmax(nmax_dir, nmax):
        target = nmax_dir if os.path.isdir(nmax_dir) else work
        if local:
            cleanup_loose_nmax(target, nmax)
        elif OUTPUT_ROOT.startswith("/lustre"):
            repack_nmax(target, nmax)
        else:
            cleanup_loose_nmax(target, nmax)

    try:
        plot_runs = _init_plot_runs()
        plot_runs.plot_directory(
            work,
            prepare_nmax=prepare_nmax,
            finalize_nmax=finalize_nmax,
            figure_leaf=directory if local else None,
        )
    finally:
        if local:
            remove_work_leaf(directory)
    return "plotted"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-existing", action="store_true", default=True)
    ap.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    ap.add_argument("--directory", default=None, help="Single leaf directory")
    ap.add_argument(
        "--index",
        type=int,
        default=None,
        help="1-based index into complete-leaf list (for SLURM arrays)",
    )
    ap.add_argument("--prune-only", action="store_true",
                    help="Only unpack, prune inactive modes, and repack (no plot)")
    ap.add_argument(
        "--force-reprune",
        action="store_true",
        help="With --prune-only: ignore existing markers and prune again",
    )
    ap.add_argument(
        "--plot-ready-only",
        action="store_true",
        help="Plot only leaves with all Nmax pruned and no figure yet",
    )
    ap.add_argument(
        "--no-require-pruned",
        action="store_true",
        help="Plot complete leaves even if .pruned markers are missing",
    )
    ap.add_argument(
        "--include-partial",
        action="store_true",
        help="Also include leaves missing one or more tarballs (prune/plot what is present)",
    )
    ap.add_argument(
        "--bh-mass",
        default=None,
        help="Only process leaves under this BH tag (e.g. BH_10, BH_1e4)",
    )
    ap.add_argument(
        "--run-tag",
        default=os.environ.get("RUN_TAG", ""),
        help="Subdirectory under output/ and figures/ (e.g. run_nodrag)",
    )
    ap.add_argument(
        "--nl-base",
        default=os.environ.get("NL_BASE_DIR", _HERE),
        help="Nlevels_Runs root (contains output/, figures/, this pipeline)",
    )
    ap.add_argument(
        "--shard-index",
        type=int,
        default=0,
        help="0-based worker index for static sharding (with --shard-count)",
    )
    ap.add_argument(
        "--shard-count",
        type=int,
        default=1,
        help="Number of parallel workers; process leaves where i %% count == index",
    )
    args = ap.parse_args()

    configure_paths(nl_base=args.nl_base, run_tag=args.run_tag)
    bh_mass = args.bh_mass

    if use_local_work():
        _log(f"Using node-local work dir under {_scratch_root()}")

    if args.directory:
        leaves = [os.path.abspath(args.directory)]
    elif args.prune_only:
        # Prune every leaf with at least one unpruned (or stale-marker) tarball.
        leaves = find_pruneable_leaves(OUTPUT_ROOT, bh_mass=bh_mass)
        if args.include_partial:
            # find_pruneable_leaves already includes partials; keep flag for clarity
            pass
    elif args.plot_ready_only:
        if args.no_require_pruned:
            leaves = find_leaves_missing_figures(
                OUTPUT_ROOT,
                bh_mass=bh_mass,
                include_partial=args.include_partial,
            )
        elif args.include_partial:
            # Present tarballs must all be freshly pruned.
            leaves = [
                leaf
                for leaf in find_leaves_missing_figures(
                    OUTPUT_ROOT, bh_mass=bh_mass, include_partial=True
                )
                if is_fully_pruned(leaf)
            ]
        else:
            leaves = find_plottable_leaves(OUTPUT_ROOT, bh_mass=bh_mass)
    else:
        leaves = find_complete_leaves(OUTPUT_ROOT, bh_mass=bh_mass)

    if args.shard_count < 1:
        raise SystemExit("--shard-count must be >= 1")
    if not (0 <= args.shard_index < args.shard_count):
        raise SystemExit("--shard-index must satisfy 0 <= index < shard-count")
    if args.shard_count > 1 and not args.directory:
        before = len(leaves)
        leaves = [
            leaf for i, leaf in enumerate(leaves) if i % args.shard_count == args.shard_index
        ]
        _log(
            f"Shard {args.shard_index}/{args.shard_count}: "
            f"{len(leaves)}/{before} leaves"
        )

    if args.index is not None:
        if args.plot_ready_only or args.prune_only:
            _log("--index ignored with --plot-ready-only/--prune-only")
        else:
            all_leaves = find_complete_leaves(OUTPUT_ROOT, bh_mass=bh_mass)
            if args.index < 1 or args.index > len(all_leaves):
                _log(f"Index {args.index} out of range (1..{len(all_leaves)})")
                sys.exit(0)
            leaves = [all_leaves[args.index - 1]]

    total_complete = (
        len(find_complete_leaves(OUTPUT_ROOT, bh_mass=bh_mass))
        if os.path.isdir(OUTPUT_ROOT)
        else 0
    )
    total_ready = (
        len(find_plottable_leaves(OUTPUT_ROOT, bh_mass=bh_mass))
        if os.path.isdir(OUTPUT_ROOT)
        else 0
    )
    _log(
        f"Processing {len(leaves)} leaf(s) "
        f"(complete: {total_complete}, plottable: {total_ready})"
    )
    stats = {"plotted": 0, "skipped": 0, "errors": 0, "pruned": 0, "not_ready": 0}

    for i, leaf in enumerate(leaves, 1):
        rel = os.path.relpath(leaf, OUTPUT_ROOT)
        _log(f"\n[{i}/{len(leaves)}] {rel}")
        t0 = _time.time()
        try:
            result = process_leaf(
                leaf,
                skip_existing=args.skip_existing and not args.prune_only,
                prune_only=args.prune_only,
                require_pruned=not args.no_require_pruned,
                force_reprune=args.force_reprune,
            )
            stats[result if result in stats else "plotted"] += 1
            _log(f"  done in {_time.time() - t0:.1f}s ({result})")
        except Exception as exc:
            stats["errors"] += 1
            _log(f"  ERROR: {exc}")
            import traceback
            traceback.print_exc()

    _log(f"\nSummary: {stats}")
    if stats["errors"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
