#!/usr/bin/env python3
"""Compare local/cluster archive inventory against the production grid."""

import argparse
import csv
import itertools
import os
import subprocess
import sys

NMAX_LIST = [3, 4, 5, 6, 7, 8, 15, 18]
MASSBH_VALS = [10.0, 1e4, 1e8]
FA_VALS = [1e18, 1e16, 1e14, 1e12]
DEFAULT_ALPHA_VALS = [0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6]

KNOWN_INCOMPLETE = {
    ("BH_10", "fa_1e12", "alpha_0.2"): [15, 18],
    ("BH_1e4", "fa_1e12", "alpha_0.2"): [15, 18],
    ("BH_10", "fa_1e16", "alpha_1.3"): [7],
    ("BH_10", "fa_1e18", "alpha_1.3"): [7],
}


def bh_tag(mass):
    if mass >= 1e7:
        exp = int(round(__import__("math").log10(mass)))
        return f"BH_1e{exp}"
    if mass >= 1e3:
        exp = int(round(__import__("math").log10(mass)))
        return f"BH_1e{exp}"
    return f"BH_{int(mass)}"


def fa_tag(fa):
    import math
    return f"fa_1e{int(round(math.log10(fa)))}"


def alpha_tag(a):
    return f"alpha_{a:g}"


def expected_leaves(bh_masses=None, alpha_vals=None):
    masses = MASSBH_VALS if bh_masses is None else bh_masses
    alphas = DEFAULT_ALPHA_VALS if alpha_vals is None else alpha_vals
    for mass, fa, alpha in itertools.product(masses, FA_VALS, alphas):
        yield bh_tag(mass), fa_tag(fa), alpha_tag(alpha)


def scan_leaf(root, bh, fa, alpha):
    leaf = os.path.join(root, bh, fa, alpha)
    present = []
    missing = []
    corrupt = []
    for n in NMAX_LIST:
        path = os.path.join(leaf, f"NL_output_Nmax_{n}.tar.gz")
        if os.path.isfile(path):
            if os.path.getsize(path) < 1024:
                corrupt.append(n)
            else:
                present.append(n)
        else:
            missing.append(n)
    if not os.path.isdir(leaf):
        status = "missing_dir"
    elif corrupt:
        status = "corrupt"
    elif len(present) == len(NMAX_LIST):
        status = "complete"
    elif present:
        status = "partial"
    else:
        status = "empty"
    return {
        "leaf": f"{bh}/{fa}/{alpha}",
        "bh": bh,
        "fa": fa,
        "alpha": alpha,
        "status": status,
        "present_nmax": present,
        "missing_nmax": missing,
        "corrupt_nmax": corrupt,
        "tarball_count": len(present),
    }


def scan_root(root, bh_masses=None, alpha_vals=None):
    if not os.path.isdir(root):
        return []
    rows = []
    for bh, fa, alpha in expected_leaves(bh_masses, alpha_vals):
        rows.append(scan_leaf(root, bh, fa, alpha))
    return rows


def list_complete_leaves(root, bh_masses=None, alpha_vals=None):
    return [
        r["leaf"]
        for r in scan_root(root, bh_masses, alpha_vals)
        if r["status"] == "complete"
    ]


def summarize(rows, label):
    complete = sum(1 for r in rows if r["status"] == "complete")
    partial = sum(1 for r in rows if r["status"] == "partial")
    empty = sum(1 for r in rows if r["status"] == "empty")
    missing_dir = sum(1 for r in rows if r["status"] == "missing_dir")
    corrupt = sum(1 for r in rows if r["status"] == "corrupt")
    tarballs = sum(r["tarball_count"] for r in rows)
    print(f"\n=== {label} ===")
    print(
        f"leaves: complete={complete} partial={partial} corrupt={corrupt} "
        f"empty={empty} missing_dir={missing_dir} tarballs={tarballs}"
    )
    for r in rows:
        if r["status"] != "complete":
            key = (r["bh"], r["fa"], r["alpha"])
            known = KNOWN_INCOMPLETE.get(key)
            note = " (known incomplete on cluster)" if known == r["missing_nmax"] else ""
            extra = ""
            if r.get("corrupt_nmax"):
                extra = f" corrupt Nmax={r['corrupt_nmax']}"
            print(
                f"  {r['leaf']}: {r['status']} "
                f"missing Nmax={r['missing_nmax']}{extra}{note}"
            )
    return complete, partial, tarballs


def remote_scan(cluster, remote_root, bh_masses, alpha_vals):
    bh_args = ""
    if bh_masses:
        bh_args = " ".join(f"--bh-mass {bh_tag(m)}" for m in bh_masses)
    alpha_args = ""
    if alpha_vals:
        alpha_args = " ".join(f"--alpha {a:g}" for a in alpha_vals)
    remote_nl = os.path.dirname(os.path.dirname(remote_root))
    script_path = f"{remote_nl}/verify_archive_inventory.py"
    try:
        subprocess.run(
            [
                "ssh",
                cluster,
                f"test -f {script_path} && python3 {script_path} --local {remote_root} {bh_args} {alpha_args} "
                f"|| echo 'remote verify script not deployed; skip'",
            ],
            check=False,
        )
        return True
    except Exception as exc:
        print(f"remote scan skipped: {exc}")
        return False


def write_csv(path, local_rows, remote_rows=None):
    remote_by_leaf = {}
    if remote_rows:
        remote_by_leaf = {r["leaf"]: r for r in remote_rows}
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "leaf",
                "local_status",
                "local_tarballs",
                "local_missing_nmax",
                "remote_status",
                "remote_tarballs",
                "remote_missing_nmax",
            ]
        )
        for lr in local_rows:
            rr = remote_by_leaf.get(lr["leaf"], {})
            w.writerow(
                [
                    lr["leaf"],
                    lr["status"],
                    lr["tarball_count"],
                    " ".join(map(str, lr["missing_nmax"])),
                    rr.get("status", ""),
                    rr.get("tarball_count", ""),
                    " ".join(map(str, rr.get("missing_nmax", []))),
                ]
            )
    print(f"\nWrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--local",
        default=None,
        help="Output root (e.g. cluster_archive/output or cluster_archive/run_nodrag/output)",
    )
    ap.add_argument(
        "--run-tag",
        default=os.environ.get("RUN_TAG", ""),
        help="Run tag; with --local pointing at cluster_archive, uses cluster_archive/{tag}/output",
    )
    ap.add_argument(
        "--remote-root",
        default="/lustre/hpc/astro/spieksma/Files_June10_cluster/Axion_SR/src/Nlevels_Runs/output",
    )
    ap.add_argument("--remote", action="store_true", help="Also scan cluster via SSH")
    ap.add_argument(
        "--cluster",
        default=os.environ.get("CLUSTER", "spieksma@astro03.hpc.ku.dk"),
    )
    ap.add_argument(
        "--bh-mass",
        action="append",
        dest="bh_masses",
        help="Limit to BH tag(s), e.g. BH_1e4 (repeatable)",
    )
    ap.add_argument(
        "--alpha",
        type=float,
        action="append",
        dest="alpha_vals",
        default=None,
        help="Limit to alpha value(s) (repeatable)",
    )
    ap.add_argument("--csv", default=None, help="Write comparison CSV")
    ap.add_argument(
        "--require-complete",
        action="store_true",
        help="Exit 1 unless all scanned leaves are complete",
    )
    ap.add_argument(
        "--list-complete",
        action="store_true",
        help="Print complete leaf paths (one per line) and exit",
    )
    args = ap.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    run_tag = (args.run_tag or "").strip()

    if args.local:
        local_root = args.local
    elif run_tag:
        local_root = os.path.join(repo, "cluster_archive", run_tag, "output")
    else:
        local_root = os.path.join(repo, "cluster_archive", "output")

    remote_root = args.remote_root
    if run_tag and remote_root.endswith("/output"):
        remote_root = os.path.join(remote_root, run_tag)

    mass_map = {"BH_10": 10.0, "BH_1e4": 1e4, "BH_1e8": 1e8}
    bh_filter = None
    if args.bh_masses:
        bh_filter = [mass_map.get(b, b) for b in args.bh_masses]

    alpha_vals = args.alpha_vals

    if args.list_complete:
        if not os.path.isdir(local_root):
            sys.exit(0)
        for leaf in list_complete_leaves(local_root, bh_filter, alpha_vals):
            print(leaf)
        sys.exit(0)

    if not os.path.isdir(local_root):
        print(f"ERROR: local root missing: {local_root}", file=sys.stderr)
        sys.exit(1)

    local_rows = scan_root(local_root, bh_filter, alpha_vals)
    lc, lp, lt = summarize(local_rows, f"local ({local_root})")

    remote_rows = None
    if args.remote:
        remote_rows = scan_root(remote_root, bh_filter, alpha_vals)
        summarize(remote_rows, f"remote ({remote_root})")
        remote_scan(args.cluster, remote_root, bh_filter, alpha_vals)

    if args.csv:
        write_csv(args.csv, local_rows, remote_rows)

    if args.require_complete:
        incomplete = [r for r in local_rows if r["status"] != "complete"]
        if incomplete:
            print(
                f"\nFAIL: {len(incomplete)} leaf(s) not complete",
                file=sys.stderr,
            )
            sys.exit(1)

    print("\nOK")


if __name__ == "__main__":
    main()
