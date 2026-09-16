# Post-process N-level evolutions (prune + plot)

After `run_Nlevels.jl` finishes a leaf (or after you compress its outputs to
`NL_output_Nmax_*.tar.gz`), use this toolkit to:

1. **Prune** inactive modes from `States` / `Modes` (shrinks large files)
2. **Plot** occupation / spin / mass panels for each Nmax
3. Optionally **repack** pruned data back into the tarballs

## Layout

| File | Role |
|------|------|
| `prune_states.py` | Drop inactive modes (`max occupation ≤ 1e-30`) |
| `plot_runs.py` | Make evolution figures for one leaf |
| `plot_output_pipeline.py` | Orchestrate unpack → prune → plot → repack |
| `compress_nmax_output.sh` | Pack loose `Time/Spin/States/Modes/MassBH` → `NL_output_Nmax_N.tar.gz` |
| `verify_archive_inventory.py` | Check which leaves have all expected tarballs |
| `postprocess_leaf.sh` | One-leaf convenience wrapper |

Expected output tree (under this directory by default):

```
Nlevels_Runs/
  output/[RUN_TAG/]BH_*/fa_*/alpha_*/
    NL_output_Nmax_3.tar.gz
    NL_output_Nmax_4.tar.gz
    ...
  figures/[RUN_TAG/]BH_*/fa_*/alpha_*/
    *.pdf / *.png
```

Override the root with `NL_BASE_DIR` or `--nl-base` if your `output/` lives elsewhere.

## Dependencies

- Python 3 with `numpy`, `pandas`, `matplotlib`
- `tar` + `gzip` (or `pigz` if available — faster)

## Quick start (one leaf)

```bash
cd src/Nlevels_Runs

# If you still have loose .dat files, compress first:
# bash compress_nmax_output.sh output/BH_10/fa_1e16/alpha_0.6 15

# Prune + plot one finished leaf:
bash postprocess_leaf.sh output/BH_10/fa_1e16/alpha_0.6
```

Or call the pipeline directly:

```bash
python3 plot_output_pipeline.py --directory output/BH_10/fa_1e16/alpha_0.6
```

## Batch modes

```bash
# Prune every leaf that still has unpruned tarballs
python3 plot_output_pipeline.py --prune-only --run-tag run_nodrag

# Plot leaves that are fully pruned but missing figures
python3 plot_output_pipeline.py --plot-ready-only --run-tag run_nodrag

# Include incomplete leaves (whatever Nmax tarballs exist)
python3 plot_output_pipeline.py --directory output/BH_10/fa_1e12/alpha_0.2 --include-partial

# Parallel workers (static sharding)
python3 plot_output_pipeline.py --prune-only --run-tag run_nodrag \
  --shard-index 0 --shard-count 4
```

## Inventory check

```bash
python3 verify_archive_inventory.py --local output/run_nodrag --list-complete
```

## Notes for shared campaigns

- Prefer keeping evolutions as `NL_output_Nmax_*.tar.gz` on shared disk; prune/plot
  can unpack to local scratch (`$SLURM_TMPDIR` / `$TMPDIR`) when under SLURM.
- Our BH-absorption rates use `BHlmax=2`; unpack `rate_sve/Rates_master_Nmax18_LvrHc.zip`
  (or `Rates.zip`) before evolving.
- Do not commit loose `*_LvrHc_.dat` or huge unpacked `States_*.dat` — keep zips/tarballs.
