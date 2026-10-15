"""Relativistic self-gravity rate tables from results_prod.dat and the Newtonian tables rate_sve/<chan>_LvrHc_SG.dat.

    Gamma_SG -> |rho|^2 Gamma_SG,   Gamma_x -> Re(rho) Gamma_x,   kappa -> Re(rho)/|rho| * kappa,   rho = A_GR / A_N
(A_SI and A_N are in phase: kappa ~ 1 in the Newtonian tables).  Output: tables/<chan>_LvrHc_SGrel.dat with the same
columns (alpha, Gamma_SI, Gamma_SG, Gamma_x, kappa), plus tables/<chan>_rho.dat (alpha, Re rho, Im rho, Re rho_NR).

usage: python3 make_tables.py results_prod.dat <rate_sve dir>
"""
import os
import sys
import numpy as np

res = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
src = sys.argv[2]
os.makedirs("tables", exist_ok=True)
for s1, s3, d in [("211", "322", "BH"), ("322", "544", "BH"), ("322", "211", "Inf"), ("544", "322", "Inf")]:
    chan = f"{s1}_{s1}_{s3}_{d}"
    tab = np.loadtxt(os.path.join(src, chan + "_LvrHc_SG.dat"))
    rows = {float(r[3]): r for r in res if (r[0], r[1], r[2]) == (s1, s3, d)}
    out, rho_out = [], []
    for al, gsi, gsg, gx, kap in tab:
        key = min(rows, key=lambda k: abs(k - al)) if rows else None
        if key is None or abs(key - al) > 1e-8 * max(al, 1):
            continue
        r = rows[key]
        rho = complex(float(r[8]), float(r[9])); rnr = complex(float(r[10]), float(r[11]))
        out.append([al, gsi, abs(rho) ** 2 * gsg, rho.real * gx, kap * rho.real / abs(rho)])
        rho_out.append([al, rho.real, rho.imag, rnr.real])
    # no header: drop-in replacement for rate_sve/<chan>_LvrHc_SG.dat (read with readdlm by load_rates.jl)
    np.savetxt(os.path.join("tables", chan + "_LvrHc_SGrel.dat"), np.array(out), delimiter="\t")
    np.savetxt(os.path.join("tables", chan + "_rho.dat"), np.array(rho_out), header="alpha  Re(A_GR/A_N)  Im(A_GR/A_N)  Re(A_NR/A_N)")
    print(f"{chan}: {len(out)} of {len(tab)} alpha points")
