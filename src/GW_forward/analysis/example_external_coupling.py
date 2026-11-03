"""
Example: couple agn_cloud_resonances.py to an external cloud evolution.

Step 1 (python):  export accretion histories        -> runs/hist_XXXX.json
Step 2 (Julia):   evolve M, a, Mc_i with the same     -> runs/evol_XXXX.csv
                  accretion prescription
Step 3 (python):  inspirals + resonances using the external evolution

Here step 2 is done by the built-in free-field CloudModel as a stand-in, so the
whole chain can be tested before plugging in the Julia output.
"""
import json
import numpy as np
import agn_cloud_resonances as A

p = A.Params(mu_eV=3e-17, l_stars=(2, 3), logM_today=(5.5, 6.5), z_obs_max=1.0,
             occ_threshold=1e-6)
rundir = "runs"

# 1. export histories
A.export_histories(20, p, rundir, seed=7)

# 2. stand-in for the Julia code (delete this block once evol_*.csv come from Julia)
cm = A.CloudModel(mu_eV=p.mu_eV, levels=((2, 1, 1), (3, 2, 2), (4, 3, 3)))
for i in range(20):
    d = json.load(open(f"{rundir}/hist_{i:04d}.json"))
    _, _, track = A.evolve_with_cloud(np.random.default_rng(0), p, cm,
                                      hist=A._hist_from_json(d), return_track=True)
    A.write_evolution_csv(track, f"{rundir}/evol_{i:04d}.csv", d["z_seed"])

# 3. inspirals + resonances conditioned on the evolved cloud
res = A.run_external(rundir, p)
for R in res:
    print(f"{R['id']:3d}  age {R['age_inj']/A.Gyr:5.2f} Gyr  alpha={R['alpha']:.3f}  "
          f"a={R['a']:.3f}  present={[A.ket_ascii(s) for s in R['live_states']]}  "
          f"-> {R['outcome']}")
A.plot_gallery(res, p, "external_gallery.png", n_panels=4)
