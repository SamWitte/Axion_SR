# Complete 1PN correction of the self-gravity amplitude

Derivation and numbers: `ongoing_tests/sg_relativistic/note_source/sg_source_note.pdf`, Sections 3.6 and 3.7.

The rate tables contain the Newtonian exchange with black-hole recoil; the overlay tables of
`ongoing_tests/sg_relativistic` contain the relativistic t/u-channel exchange. This directory computes the first
post-Newtonian correction analytically, **in the variables of the tables** (Klein-Gordon mode functions, Poisson
equation solved in the Boyer-Lindquist radius), so that it can be compared with the relativistic result term by term.

## The source

In harmonic coordinates (`x = r_BL - GM`, `V_B = GM/x`), for the exchange `b -> 3` acting on level `a`:

    S = -2 mu^2 [1 + 2 (w_a w_4 - mu^2)/mu^2 + 4 V_B(x)] v phi_a  +  gravitomagnetic terms
    v = G int 2 mu^2 [1 + 2 (w_b w_3 - mu^2)/mu^2 + 6 V_B(y)] phi_b phi_3^* [1/|x-y| - w_g^2 |x-y| / 2]

and, for transitions with `|m_b - m_3| = 1`, the recoil term in the frame of the hole

    S_rec = 2 mu^2 [1 + (eps_a + eps_4)/mu + 5/2 V_B(x)] (a . x) phi_a ,   a = grad v(0) + 4 d/dt v_i(0)

plus the pair (s-channel) term, which is local: `1/f_a^2 -> 1/f_a^2 + 12 pi G` (included by `load_rates.jl` whenever
`self_grav = true`; `NL_S_CHANNEL=0` leaves it out).

## Code

| file | what it computes |
|---|---|
| `ongoing_tests/sg_relativistic/exchange_gen.jl` | Newtonian amplitudes with Kerr states; switches `SGX_NEWTON_SHIFT` (radial coordinate), `SGX_NEWTON_WX`, `SGX_NEWTON_WY` (weights of the hole's potential), `SGX_NEWTON_WRX`, `SGX_NEWTON_RECRET` (recoil weight, retardation in the acceleration of the hole). All off by default. |
| `ongoing_tests/sg_vs_si/sg_1pn.jl`, `sg_1pn_gen.jl` | gravitomagnetism and retardation (and the stress, kinetic and G^2 pieces in canonical Schroedinger variables), hydrogenic states; `sg_1pn_gen.jl` allows different incoming levels |
| `ongoing_tests/sg_vs_si/sg_1pn_recoil.jl` | Newtonian recoil term and the gravitomagnetic induction term of the acceleration of the hole |
| `pieces_job.sh` | batch job that runs all of the above for a list of points |
| `complete_1pn.py` | assembles the pieces: coefficients of alpha^2 of `(A - A_NR)/A_N`, term by term, with the pair term, and optionally next to the relativistic exchange |

    sbatch --export=ALL,WORK=<work dir>,JOBS=<job file>,SRC_TESTS=<repository>/src/ongoing_tests pieces_job.sh
    python3 complete_1pn.py <work dir> [--rel <dir of relativistic points>]

`complete_1pn.py` handles identical incoming levels (`N1 x N1 -> N3`); for different ones the two exchanges carry
different frequency factors and are not separated in the stored amplitudes.

## The pieces of `sg_1pn.jl` are in different variables

`sg_1pn.jl` reports corrections to a Newtonian amplitude written with canonical Schroedinger wavefunctions
`psi = sqrt(2 mu) (1 + 2 V_B + eps/2mu) phi` in harmonic coordinates. Its stress, kinetic and G^2 pieces add up to
`3/2 sum(eps)/mu + 2 V_B(x) + 2 V_B(y)` (the Einstein-Infeld-Hoffmann weights), which is the same physics as the
weights above; but their sum is **not** the correction of the tabulated rate. To apply them to the tables one has to
add the change of radial coordinate, `2 V_B(x) + 4 V_B(y)` and `1/2 sum(eps)/mu`, which is what `complete_1pn.py`
does implicitly by working in the variables of the tables.

## Result (spin 0.95, coefficients of alpha^2 of `(A - A_NR)/A_N`, without the pair term)

| channel | complete 1PN | relativistic exchange | pair term |
|---|---|---|---|
| 322 x 322 -> 544 x BH | +0.040 | +0.039 | +0.022 |
| 544 x 544 -> 322 x Inf | +0.042 | +0.041 | +0.011 |
| 211 x 211 -> 322 x BH | -1.09 | -1.06 | +0.061 |
| 322 x 322 -> 211 x Inf | -1.50 | -1.48 | +0.031 |

Without recoil the remainder is of order alpha^4 (-0.17 and -0.10 alpha^4). With recoil it is 0.03-0.05 at
alpha = 0.07, partly an effect of the spin of the hole (it is smaller at spin 0.3 and 0.6); 0.01-0.03 is left at
zero spin and is not explained.
