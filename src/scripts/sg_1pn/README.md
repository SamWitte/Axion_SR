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
| `ongoing_tests/sg_relativistic/exchange_gen.jl` | Newtonian amplitudes with Kerr states; switches `SGX_NEWTON_SHIFT` (radial coordinate), `SGX_NEWTON_WX`, `SGX_NEWTON_WY` (weights of the hole's potential), `SGX_NEWTON_WRX`, `SGX_NEWTON_RECRET` (recoil weight, retardation in the acceleration of the hole), `SGX_NEWTON_EFAC` (frequency factors, per exchange, with the Kerr frequencies). All off by default. |
| `ongoing_tests/sg_vs_si/sg_1pn.jl`, `sg_1pn_gen.jl` | gravitomagnetism and retardation (and the stress, kinetic and G^2 pieces in canonical Schroedinger variables), hydrogenic states; `sg_1pn_gen.jl` allows different incoming levels |
| `ongoing_tests/sg_vs_si/sg_1pn_recoil.jl` | Newtonian recoil term and the gravitomagnetic induction term of the acceleration of the hole, any incoming levels |
| `pieces_job.sh` | batch job that runs all of the above for a list of points |
| `complete_1pn.py` | assembles the pieces: every term divided by `A_NR` (Newtonian exchange + Newtonian recoil), their total and the change of the rate; optionally next to the relativistic exchange |

    sbatch --export=ALL,WORK=<work dir>,JOBS=<job file>,SRC_TESTS=<repository>/src/ongoing_tests pieces_job.sh
    python3 complete_1pn.py <work dir> [--rel <dir of relativistic points>]

Terms printed by `complete_1pn.py`: gravitomagnetic, hole potential (direct exchange in the harmonic radius with the
weights 4 V_B and 6 V_B), binding energy (frequency factors of the direct exchange), retardation, recoil (the complete
1PN correction of the recoil term), s-channel, total.

## The pieces of `sg_1pn.jl` are in different variables

`sg_1pn.jl` reports corrections to a Newtonian amplitude written with canonical Schroedinger wavefunctions
`psi = sqrt(2 mu) (1 + 2 V_B + eps/2mu) phi` in harmonic coordinates, with the hole held fixed. Its stress, kinetic and
G^2 pieces add up to `3/2 sum(eps)/mu + 2 V_B(x) + 2 V_B(y)` (the Einstein-Infeld-Hoffmann weights), which is the same
physics as the weights above; but their sum is **not** the correction of the tabulated rate: the Newtonian amplitude
they refer to differs from the tabulated one at order alpha^2, and the 1PN correction of the recoil term is missing.
Only its gravitomagnetic and retardation pieces are used here.

## Result (spin 0.95)

Change of the rate relative to the tables at alpha = 0.1, in per cent, complete 1PN against the relativistic exchange
(both with the s-channel term):

| channel | complete 1PN | relativistic |
|---|---|---|
| 211 x 211 -> 322 x BH | +3.33 | +3.12 |
| 322 x 322 -> 211 x Inf | +5.26 | +4.99 |
| 322 x 322 -> 544 x BH | +0.12 | +0.12 |
| 544 x 544 -> 322 x Inf | +0.11 | +0.10 |
| 211 x 211 -> 422 x BH | +3.52 | +3.28 |
| 422 x 422 -> 211 x Inf | +9.57 | +8.87 |
| 322 x 422 -> 211 x Inf | +6.22 | +5.84 |
| 322 x 544 -> 211 x Inf | -17.6 | -16.5 |
| 422 x 544 -> 211 x Inf | +4.5 | -10.7 |
| 422 x 844 -> 211 x Inf | +1.6 | -9.4 |
| 322 x 844 -> 211 x Inf | +1.2 | -4.2 |

In the first eight channels the two agree up to terms of higher order in alpha (the remainder is of order alpha^4
without recoil; with recoil it is partly an effect of the spin of the hole, see the note). In the last three, where the
level acted on by the dipole exchange is 544 or 844 and the exchanged transition is 422 -> 211 or 322 -> 211, they
disagree at order alpha^2 (with the hole held fixed they agree; the difference is in the recoil part). Which of the
two calculations is wrong there is not known (19 November 2026).
