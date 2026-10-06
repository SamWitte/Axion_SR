#!/usr/bin/env python3
"""
agn_cloud_resonances.py
=======================

Monte Carlo baseline: which gravitational-atom resonance does an AGN-disk
companion encounter first?

Per realization
---------------
1. SMBH history. Seed at z_seed, grow by episodic accretion at fixed radiative
   efficiency eps (Eddington ratio lam_on while ON), conditioned so that
   M(z=0) = M_today. Spin follows Bardeen thin-disk accretion (ISCO specific
   energy & angular momentum), either coherent or chaotic (random disk sense
   per episode).
2. Injection. Pick a time inside an active AGN episode within the observing
   window z < z_obs_max. M and a are frozen there (inspiral << Salpeter time).
3. Companion. Seed at r0 drawn from dN/dln r ~ r^r0_slope on [r0_min, r0_max]
   (units of r_g), migrate via Type I / Type II torques in a steady alpha-disk
   fed at the SMBH accretion rate, plus GW reaction. When the episode ends the
   disk is removed and the orbit evolves by GWs only.
4. Resonances. For every populated state in Params.states, all FINE (same n,
   different l) and HYPERFINE (same n, l, different m) transitions allowed by
   the tidal selection rules for multipoles l_stars are computed from the
   hydrogenic + fine + hyperfine spectrum (fixed boson mass mu):
   Omega_res = Delta(omega) / Delta(m). Bohr (Delta n != 0) transitions are
   excluded. Every crossing time is recorded (see crossing_table), and
   plot_gallery shows f_orb(t) against all resonance lines.

Baseline: superradiance OFF (the spectrum is used kinematically only).
Generalized: pass a CloudModel to evolve M, a and level masses self-
consistently with superradiant extraction, GW annihilation, and any user-
supplied extra rates (e.g. axion self-interactions).

Conventions: everything in the BH-spin frame. a >= 0 is the dimensionless
spin; cloud m is measured along the spin; sigma = +1 (-1) for a companion orbit
co- (counter-) rotating with the spin. CHECK the sign convention of the
resonance condition against your reference before relying on co/counter labels.

Units: cgs internally; GM/c^3 and r_g = GM/c^2 for dimensionless quantities.
"""

from dataclasses import dataclass
from math import factorial
import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------- constants
G = 6.674e-8
c = 2.998e10
Msun = 1.989e33
yr = 3.156e7
Myr = 1e6 * yr
Gyr = 1e9 * yr
hbar_eVs = 6.582119569e-16
sigma_T = 6.652e-25
m_p = 1.6726e-24
T_EDD = sigma_T * c / (4 * np.pi * G * m_p)     # ~4.5e8 yr
A_MAX = 0.998                                    # Thorne limit

# flat LCDM (Planck-like)
H0 = 67.7e5 / 3.0857e24
OM, OL = 0.31, 0.69


def t_of_z(z):
    return 2 / (3 * H0 * np.sqrt(OL)) * np.arcsinh(np.sqrt(OL / OM) * (1 + z) ** -1.5)


T0 = t_of_z(0.0)


def z_of_t(t):
    x = np.sinh(1.5 * H0 * np.sqrt(OL) * t) / np.sqrt(OL / OM)
    return x ** (-2 / 3) - 1


# ---------------------------------------------------------------- Kerr helpers
def r_isco(a):
    """ISCO radius [r_g]; a signed (a<0 = retrograde orbit)."""
    z1 = 1 + (1 - a**2) ** (1 / 3) * ((1 + a) ** (1 / 3) + (1 - a) ** (1 / 3))
    z2 = np.sqrt(3 * a**2 + z1**2)
    return 3 + z2 - np.sign(a) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def l_isco(a):
    r = r_isco(a)
    return (r**2 - 2 * a * np.sqrt(r) + a**2) / (r**0.75 * np.sqrt(r**1.5 - 3 * r**0.5 + 2 * a))


def e_isco(a):
    r = r_isco(a)
    return (r**1.5 - 2 * r**0.5 + a) / (r**0.75 * np.sqrt(r**1.5 - 3 * r**0.5 + 2 * a))


def dadlnM_acc(a, s):
    """Bardeen spin evolution, disk sense s=+-1 relative to the spin axis.
    Uses E_isco (not 1-eps) so spin has the physical a->1 limit; eps only sets
    the growth timescale."""
    ad = np.clip(s * a, -A_MAX, A_MAX)
    return s * l_isco(ad) / e_isco(ad) - 2 * a


# Precompute the prograde spin flow a(x), x = ln M (autonomous ODE).
def _build_spin_flow():
    ev = lambda x, a: a[0] - A_MAX
    ev.terminal, ev.direction = True, 1
    sol = solve_ivp(lambda x, a: [dadlnM_acc(a[0], 1)], (0, 20), [-A_MAX],
                    events=ev, max_step=0.01, rtol=1e-9, atol=1e-12)
    return sol.t, sol.y[0]


_X_FLOW, _A_FLOW = _build_spin_flow()


def spin_after(a0, s, dlnM):
    """Spin after growing by dlnM with disk sense s (exact flow map)."""
    ad = s * a0
    x1 = np.interp(ad, _A_FLOW, _X_FLOW) + dlnM
    ad1 = A_MAX if x1 >= _X_FLOW[-1] else np.interp(x1, _X_FLOW, _A_FLOW)
    return s * ad1


# ---------------------------------------------------------------- parameters
@dataclass
class Params:
    # boson / cloud spectrum (baseline uses it kinematically)
    mu_eV: float = 1.3e-17              # alpha ~ 0.1 at 1e6 Msun
    # states known to be populated, |n l m> with m along the BH spin
    states: tuple = ((2, 1, 1), (3, 2, 2), (3, 1, 1), (3, 0, 0))
    kinds: tuple = ("fine", "hyperfine")   # Bohr (Delta n != 0) excluded
    final_in_set_only: bool = False     # True: only transitions between listed states
    l_stars: tuple = (2,)               # tidal multipoles kept (2 = quadrupole)
    equatorial: bool = True             # companion orbit in the spin plane:
                                        # tidal field has only l*+m* even
    require_superradiant: bool = True   # drop states with omega > m Omega_H
    occ_threshold: float = 1e-6         # external mode: state present if Mc/M >= this
    # SMBH growth
    eps: float = 0.1
    lam_on: float = 0.3
    t_ep: float = 1e7 * yr
    z_seed: float = 15.0
    logM_seed: tuple = (2.0, 5.0)       # log-uniform range [Msun]
    logM_today: tuple = (5.0, 7.0)
    a_seed: tuple = (0.0, 0.5)
    accretion: str = "coherent"         # or "chaotic"
    z_obs_max: float = 3.0
    # disk
    alpha_ss: float = 0.01
    h: float = 0.01
    # companion
    m_comp_Msun: float = 10.0
    r0_min: float = 1e3                 # [r_g]
    r0_max: float = 1e5
    r0_slope: float = 0.0               # dN/dlnr ~ r^slope (0 = log-uniform)
    t_max_inspiral: float = 1 * Gyr
    n_grid: int = 3000


# ---------------------------------------------------------------- SMBH history
def make_history(rng, p):
    """Episodic growth conditioned on M_today. Returns dict or None."""
    t_seed = t_of_z(p.z_seed)
    T = T0 - t_seed
    lnMs = np.log(10 ** rng.uniform(*p.logM_seed) * Msun)
    lnMt = np.log(10 ** rng.uniform(*p.logM_today) * Msun)
    if lnMt <= lnMs:
        return None
    lam_bar = p.eps / (1 - p.eps) * T_EDD * (lnMt - lnMs) / T
    if lam_bar >= p.lam_on:
        return None
    n_ep = max(1, int(round(T * lam_bar / (p.lam_on * p.t_ep))))
    t_ep = T * lam_bar / (p.lam_on * n_ep)        # exact duty cycle
    P = T / n_ep
    starts = t_seed + np.arange(n_ep) * P + rng.uniform(0, P - t_ep, n_ep)
    dlnM_ep = (lnMt - lnMs) / n_ep
    s = np.ones(n_ep) if p.accretion == "coherent" else rng.choice([-1.0, 1.0], n_ep)
    a = np.empty(n_ep + 1)
    a[0] = rng.uniform(*p.a_seed)
    for k in range(n_ep):
        a[k + 1] = spin_after(a[k], s[k], dlnM_ep)
    return dict(starts=starts, t_ep=t_ep, s=s, a=a, lnMs=lnMs, lnMt=lnMt,
                dlnM_ep=dlnM_ep, lam_bar=lam_bar, n_ep=n_ep)


def draw_injection(rng, hist, p):
    """Uniform in active time within the observing window."""
    t_lo = max(t_of_z(p.z_obs_max), hist["starts"][0])
    lo = np.maximum(hist["starts"], t_lo)
    hi = np.minimum(hist["starts"] + hist["t_ep"], T0)
    w = np.clip(hi - lo, 0, None)
    if w.sum() == 0:
        return None
    k = rng.choice(len(w), p=w / w.sum())
    t_inj = rng.uniform(lo[k], hi[k])
    frac = (t_inj - hist["starts"][k]) / hist["t_ep"]
    lnM = hist["lnMs"] + (k + frac) * hist["dlnM_ep"]
    a = spin_after(hist["a"][k], hist["s"][k], frac * hist["dlnM_ep"])
    t_disk = hist["starts"][k] + hist["t_ep"] - t_inj
    return dict(t_inj=t_inj, M=np.exp(lnM), a=a, s=hist["s"][k], t_disk=t_disk)


# ---------------------------------------------------------------- disk + inspiral
def rdot(r, M, m, p, disk_on=True):
    """Radial velocity [cm/s] (negative = inward)."""
    Om = np.sqrt(G * M / r**3)
    beta = 64 / 5 * G**3 * M * m * (M + m) / c**5
    v = -beta / r**3
    if not disk_on:
        return v
    q, h, al = m / M, p.h, p.alpha_ss
    Mdot0 = p.lam_on * M / (p.eps * T_EDD)
    nu = al * h**2 * r**2 * Om
    Sigma = Mdot0 / (3 * np.pi * nu)
    crida = 0.75 * h * (3 / q) ** (1 / 3) + 50 * al * h**2 / q
    C_I = 1.364 + 0.541 * 0.5                       # Tanaka+02, Sigma ~ r^-1/2
    v_I = -2 * C_I * q * h**-2 * Sigma * r**3 * Om / M
    v_II = -1.5 * nu / r * np.minimum(1.0, 4 * np.pi * Sigma * r**2 / m)
    return v + np.where(crida > 1, v_I, v_II)


def inspiral(M, m, r0, r_stop, t_disk, p):
    """t(r) on a log-r grid; disk on for t < t_disk, GW-only after.
    Returns t [s], r [cm], stalled flag (t_max reached before r_stop)."""
    beta = 64 / 5 * G**3 * M * m * (M + m) / c**5
    x = np.linspace(np.log(r0), np.log(r_stop), p.n_grid)
    r = np.exp(x)
    dtdx = r / np.abs(rdot(r, M, m, p, True))
    t = np.concatenate([[0], np.cumsum(0.5 * (dtdx[1:] + dtdx[:-1]) * np.abs(np.diff(x)))])
    if t[-1] > t_disk:                                 # disk switches off
        j = np.searchsorted(t, t_disk)
        x_off = x[j - 1] + (x[j] - x[j - 1]) * (t_disk - t[j - 1]) / (t[j] - t[j - 1])
        r_off = np.exp(x_off)
        r = np.concatenate([r[:j], [r_off], r[j:]])
        t = np.concatenate([t[:j], [t_disk], t_disk + (r_off**4 - r[j + 1:]**4) / (4 * beta)])
    stalled = t[-1] > p.t_max_inspiral
    if stalled:
        j = np.searchsorted(t, p.t_max_inspiral)
        r4 = r[j - 1]**4 - 4 * beta * (p.t_max_inspiral - t[j - 1]) if t[j - 1] >= t_disk \
            else np.interp(p.t_max_inspiral, t, r)**4
        r = np.concatenate([r[:j], [max(r4, r_stop**4) ** 0.25]])
        t = np.concatenate([t[:j], [p.t_max_inspiral]])
    return t, r, stalled


# ---------------------------------------------------------------- cloud spectrum
def eps_level(n, l, m, alpha, a):
    """(omega_nlm/mu - 1): hydrogenic + fine + hyperfine."""
    f = -6 / (2 * l + 1) + 2 / n
    hl = 16 / (2 * l * (2 * l + 1) * (2 * l + 2)) if l > 0 else 0.0
    return (-alpha**2 / (2 * n**2) - alpha**4 / (8 * n**4)
            + f * alpha**4 / n**3 + hl * a * m * alpha**5 / n**3)


def tidal_allowed(l, lp, dm, l_stars):
    return any(abs(l - lp) <= ls <= l + lp and (l + lp + ls) % 2 == 0 and abs(dm) <= ls
               for ls in l_stars)


def transitions(alpha, a, p):
    """Fine (same n, different l) and hyperfine (same n, l; different m)
    transitions out of every state in p.states, passing the tidal selection
    rules for p.l_stars. Returns list of dicts with
      init, final, kind, Omega (c^3/GM, signed: >0 co-, <0 counter-rotating),
      r_res (r_g), l_star (lowest multipole that allows it)."""
    out = []
    for (n, l, m) in p.states:
        e0 = eps_level(n, l, m, alpha, a)
        for l2 in range(n):
            for m2 in range(-l2, l2 + 1):
                dm = m2 - m
                kind = "fine" if l2 != l else "hyperfine"
                if dm == 0 or kind not in p.kinds:
                    continue
                if p.final_in_set_only and (n, l2, m2) not in p.states:
                    continue
                ls = [L for L in p.l_stars if abs(l - l2) <= L <= l + l2
                      and (l + l2 + L) % 2 == 0 and abs(dm) <= L
                      and (not p.equatorial or (L + dm) % 2 == 0)]
                if not ls:
                    continue
                Om = alpha * (eps_level(n, l2, m2, alpha, a) - e0) / dm
                if Om == 0:
                    continue
                out.append(dict(init=(n, l, m), final=(n, l2, m2), kind=kind,
                                Omega=Om, r_res=abs(Om) ** (-2 / 3), l_star=min(ls)))
    return out


def ket(s):
    return "|" + "".join(str(x) if x >= 0 else "̅" + str(-x) for x in s) + "⟩"


def ket_ascii(s):
    return "".join(str(x) if x >= 0 else "m" + str(-x) for x in s)


def alpha_of(M, mu_eV):
    return G * M / c**3 * mu_eV / hbar_eVs


# ---------------------------------------------------------------- one realization
def realize(rng, p, cloud_model=None):
    if cloud_model is None:
        hist = make_history(rng, p)
        if hist is None:
            return None
        inj = draw_injection(rng, hist, p)
    else:
        hist, inj = evolve_with_cloud(rng, p, cloud_model)
    if inj is None:
        return None
    return realize_from_injection(rng, p, inj)


def realize_from_injection(rng, p, inj, populated=None):
    """Companion inspiral + resonance bookkeeping for a given injection.

    populated: optional callable t_cosmic[s] -> tuple of states present at
    that time (e.g. from an external cloud evolution). If given, it replaces
    p.states and the superradiance filter: a transition counts as crossed only
    if its initial state is populated at the crossing time."""
    M, a_sgn = inj["M"], inj["a"]
    a = abs(a_sgn)
    sigma = 1.0 if a_sgn == 0 else float(np.sign(inj["s"] * a_sgn))
    rg = G * M / c**2
    tg = G * M / c**3
    m = p.m_comp_Msun * Msun
    u = rng.uniform()
    if p.r0_slope == 0:
        r0 = p.r0_min * (p.r0_max / p.r0_min) ** u
    else:
        k = p.r0_slope
        r0 = (p.r0_min**k + u * (p.r0_max**k - p.r0_min**k)) ** (1 / k)
    r_stop = r_isco(sigma * min(a, A_MAX))
    t, r, stalled = inspiral(M, m, r0 * rg, r_stop * rg, inj["t_disk"], p)
    alpha = alpha_of(M, p.mu_eV)
    r_end = r[-1] / rg
    from dataclasses import replace
    if populated is not None:
        t_samp = inj["t_inj"] + np.linspace(t[0], t[-1], 200)
        live = tuple(sorted({st for ts in t_samp for st in populated(ts)}))
        sr_ok = {st: _sr_condition(st, alpha, a) for st in live}
        trans = transitions(alpha, a, replace(p, states=live)) if live else []
    else:
        sr_ok = {st: _sr_condition(st, alpha, a) for st in p.states}
        live = tuple(st for st in p.states if sr_ok[st]) if p.require_superradiant else p.states
        trans = transitions(alpha, a, replace(p, states=live)) if live else []
    for T in trans:
        T["f_res"] = abs(T["Omega"]) / tg / (2 * np.pi)          # Hz
        T["matches_orbit"] = np.sign(T["Omega"]) == sigma
        T["ahead"] = T["r_res"] < r0
        reached = bool(T["matches_orbit"] and T["ahead"] and r_end <= T["r_res"])
        T["t_cross"] = (np.interp(-np.log(T["r_res"] * rg), -np.log(r), t)
                        if reached else np.nan)
        T["init_populated"] = (True if populated is None or not reached
                               else T["init"] in populated(inj["t_inj"] + T["t_cross"]))
        T["crossed"] = reached and T["init_populated"]
    crossed = sorted([T for T in trans if T["crossed"]], key=lambda T: T["t_cross"])
    first = crossed[0] if crossed else None
    if not live:
        outcome = "no_SR_state"
    elif first is not None:
        outcome = first["kind"]
    elif any(T["matches_orbit"] and T["ahead"] for T in trans):
        outcome = "stalled"
    else:
        outcome = "none_ahead"
    Omega = np.sqrt(G * M / r**3)
    t_seed = t_of_z(p.z_seed)
    return dict(age_inj=inj["t_inj"] - t_seed, live_states=live, M=M, a=a, sigma=sigma,
                alpha=alpha, r0=r0, r_end=r_end, stalled=stalled, outcome=outcome,
                first=first, crossed=crossed, transitions=trans, t=t,
                f_orb=Omega / (2 * np.pi), sr_ok=sr_ok, t_inj=inj["t_inj"],
                cloud=inj.get("cloud"))


def _sr_condition(state, alpha, a):
    """Is |nlm> superradiantly unstable (omega < m Omega_H)?"""
    n, l, m = state
    w = alpha * (1 + eps_level(n, l, m, alpha, a))
    return bool(m > 0 and w < m * a / (2 * (1 + np.sqrt(1 - a**2))))


def crossing_table(results):
    """Flat list of rows (one per crossing) for saving / inspection."""
    rows = []
    for i, R in enumerate(results):
        for T in R["crossed"]:
            rows.append(dict(realization=i, M_Msun=R["M"] / Msun, a=R["a"], alpha=R["alpha"],
                             sigma=R["sigma"], init=ket_ascii(T["init"]),
                             final=ket_ascii(T["final"]), kind=T["kind"], l_star=T["l_star"],
                             r_res_rg=T["r_res"], f_res_Hz=T["f_res"],
                             t_cross_yr=T["t_cross"] / yr,
                             age_cross_Gyr=(R["age_inj"] + T["t_cross"]) / Gyr))
    return rows


def run_monte_carlo(N, p, seed=0, cloud_model=None):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < N:
        R = realize(rng, p, cloud_model)
        if R is not None:
            out.append(R)
    return out


# ======================================================================
# Generalized: self-consistent BH + cloud evolution (superradiance ON)
# ======================================================================
@dataclass
class CloudModel:
    """levels: tracked |n l m> (m along the spin). extra_rates(t, M, a, Mc, alpha)
    may return (dMc_dt array [g/s], dM_dt [g/s], dJ_dt [g cm^2/s^2]) for
    processes such as axion self-interactions (e.g. 211x211 -> 322 x BH,
    322x322 -> 211 x inf). Coefficients are NOT supplied here: take them
    from the literature and verify them."""
    mu_eV: float = 1.3e-17
    levels: tuple = ((2, 1, 1), (3, 2, 2))
    gw_annihilation_211: bool = True
    extra_rates: object = None
    floor_frac: float = 1e-60


def gamma_sr(n, l, m, alpha, a):
    """Superradiant field-amplitude rate x GM/c^3 (Detweiler, with the factor-2
    correction as in Baumann et al. 2019). Occupation grows at 2*Gamma."""
    rp = 1 + np.sqrt(1 - a**2)
    Om_H = a / (2 * rp)
    w = alpha * (1 + eps_level(n, l, m, alpha, a))
    C = (2 ** (4 * l + 1) * factorial(n + l) / (n ** (2 * l + 4) * factorial(n - l - 1))
         * (factorial(l) / (factorial(2 * l) * factorial(2 * l + 1))) ** 2)
    g = np.prod([k**2 * (1 - a**2) + (a * m - 2 * rp * w) ** 2 for k in range(1, l + 1)])
    return 2 * rp * C * g * (m * Om_H - w) * alpha ** (4 * l + 5)


def _rhs(t, y, s_acc, p, cm):
    lnM, a = y[0], np.clip(y[1], -A_MAX, A_MAX)
    lnMc = np.minimum(y[2:], lnM)          # guard against trial-step overflow
    M = np.exp(lnM)
    tg = G * M / c**3
    alpha = alpha_of(M, cm.mu_eV)
    Mc = np.exp(lnMc)
    dlnM = dlna_acc = 0.0
    if s_acc != 0:
        dlnM = p.lam_on * (1 - p.eps) / (p.eps * T_EDD)
        dlna_acc = dadlnM_acc(a, s_acc) * dlnM
    dlnMc = np.zeros_like(Mc)
    dM_sr = dJ_sr = 0.0
    for i, (n, l, m) in enumerate(cm.levels):
        G_i = gamma_sr(n, l, m, alpha, abs(a)) / tg
        w_hat = alpha * (1 + eps_level(n, l, m, alpha, abs(a)))
        dlnMc[i] = 2 * G_i
        dM_sr -= 2 * G_i * Mc[i]
        dJ_sr -= (m / w_hat) * 2 * G_i * Mc[i] * (G * M / c)        # g cm^2/s
        if cm.gw_annihilation_211 and (n, l, m) == (2, 1, 1):
            C_gw = (484 + 9 * np.pi**2) / 23040
            dlnMc[i] -= C_gw * (Mc[i] / M) * alpha**14 / tg
    if cm.extra_rates is not None:
        dMc_x, dM_x, dJ_x = cm.extra_rates(t, M, a, Mc, alpha)
        dlnMc += np.asarray(dMc_x) / Mc
        dM_sr += dM_x
        dJ_sr += dJ_x
    floor = np.log(cm.floor_frac * M)
    dlnMc = np.where((lnMc <= floor) & (dlnMc < 0), 0.0, dlnMc)
    dlnM_tot = dlnM + dM_sr / M
    da = dlna_acc + c * dJ_sr / (G * M**2) - 2 * a * dM_sr / M
    if abs(y[1]) >= A_MAX and np.sign(da) == np.sign(y[1]):
        da = 0.0                                                   # Thorne cap
    return np.concatenate([[dlnM_tot, da], dlnMc])


def evolve_with_cloud(rng, p, cm, hist=None, return_track=False):
    """Integrate the full history with superradiance. Mass conditioning ignores
    the (small) SR mass loss."""
    if hist is None:
        hist = None
        while hist is None:
            hist = make_history(rng, p)
    t_seed = t_of_z(p.z_seed)
    ep_s, ep_e = hist["starts"], hist["starts"] + hist["t_ep"]
    edges = np.sort(np.concatenate([[t_seed, T0], ep_s, ep_e]))
    M0 = np.exp(hist["lnMs"])
    y = np.concatenate([[hist["lnMs"], hist["a"][0]],
                        np.full(len(cm.levels), np.log(1e3 * cm.floor_frac * M0))])
    ts, ys = [], []
    for t1, t2 in zip(edges[:-1], edges[1:]):
        if t2 <= t1:
            continue
        k = np.searchsorted(ep_s, 0.5 * (t1 + t2)) - 1
        on = k >= 0 and 0.5 * (t1 + t2) < ep_e[k]
        s_acc = hist["s"][k] if on else 0
        sol = solve_ivp(_rhs, (t1, t2), y, args=(s_acc, p, cm), method="Radau",
                        rtol=1e-6, atol=1e-8, dense_output=True)
        ts.append(sol.t)
        ys.append(sol.y)
        y = sol.y[:, -1].copy()
        y[1] = np.clip(y[1], -A_MAX, A_MAX)
    T = np.concatenate(ts)
    Y = np.concatenate(ys, axis=1)
    track = dict(t=T, M=np.exp(Y[0]), a=Y[1], Mc=np.exp(Y[2:]), levels=cm.levels)
    inj = draw_injection(rng, hist, p)
    if inj is not None:
        i = np.searchsorted(T, inj["t_inj"])
        inj["M"], inj["a"] = track["M"][i], track["a"][i]
        inj["cloud"] = {lv: track["Mc"][j, i] / track["M"][i] for j, lv in enumerate(cm.levels)}
    return (hist, inj, track) if return_track else (hist, inj)


# ======================================================================
# Coupling to an external cloud evolution (e.g. a Julia code)
# ======================================================================
#
# Workflow
#   1. export_histories(N, p, "runs/")  -> runs/hist_0000.json, ...
#      Each JSON fixes the accretion history exactly (seed mass & spin, fixed
#      eps, Eddington ratio while ON, episode list with disk sense s). The same
#      content is written as runs/hist_0000.txt for codes without a JSON reader
#      ("# key = value" header lines, then rows: start_age_yr end_age_yr s).
#      The Julia side is src/GW_forward/analysis/evolve_agn_history.jl, which
#      also writes runs/evol_0000_meta.json (mu_eV, f_a, Nmax, status); if
#      present, run_external checks mu_eV against p.mu_eV.
#   2. Your code reads each JSON, evolves M, a and the level populations
#      under the SAME accretion prescription, and writes runs/evol_0000.csv:
#          age_yr, M_Msun, a, Mc_2_1_1, Mc_3_2_2, Mc_2_1_m1, ...
#      Mc_* = cloud mass of that level divided by M (dimensionless);
#      m < 0 written as m1, m2, ... ; age measured from seed formation.
#      Rows must resolve the times you care about (interpolation is linear
#      in age for M, a and in log for Mc/M).
#   3. run_external("runs/", p) draws injections and inspirals from the same
#      histories, takes M and a from your evolution at the injection time, and
#      treats a state as present only while Mc/M >= p.occ_threshold.
#
# Accretion prescription to reproduce (written into every JSON as well):
#   ON  : dlnM/dt = lam_on (1 - eps) / (eps t_Edd)          [fixed eps]
#         da/dlnM = s l_isco(s a) / E_isco(s a) - 2 a        [Bardeen, s = +-1]
#   OFF : no accretion.   Spin capped at 0.998.

import json
import os


def export_histories(N, p, outdir, seed=0):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    t_seed = t_of_z(p.z_seed)
    i = 0
    while i < N:
        hist = make_history(rng, p)
        if hist is None:
            continue
        d = dict(
            id=i, rng_seed=int(rng.integers(2**63)),
            z_seed=p.z_seed, age_today_yr=(T0 - t_seed) / yr,
            M_seed_Msun=float(np.exp(hist["lnMs"]) / Msun),
            M_today_target_Msun=float(np.exp(hist["lnMt"]) / Msun),
            a_seed=float(hist["a"][0]), eps=p.eps, lam_on=p.lam_on,
            t_Edd_yr=T_EDD / yr,
            dlnM_dt_on_per_yr=p.lam_on * (1 - p.eps) / (p.eps * T_EDD) * yr,
            spin_cap=A_MAX,
            episodes=[dict(start_age_yr=float((t0 - t_seed) / yr),
                           end_age_yr=float((t0 + hist["t_ep"] - t_seed) / yr),
                           s=int(sk))
                      for t0, sk in zip(hist["starts"], hist["s"])],
            prescription="ON: dlnM/dt = lam_on(1-eps)/(eps t_Edd); "
                         "da/dlnM = s*l_isco(s*a)/E_isco(s*a) - 2a. OFF: none.")
        with open(os.path.join(outdir, f"hist_{i:04d}.json"), "w") as fh:
            json.dump(d, fh, indent=1)
        with open(os.path.join(outdir, f"hist_{i:04d}.txt"), "w") as fh:
            for k, v in d.items():
                if k not in ("episodes", "prescription"):
                    fh.write(f"# {k} = {v!r}\n")
            fh.write("# columns: start_age_yr end_age_yr s\n")
            for e in d["episodes"]:
                fh.write(f"{e['start_age_yr']!r} {e['end_age_yr']!r} {e['s']}\n")
        i += 1


def _hist_from_json(d):
    t_seed = t_of_z(d["z_seed"])
    starts = np.array([e["start_age_yr"] for e in d["episodes"]]) * yr + t_seed
    t_ep = (d["episodes"][0]["end_age_yr"] - d["episodes"][0]["start_age_yr"]) * yr
    s = np.array([e["s"] for e in d["episodes"]], float)
    lnMs, lnMt = np.log(d["M_seed_Msun"] * Msun), np.log(d["M_today_target_Msun"] * Msun)
    dlnM_ep = (lnMt - lnMs) / len(starts)
    a = np.empty(len(starts) + 1)
    a[0] = d["a_seed"]
    for k in range(len(starts)):
        a[k + 1] = spin_after(a[k], s[k], dlnM_ep)
    return dict(starts=starts, t_ep=t_ep, s=s, a=a, lnMs=lnMs, lnMt=lnMt,
                dlnM_ep=dlnM_ep, n_ep=len(starts))


def _parse_state(col):
    n, l, m = col[3:].split("_")
    return (int(n), int(l), -int(m[1:]) if m.startswith("m") else int(m))


def load_evolution(path, z_seed):
    """Read an external evolution CSV -> interpolators in cosmic time [s]."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    cols = data.dtype.names
    t = data["age_yr"] * yr + t_of_z(z_seed)
    states = {_parse_state(cname): np.log(np.maximum(data[cname], 1e-300))
              for cname in cols if cname.startswith("Mc_")}
    return dict(
        t=t,
        M=lambda tt: np.exp(np.interp(tt, t, np.log(data["M_Msun"] * Msun))),
        a=lambda tt: np.interp(tt, t, data["a"]),
        frac=lambda tt: {st: float(np.exp(np.interp(tt, t, lf))) for st, lf in states.items()},
    )


def run_external(rundir, p, check_mass=0.2, n_inj=1):
    """Inspirals + resonances using external M(t), a(t), Mc_i(t).
    n_inj companions are dropped into each history (independent injection
    times and initial radii; the first one is the same as with n_inj=1)."""
    results = []
    for fn in sorted(f for f in os.listdir(rundir) if f.startswith("hist_") and f.endswith(".json")):
        idx = fn[5:9]
        evol_path = os.path.join(rundir, f"evol_{idx}.csv")
        if not os.path.exists(evol_path):
            continue
        meta_path = os.path.join(rundir, f"evol_{idx}_meta.json")
        if os.path.exists(meta_path):
            meta = json.load(open(meta_path))
            if abs(meta["mu_eV"] / p.mu_eV - 1) > 1e-6:
                raise ValueError(f"{meta_path}: evolution used mu = {meta['mu_eV']} eV, "
                                 f"but p.mu_eV = {p.mu_eV}")
            if meta.get("status", "ok") != "ok":
                print(f"[warn] {idx}: external evolution status '{meta['status']}'; skipped")
                continue
        d = json.load(open(os.path.join(rundir, fn)))
        hist = _hist_from_json(d)
        ev = load_evolution(evol_path, d["z_seed"])
        rng = np.random.default_rng(d["rng_seed"])
        populated = lambda tt, ev=ev: tuple(st for st, f in ev["frac"](tt).items()
                                            if f >= p.occ_threshold)
        for k in range(n_inj):
            inj = draw_injection(rng, hist, p)
            if inj is None:
                break
            M_py = inj["M"]
            inj["M"], inj["a"] = float(ev["M"](inj["t_inj"])), float(ev["a"](inj["t_inj"]))
            if abs(np.log(inj["M"] / M_py)) > check_mass:
                print(f"[warn] {idx}: external M differs from accretion model by "
                      f"{inj['M'] / M_py:.2f}x at injection (different accretion history?)")
            R = realize_from_injection(rng, p, inj, populated=populated)
            R["id"] = int(idx)
            R["inj"] = k
            R["cloud"] = ev["frac"](inj["t_inj"])
            results.append(R)
    return results


def write_evolution_csv(track, path, z_seed):
    """Write a CloudModel track in the external format (stand-in / test)."""
    age = (track["t"] - t_of_z(z_seed)) / yr
    cols = ["age_yr", "M_Msun", "a"] + [
        "Mc_" + "_".join(str(x) if x >= 0 else f"m{-x}" for x in lv) for lv in track["levels"]]
    arr = np.column_stack([age, track["M"] / Msun, track["a"]] +
                          [track["Mc"][j] / track["M"] for j in range(len(track["levels"]))])
    np.savetxt(path, arr, delimiter=",", header=",".join(cols), comments="")


# ======================================================================
# Plotting / demo
# ======================================================================
KIND_STYLE = {"fine": dict(color="C2", ls="-"), "hyperfine": dict(color="C3", ls="--")}
OUTCOME_COL = {"hyperfine": "C3", "fine": "C2", "stalled": "0.5",
               "none_ahead": "0.8", "no_SR_state": "0.3"}


def _group_lines(trans, rtol=1e-6):
    """Merge transitions at the same resonant frequency (e.g. Delta m=-1,-2
    hyperfine lines are degenerate at this order)."""
    groups = []
    for T in sorted(trans, key=lambda T: T["f_res"]):
        if groups and abs(T["f_res"] / groups[-1][0]["f_res"] - 1) < rtol:
            groups[-1].append(T)
        else:
            groups.append([T])
    return groups


def plot_track_with_resonances(R, ax, show_opposite=False, label_size=7):
    """One realization: f_orb(t) with every fine / hyperfine line lying ahead
    of the companion; crossings marked; labels staggered in the right margin."""
    from matplotlib.transforms import blended_transform_factory
    tt = R["t"] / yr                     # yr after injection; axis labelled by BH age
    f_lo, f_hi = R["f_orb"][0], R["f_orb"][-1]
    trans = [T for T in R["transitions"] if T["matches_orbit"] or show_opposite]
    groups = [g for g in _group_lines(trans) if f_lo <= g[0]["f_res"] <= f_hi]
    f_top = max([g[0]["f_res"] for g in groups], default=f_lo)
    ylo, yhi = 0.6 * f_lo, min(f_hi, 3 * f_top)
    ax.semilogy(tt, R["f_orb"], color="k", lw=1.2)
    ax.set_ylim(ylo, yhi)
    # stagger label heights (log space) to avoid overlaps
    span = np.log10(yhi / ylo)
    gap = min(0.07 * span, span / (len(groups) + 1)) if groups else 0
    ylab, prev = [], -np.inf
    for g in groups:
        y = max(np.log10(g[0]["f_res"]), prev + gap)
        ylab.append(y)
        prev = y
    if ylab and ylab[-1] > np.log10(yhi):                 # shift block down if needed
        ylab = list(np.array(ylab) - (ylab[-1] - np.log10(yhi)))
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    for g, yl in zip(groups, ylab):
        f = g[0]["f_res"]
        st = dict(KIND_STYLE[g[0]["kind"]])
        if not g[0]["matches_orbit"]:
            st.update(alpha=0.35, ls=":")
        ax.axhline(f, lw=0.8, **st)
        lab = ", ".join(f"{ket(T['init'])}\u2192{ket(T['final'])}" for T in g)
        ax.annotate(lab, xy=(1.0, f), xycoords=tr, xytext=(1.03, 10**yl), textcoords=tr,
                    fontsize=label_size, color=st["color"], va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=st["color"]),
                    annotation_clip=False)
        for T in g:
            if T["crossed"]:
                ax.plot(T["t_cross"] / yr, f, "o", ms=5, color=st["color"],
                        mec="k", mew=0.5, zorder=5)
    z_inj = z_of_t(R["t_inj"])
    ax.set_title(f"age {R['age_inj']/Gyr:.2f} Gyr (z={z_inj:.2f}), M={R['M']/Msun:.1e} M\u2609, "
                 f"a={R['a']:.3f}, \u03b1={R['alpha']:.3f}, "
                 f"{'co' if R['sigma'] > 0 else 'counter'}-rot.", fontsize=9)


def plot_gallery(results, p, fname, n_panels=6, require_crossing=True):
    import matplotlib.pyplot as plt
    pick = [R for R in results if (R["crossed"] or not require_crossing)][:n_panels]
    if not pick:
        print("plot_gallery: no realization with a crossing; try require_crossing=False")
        return
    nc = 1
    nr = len(pick)
    fig, axs = plt.subplots(nr, nc, figsize=(9, 3.4 * nr), squeeze=False)
    for ax, R in zip(axs.flat, pick):
        plot_track_with_resonances(R, ax)
        ax.set_xlabel(f"black-hole age \u2212 {R['age_inj']/Gyr:.4f} Gyr  [yr]", fontsize=8)
        ax.set_ylabel(r"$f_{\rm orb}$ [Hz]", fontsize=8)
    for ax in axs.flat[len(pick):]:
        ax.axis("off")
    axs.flat[0].plot([], [], **KIND_STYLE["fine"], label="fine")
    axs.flat[0].plot([], [], **KIND_STYLE["hyperfine"], label="hyperfine")
    axs.flat[0].legend(fontsize=8, loc="upper left")
    fig.suptitle(f"States: {', '.join(ket(s) for s in p.states)}   "
                 f"μ={p.mu_eV:.1e} eV,  l*∈{p.l_stars}", fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.78, 0.98))
    fig.savefig(fname, dpi=140)


def plot_baseline(results, p, fname, n_tracks=40):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw=dict(width_ratios=[2, 1]))
    for R in results[:n_tracks]:
        tt = (R["age_inj"] + R["t"]) / Gyr
        ax[0].semilogy(tt, R["f_orb"], color=OUTCOME_COL[R["outcome"]], lw=0.8, alpha=0.7)
        for T in R["crossed"]:
            ax[0].plot((R["age_inj"] + T["t_cross"]) / Gyr, T["f_res"], "o", ms=3,
                       color=KIND_STYLE[T["kind"]]["color"])
    ax[0].set_xlabel("black-hole age [Gyr]")
    ax[0].set_ylabel(r"$f_{\rm orb}=\Omega/2\pi$ [Hz]")
    ax[0].set_title(f"Orbital tracks (dots: all fine/hyperfine crossings), μ={p.mu_eV:.1e} eV")
    kinds = ["hyperfine", "fine", "stalled", "none_ahead", "no_SR_state"]
    frac = [np.mean([R["outcome"] == k for R in results]) for k in kinds]
    ax[1].bar(kinds, frac, color=[OUTCOME_COL[k] for k in kinds])
    ax[1].set_ylabel("fraction")
    ax[1].set_title(f"First resonance ({len(results)} realizations)")
    ax[1].tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(fname, dpi=140)


def plot_cloud_track(track, fname):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
    t = track["t"] / Gyr
    ax[0].semilogy(t, track["M"] / Msun)
    ax[0].set_ylabel(r"$M\ [M_\odot]$")
    ax[1].plot(t, track["a"])
    ax[1].set_ylabel("spin a")
    for j, lv in enumerate(track["levels"]):
        ax[2].semilogy(t, np.maximum(track["Mc"][j] / track["M"], 1e-12), label=str(lv))
    ax[2].set_ylabel(r"$M_c/M$")
    ax[2].set_ylim(1e-12, 1)
    ax[2].legend()
    ax[2].set_xlabel("cosmic time [Gyr]")
    fig.tight_layout()
    fig.savefig(fname, dpi=140)


if __name__ == "__main__":
    import csv
    p = Params(mu_eV=1e-16, states=((2, 1, 1), (3, 2, 2), (3, 1, 1), (3, 0, 0), (4, 3, 3)),
               l_stars=(2, 3))
    res = run_monte_carlo(2000, p, seed=1)
    print("outcome fractions:")
    for k in ["hyperfine", "fine", "stalled", "none_ahead", "no_SR_state"]:
        print(f"  {k:16s} {np.mean([R['outcome'] == k for R in res]):.3f}")
    rows = crossing_table(res)
    with open("crossings.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} crossings written to crossings.csv")
    plot_baseline(res, p, "baseline_first_resonance.png")
    plot_gallery(res, p, "resonance_gallery.png", n_panels=6)

    # generalized demo: one history with superradiance on (free field)
    rng = np.random.default_rng(3)
    p2 = Params(logM_seed=(3, 3), logM_today=(6, 6))
    hist, inj, track = evolve_with_cloud(rng, p2, CloudModel(mu_eV=p2.mu_eV), return_track=True)
    plot_cloud_track(track, "cloud_history_demo.png")
