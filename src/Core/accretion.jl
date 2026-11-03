"""
    accretion.jl

Thin-disc baryonic accretion onto the BH, following arXiv:2510.19443 (Sec. 2.1)
and Bardeen (1970). Gas reaches the BH from the prograde ISCO, so

    dM/dt = Mdot_acc,      dJ/dt = (L_isco / E_isco) Mdot_acc,

with L_isco, E_isco the specific angular momentum and energy at the ISCO
(per unit rest mass). For a = J c / (G M^2),

    da/dt = (Mdot_acc / M) [ l(a) - 2a ],    l = c L_isco / (G M E_isco) .

Mdot_acc = f_Edd * Mdot_Edd(M_0) is held constant in time (M_0 = initial BH
mass), or (solve_system acc_dlnM_dt) proportional to M(t) for a fixed Eddington
ratio. Mdot_Edd = L_Edd / (eta c^2) = M / (eta t_Edd), t_Edd = sigma_T c / (4 pi G m_p).
The paper does not state eta; eta = 0.1 is the usual choice. As in the paper,
spin-up stops at the Thorne limit maxSpin = 0.998 (the mass keeps growing); it is
switched off smoothly over the last ACC_SPIN_TAPER below it. A hard switch at
maxSpin makes the ODE crawl once a cloud also spins the BH down: every step then
crosses the limit and is reset by the solver's spin callback.
"""

const ACC_SPIN_TAPER = 1e-3

# Eddington time sigma_T c / (4 pi G m_p) in years (cgs constants)
const T_EDD_YR = 6.6524587e-25 * 2.99792458e10 / (4π * 6.674e-8 * 1.67262192e-24) / YEAR_IN_SECONDS

"""Eddington accretion rate [M_sun/yr] for BH mass M [M_sun] and radiative efficiency eta."""
eddington_rate(M; eta=0.1) = M / (eta * T_EDD_YR)

"""Prograde ISCO radius in units of G M / c^2 (Bardeen, Press & Teukolsky 1972)."""
function kerr_isco(a)
    z1 = 1 + cbrt(1 - a^2) * (cbrt(1 + a) + cbrt(1 - a))
    z2 = sqrt(3a^2 + z1^2)
    return 3 + z2 - sqrt((3 - z1) * (3 + z1 + 2z2))
end

"""c L_isco / (G M E_isco): angular momentum per unit accreted mass-energy, in units of G M / c."""
function isco_l_over_e(a)
    x = kerr_isco(a)
    L = 2 / (3sqrt(3)) * (1 + 2sqrt(3x - 2))
    E = sqrt(1 - 2 / (3x))
    return L / E
end

"""
    accretion_rhs(M, a, Mdot_acc) -> (dM/dt, da/dt)

Accretion contribution to the BH mass [M_sun/yr] and spin [1/yr]. Spin-up is
switched off at a >= maxSpin.
"""
function accretion_rhs(M, a, Mdot_acc)
    a_c = clamp(a, 0.0, maxSpin)
    dadt = Mdot_acc / M * (isco_l_over_e(a_c) - 2a_c)
    if dadt > 0
        dadt *= clamp((maxSpin - a) / ACC_SPIN_TAPER, 0.0, 1.0)
    end
    return Mdot_acc, dadt
end
