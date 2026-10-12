"""
Feasibility model: packed-bed PFR scale-up for Cumene (UOP Q-Max benzene alkylation).

Reaction: C6H6(l) + C3H6(dissolved) -> C6H5CH(CH3)2(l), pseudo-first-order in propylene, plus a
DIPB transalkylation side reaction (kinetics sourced but not yet modeled -- see Plan.md).

Scoping finding (2026-08-17, confirmed by direct numerical check before further build effort per
user instruction): this reaction's exotherm (-110.6 kJ/mol) is much milder than Aniline's (-552.9
kJ/mol), and Leva's packed-bed correlation gives U in the thousands of W/m2K (vs Aniline's tens to
hundreds) -- so this reactor is NOT realistically thermal-runaway-limited at practical diameters;
it stays well-cooled throughout. What DOES become prohibitive is pressure drop (Ergun), which rises
by >1000x over the same diameter range that leaves peak T only mildly elevated. So unlike Aniline
(single thermal criterion), Cumene's recommended diameter is governed by TWO criteria of different
character: a thermal Semenov/F-K stability bound (same closed-form method as Aniline, for
completeness/consistency) and a pressure-drop ceiling (the one that actually binds in practice).
User decision 2026-08-17: show both, recommend whichever is smaller (binds first) -- same pattern
as Aniline's dual-criteria approach, generalized.

Data sources (see "Reactor Scaleup Tool - Plan.md", "Reaction Data" section, confirmed 2026-08-14):
  - dH_rxn = -110.6 kJ/mol (propylene basis): Hess's Law from NIST WebBook liquid-phase heats of
    formation (benzene, cumene, propylene).
  - Kinetics: rate = K1*C_propylene (pseudo-first-order in propylene), K1 = 6510*exp(-52564/RT),
    Ea=52.6 kJ/mol, A=6510 1/s (E given in kJ/kmol by the source, numerically identical to J/mol).
    Source: Ben Hamida, Preprints.org 202601.1995.v1 (2026), an Aspen HYSYS PFR simulation study
    of the real Q-Max process (170-180 C, 3540-3600 kPa, beta-zeolite Si/Al=25, 2mm pellets,
    liquid-phase fixed-bed, benzene:propylene ~0.75-0.95:1).
  - MODELING DECISION (not in the source, stated explicitly): K1 is treated as an already
    bed-effective (lumped) pseudo-homogeneous rate constant, since it was fitted/used directly in
    a reactor-scale PFR simulation, not measured as intrinsic single-pellet kinetics -- no separate
    catalyst effectiveness factor (Thiele modulus) is layered on top, to avoid double-counting
    whatever internal mass-transfer limitation is already implicit in the fitted correlation.
    Flagged as an open item if a more rigorous intrinsic-kinetics treatment is wanted later.
  - Benzene:propylene ratio = 0.85 (mid of the source's stated 0.75-0.95:1) taken at face value.
  - UNSOURCED, reasonable-engineering-judgment assumptions (same role as Aniline's Tc=15C jacket
    choice): target conversion 90%, coolant Tc = T_feed - 20K, allowable pressure drop = 10% of
    operating pressure (a commonly-cited packed-bed design heuristic -- not verified here against
    a specific page in this project's own library; flag if it needs firming up).

Heat transfer: Leva's correlation for wall-to-bed heat transfer in a packed tube (widely cited,
e.g. Coulson & Richardson Vol 2 particulate-solids chapters), Nu_w = h_w*dp/kf = 2.03*Re_p^0.8 *
exp(-6*dp/D), valid for Re_p > ~20. Re_p = rho*u_superficial*dp/mu (particle-diameter-based).

Pressure drop: Ergun equation, superficial-velocity basis.

Thermal stability: same exact Semenov/Frank-Kamenetskii closed-form as aniline_semenov_runaway.py
(R*T*^2 - Ea*T* + Ea*Tc = 0), reusing this project's own validated methodology.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from thermo import ChemicalConstantsPackage

# --- Reaction / kinetics (Ben Hamida 2026, Aspen HYSYS PFR sim study fit) ---
A_arr = 6510.0        # 1/s
Ea = 52564.0          # J/mol (source's kJ/kmol == J/mol numerically)
R = 8.314             # J/mol/K
dH_rxn = -110_600.0   # J/mol propylene (Hess's Law, NIST liquid-phase values)

def k_rate(T):
    return A_arr * np.exp(-Ea / (R * T))

# --- Bulk liquid = benzene solvent; propylene is the dissolved, limiting reactant ---
constants, correlations = ChemicalConstantsPackage.from_IDs(['benzene'])
T_feed = 448.15   # K (175 C), midpoint of Q-Max's stated 170-180 C range
P_op = 3_570_000  # Pa (~3570 kPa), midpoint of stated 3540-3600 kPa
rho_feed = constants.MWs[0] / 1000 / correlations.VolumeLiquids[0](T_feed, P_op)  # kg/m3
Cp_feed_molar = correlations.HeatCapacityLiquids[0](T_feed)  # J/mol/K
Cp_feed_mass = Cp_feed_molar / (constants.MWs[0] / 1000)     # J/kg/K
mu_feed = correlations.ViscosityLiquids[0](T_feed, P_op)          # Pa.s
kcond_feed = correlations.ThermalConductivityLiquids[0](T_feed, P_op)  # W/m/K

BC_ratio = 0.85  # benzene:propylene molar ratio, mid of source's stated 0.75-0.95:1
x_propylene = 1.0 / (1.0 + BC_ratio)
C_total = rho_feed / (constants.MWs[0] / 1000)  # mol/m3, approximated with pure benzene's molar density
C0 = x_propylene * C_total  # mol/m3 propylene

# --- Packed bed geometry (Q-Max: 2mm beta-zeolite pellets, Si/Al=25) ---
dp = 0.002       # m
eps_bed = 0.40   # typical packed-bed voidage, unsourced standard assumption
h_jacket = 8112.0  # W/m2K -- reused from Aniline's boiling-jacket film coefficient (same
                    # representative-boiling-regime simplification, not literal-fluid claim)

def U_overall(D, u_s):
    """Leva's packed-bed wall heat-transfer correlation. Returns (U_overall, Re_p).

    Leva's correlation is only valid for Re_p > ~20 (Coulson & Richardson Vol 2). Confirmed by
    direct check (2026-08-17): extrapolating the Re_p^0.8 power law below that floor produces a
    sharp, unphysical jump (e.g. a spurious ~400C "runaway" at D=10mm, L/D=25, Re_p=4.98, that
    vanishes the instant Re_p crosses 20 by D=50mm) -- a correlation-validity artifact, not real
    reactor physics (real low-Re packed beds still have meaningful particle-contact/conduction
    heat transfer, not captured by extrapolating a turbulent-leaning power law toward zero). Fixed
    by flooring the value used IN the correlation at Re_p=20 -- the correlation's own boundary,
    not an extrapolation past it -- while still returning the true Re_p for diagnostics/display.
    """
    Re_p = rho_feed * u_s * dp / mu_feed
    Re_p_eff = max(Re_p, 20.0)
    Nu_w = 2.03 * (Re_p_eff ** 0.8) * np.exp(-6 * dp / D)
    h_process = Nu_w * kcond_feed / dp
    return 1.0 / (1.0 / h_process + 1.0 / h_jacket), Re_p

def ergun_dPdz(u_s):
    """Ergun equation, Pa/m (superficial velocity basis)."""
    term1 = 150 * mu_feed * u_s * (1 - eps_bed) ** 2 / (dp ** 2 * eps_bed ** 3)
    term2 = 1.75 * rho_feed * u_s ** 2 * (1 - eps_bed) / (dp * eps_bed ** 3)
    return term1 + term2

# --- Coolant + base-case sizing ---
Tc = T_feed - 20.0   # K -- unsourced, reasonable-choice jacket coolant (analogous to Aniline's Tc)
X_target = 0.90      # unsourced, reasonable per-pass conversion target
k_feed = k_rate(T_feed)
tau_base = -np.log(1 - X_target) / k_feed

# --- Pressure-drop ceiling: 10% of operating pressure, common packed-bed design heuristic ---
dP_allowable = 0.10 * P_op  # Pa


def integrate_pfr(D, L_over_D, tau):
    L = L_over_D * D
    u_s = L / tau
    U, Re_p = U_overall(D, u_s)

    def odes(z, y):
        C, T = y
        k = k_rate(T)
        rate = k * C
        dCdz = -rate / u_s
        gen = -dH_rxn * rate
        remove = (4 * U / D) * (T - Tc)
        dTdz = (gen - remove) / (rho_feed * Cp_feed_mass * u_s)
        return [dCdz, dTdz]

    sol = solve_ivp(odes, [0, L], [C0, T_feed], max_step=L / 500, dense_output=False)
    dP = ergun_dPdz(u_s) * L
    return sol, u_s, U, Re_p, L, dP


# --- Exact Semenov/FK critical hot-spot temperature (same closed form as Aniline) ---
disc = Ea**2 - 4 * R * Ea * Tc
T_star = (Ea - np.sqrt(disc)) / (2 * R)


def k1_crit_at(C):
    g = -dH_rxn * A_arr * np.exp(-Ea / (R * T_star)) * C
    return g * Ea / (R * T_star**2)


if __name__ == "__main__":
    print(f"Feed T = {T_feed} K ({T_feed-273.15:.0f} C), P = {P_op/1000:.0f} kPa")
    print(f"Benzene(l): rho={rho_feed:.1f} kg/m3, Cp={Cp_feed_mass:.1f} J/kg/K, mu={mu_feed*1e3:.3f} mPa.s, "
          f"k={kcond_feed:.3f} W/m/K")
    print(f"x_propylene = {x_propylene:.3f}, C0(propylene) = {C0:.1f} mol/m3")
    print(f"k(T_feed) = {k_feed:.3e} 1/s, tau (isothermal estimate) = {tau_base:.2f} s ({tau_base/60:.2f} min)")
    print(f"Semenov T* = {T_star-273.15:.1f} C, dT_crit = {T_star-Tc:.1f} K above Tc={Tc-273.15:.1f} C")
    print(f"Allowable dP (10% of {P_op/1000:.0f} kPa) = {dP_allowable/1000:.1f} kPa\n")

    # D_MIN = 50mm: enforces D/dp >= 25 (dp=2mm), a standard packed-bed design rule of thumb
    # (commonly D/dp >= 20) to avoid wall-channeling/flow-maldistribution -- confirmed 2026-08-17
    # this also matters numerically: Leva's own exp(-6*dp/D) wall-proximity term correctly predicts
    # a severe heat-transfer penalty as D approaches dp, which combined with the Re_p<20 floor
    # produced spurious small-D "instability" at D<20mm that isn't a design anyone would build.
    D_values = np.linspace(0.05, 0.30, 60)
    LOD_test = [25, 50, 100, 150, 200]

    print("=== D_crit by criterion, across L/D ===")
    for LOD in LOD_test:
        peakT, dPs, k1_actual, k1_crit = [], [], [], []
        for D in D_values:
            sol, u_s, U, Re_p, L, dP = integrate_pfr(D, LOD, tau_base)
            peakT.append(sol.y[1].max())
            dPs.append(dP)
            k1_actual.append(4 * U / D)
            k1_crit.append(k1_crit_at(C0))  # base case, no recycle -> C_in = C0
        peakT = np.array(peakT)
        dPs = np.array(dPs)
        margin = np.array(k1_actual) - np.array(k1_crit)

        idx_p = int(np.argmax(dPs > dP_allowable)) if np.any(dPs > dP_allowable) else None
        Dp = D_values[idx_p] if idx_p else None
        idx_s = int(np.argmax(margin < 0)) if np.any(margin < 0) else None
        Ds = D_values[idx_s] if idx_s else None

        print(f"L/D={LOD:4d}: D_crit_pressure = {f'{Dp*1000:.1f} mm' if Dp else 'not reached (300mm)'}"
              f"   D_crit_semenov = {f'{Ds*1000:.1f} mm' if Ds else 'not reached (300mm)'}")
