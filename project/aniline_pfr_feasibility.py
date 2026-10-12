"""
Feasibility model: jacketed PFR scale-up for MNB hydrogenation (DuPont/KBR aniline process).

Reaction: C6H5NO2(l) + 3 H2(g) -> C6H5NH2(l) + 2 H2O(l), pseudo-first-order in nitrobenzene.

Purpose: prove out the reactor-scaleup-tool physics on the first (flagship) reaction before
building the web tool or adding more reactions. Answers: as the tube diameter D grows (at fixed
L/D and fixed nominal residence time), when does the reactor hot spot run away?

Data sources (see "Reactor Scaleup Tool - Plan.md" for full citations):
  - dH_rxn = -552.9 kJ/mol: Hess's Law from NIST WebBook liquid-phase heats of formation.
  - k(T) = A*exp(-Ea/RT), A = 1.93 1/s, Ea = 20.6 kJ/mol: fitted here directly from the 4 real
    (T, k) pairs reported in Prekob et al., Nanomaterials 11(5):1172 (2021), Table 2, Pd/GCF
    catalyst (0.10 g catalyst / 200 mL of 0.25 mol/dm3 nitrobenzene in methanol, 20 bar H2,
    283-323 K). This is a LAB-SCALE, METHANOL-DILUTED system, not DuPont's neat/proprietary
    industrial catalyst -- used here as the best verified public kinetics for this reaction
    class. Running the model within/near 283-323 K (this script's operating point) avoids
    extrapolating the Arrhenius fit outside its validated range.
  - Liquid density/Cp for nitrobenzene, aniline, water: thermo library (Poling/DIPPR correlations).
  - Confirmed from the Meyers handbook: real DuPont/KBR reactor EXPORTS steam (0.15 kg/kg aniline)
    -- i.e. the real industrial reactor runs hot enough to raise steam (boiling/evaporative jacket)
    on a NEAT, concentrated feed. That's a different regime from the dilute lab kinetics used here,
    so this script instead models Prekob's actual system (dilute nitrobenzene in methanol) with a
    sensible-heat coolant jacket, consistent with their chiller-cooled autoclave setup -- an
    explicit, flagged simplification, not a claim about DuPont's real feed concentration or cooling
    scheme. (A first attempt at this model using a neat/undiluted nitrobenzene feed produced an
    immediate, nonphysical runaway even at the smallest lab diameter -- the exotherm from
    undiluted nitrobenzene at these rate constants vastly exceeds what any simple single-pass
    jacketed tube can remove. That is itself a real finding, not a bug: it is exactly why real
    industrial hydrogenations this exothermic are run diluted or with product recycle, never neat
    single-pass -- see the writeup accompanying this script.)

Model: ideal 1D PFR (no axial dispersion -- appropriate for a well-mixed liquid-phase tubular
reactor, unlike Zeta's rotary-kiln solids bed), constant properties evaluated at feed conditions,
single jacket-side coolant at fixed Tc.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from thermo import ChemicalConstantsPackage

# --- Reaction / kinetics (fitted from Prekob et al. 2021, Table 2, Pd/GCF) ---
A_arr = 1.9325       # 1/s
Ea = 20604.0         # J/mol
R = 8.314            # J/mol/K
dH_rxn = -552_900.0  # J/mol nitrobenzene (Hess's Law, NIST liquid-phase values)

def k_rate(T):
    return A_arr * np.exp(-Ea / (R * T))

# --- Bulk liquid = methanol solvent (Prekob's own system); nitrobenzene is the dilute solute ---
constants, correlations = ChemicalConstantsPackage.from_IDs(['methanol'])
T_feed = 313.15  # K (40 C) -- within the 283-323 K validated kinetics range
rho_feed = constants.MWs[0] / 1000 / correlations.VolumeLiquids[0](T_feed, 101325)  # kg/m3
Cp_feed_molar = correlations.HeatCapacityLiquids[0](T_feed)  # J/mol/K
Cp_feed_mass = Cp_feed_molar / (constants.MWs[0] / 1000)     # J/kg/K
C0 = 250.0  # mol/m3 nitrobenzene, matching Prekob et al.'s 0.25 mol/dm3

print(f"Feed T = {T_feed} K, solvent(methanol) rho = {rho_feed:.1f} kg/m3, "
      f"Cp = {Cp_feed_mass:.1f} J/kg/K, C0(nitrobenzene) = {C0:.1f} mol/m3")

# --- Heat transfer: Dittus-Boelter process-side h + assumed jacket-side h, combine to U ---
mu_feed = correlations.ViscosityLiquids[0](T_feed, 101325)          # Pa.s
kcond_feed = correlations.ThermalConductivityLiquids[0](T_feed, 101325)  # W/m/K
# Boiling/evaporative jacket (motivated by the confirmed real DuPont/KBR steam export, Table
# 2.1.1): h from ChemE Rules of Thumb Table 1 "Film Resistances", Evaporation/Water row,
# R=0.0007 hr.ft2.F/Btu -> h=1/R=1428.6 Btu/hr/ft2/F = 8112 W/m2K. Note: at this reactor's actual
# ~40-100 C operating range, a jacket boiling at Tc=15 C physically implies a vacuum-operated
# water loop or a refrigerant, not atmospheric steam generation like DuPont's real (much hotter)
# process -- using water's evaporative film coefficient as the representative boiling-regime
# number, not a claim about the literal jacket fluid.
h_jacket = 8112.0  # W/m2/K

def U_overall(D, u):
    """Nu(Re) blends smoothly across the laminar-turbulent transition (Re 2300-4000) instead of
    switching discontinuously -- a hard if/else here produced sharp, unphysical jumps in peak T
    vs D in the recycle-loop sweep (traced and confirmed as a correlation artifact, not real
    reactor bistability, via forward/backward continuation -- see aniline_pfr_recycle.py). Real
    transitional flow doesn't jump; the correlation shouldn't either."""
    Re = rho_feed * u * D / mu_feed
    Pr = Cp_feed_mass * mu_feed / kcond_feed

    Nu_lam = 3.66  # laminar, constant wall T

    # Gnielinski (turbulent, and smoother/more accurate near the transition than Dittus-Boelter)
    f = (0.790 * np.log(max(Re, 1.0)) - 1.64) ** -2
    Nu_turb = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
    Nu_turb = max(Nu_turb, Nu_lam)  # guard outside Gnielinski's own valid range

    Re_lo, Re_hi = 2300.0, 4000.0  # standard transitional-regime bracket
    if Re <= Re_lo:
        Nu = Nu_lam
    elif Re >= Re_hi:
        Nu = Nu_turb
    else:
        x = (Re - Re_lo) / (Re_hi - Re_lo)
        w = 3 * x**2 - 2 * x**3  # cubic smoothstep: w(0)=0, w(1)=1, zero slope at both ends
        Nu = (1 - w) * Nu_lam + w * Nu_turb

    h_process = Nu * kcond_feed / D
    return 1.0 / (1.0 / h_process + 1.0 / h_jacket), Re

# --- Coolant ---
Tc = 288.15  # K (15 C), sensible-heat jacket, consistent with a chilled/cooling-water loop

# --- Base-case sizing: pick target conversion, get residence time at feed T ---
X_target = 0.95
k_feed = k_rate(T_feed)
tau_base = -np.log(1 - X_target) / k_feed
print(f"k(T_feed) = {k_feed:.3e} 1/s, tau (isothermal estimate) = {tau_base:.1f} s "
      f"({tau_base/60:.1f} min)")


def integrate_pfr(D, L_over_D, tau):
    """Integrate the 1D PFR energy+mass balance for a given diameter, L/D, and residence time."""
    L = L_over_D * D
    u = L / tau
    U, Re = U_overall(D, u)

    def odes(z, y):
        C, T = y
        k = k_rate(T)
        rate = k * C  # mol/m3/s
        dCdz = -rate / u
        gen = -dH_rxn * rate  # W/m3 (positive, exothermic)
        remove = (4 * U / D) * (T - Tc)  # W/m3
        dTdz = (gen - remove) / (rho_feed * Cp_feed_mass * u)
        return [dCdz, dTdz]

    sol = solve_ivp(odes, [0, L], [C0, T_feed], max_step=L / 500, dense_output=False)
    return sol, u, U, Re, L


# --- Sanity check: base-case lab diameter ---
D0 = 0.01  # m, 1 cm ID -- lab/bench scale
L_over_D = 100.0
sol0, u0, U0, Re0, L0 = integrate_pfr(D0, L_over_D, tau_base)
X_exit = 1 - sol0.y[0, -1] / C0
T_peak = sol0.y[1].max()
print(f"\nBase case D={D0*1000:.1f} mm, L={L0:.2f} m, u={u0*100:.2f} cm/s, Re={Re0:.0f}, "
      f"U={U0:.0f} W/m2K")
print(f"Exit conversion = {X_exit*100:.1f}%, peak T = {T_peak-273.15:.1f} C "
      f"(feed {T_feed-273.15:.1f} C, coolant {Tc-273.15:.1f} C)")

# --- Scale-up sweep: hold L/D and tau fixed, increase D, track peak T and mass throughput ---
D_values = np.linspace(0.005, 0.20, 60)  # 5 mm to 20 cm
peak_T = []
exit_X = []
mdot = []  # kg/s product-equivalent throughput (based on inlet MNB mass flow)
Re_list = []
for D in D_values:
    sol, u, U, Re, L = integrate_pfr(D, L_over_D, tau_base)
    peak_T.append(sol.y[1].max())
    exit_X.append(1 - sol.y[0, -1] / C0)
    area = np.pi * D**2 / 4
    Q = u * area  # m3/s
    mdot.append(Q * rho_feed)
    Re_list.append(Re)

peak_T = np.array(peak_T)
exit_X = np.array(exit_X)
mdot = np.array(mdot)

# Runaway criterion: peak T exceeds the bulk solvent's (methanol) normal boiling point -- physically
# the reactor would start flashing/boiling, a hard physical limit, not just a soft safety margin.
Tb_MNB = constants.Tbs[0]
runaway_idx = int(np.argmax(peak_T > Tb_MNB)) if np.any(peak_T > Tb_MNB) else None

print(f"\nMethanol normal boiling point: {Tb_MNB-273.15:.1f} C")
if runaway_idx is not None and runaway_idx > 0:
    D_crit = D_values[runaway_idx]
    mdot_crit = mdot[runaway_idx]
    print(f"Peak T crosses boiling point between D={D_values[runaway_idx-1]*1000:.1f} mm and "
          f"{D_values[runaway_idx]*1000:.1f} mm")
    print(f"Critical diameter ~ {D_crit*1000:.1f} mm -> max safe MNB throughput "
          f"~ {mdot_crit*3600:.2f} kg/h (at L/D={L_over_D:.0f}, tau={tau_base:.0f} s)")
elif runaway_idx == 0:
    print("Runs away even at the smallest swept diameter -- base case is already infeasible at "
          "this L/D and residence time; need a smaller D, shorter tau, or better cooling.")
else:
    print("No runaway found in the swept diameter range -- cooling keeps up throughout.")

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

axes[0].plot(D_values * 1000, peak_T - 273.15, 'r-', lw=2)
axes[0].axhline(Tb_MNB - 273.15, color='k', ls='--', label='MNB boiling point')
axes[0].set_xlabel('Tube diameter D (mm)')
axes[0].set_ylabel('Peak reactor T (C)')
axes[0].set_title('Hot-spot growth with scale-up')
axes[0].legend()

axes[1].plot(D_values * 1000, exit_X * 100, 'b-', lw=2)
axes[1].set_xlabel('Tube diameter D (mm)')
axes[1].set_ylabel('Exit conversion (%)')
axes[1].set_title('Conversion vs. diameter')

axes[2].plot(D_values * 1000, mdot * 3600, 'g-', lw=2)
if runaway_idx:  # note: intentionally excludes idx 0 (immediate-runaway case) from the vline
    axes[2].axvline(D_values[runaway_idx] * 1000, color='k', ls='--', label='runaway limit')
    axes[2].legend()
axes[2].set_xlabel('Tube diameter D (mm)')
axes[2].set_ylabel('MNB throughput (kg/h)')
axes[2].set_title('Max safe scale-up throughput')

plt.tight_layout()
plt.savefig('aniline_pfr_feasibility.png', dpi=130)
print("\nSaved plot to aniline_pfr_feasibility.png")
plt.show()
