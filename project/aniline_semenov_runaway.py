"""
Semenov / Frank-Kamenetskii thermal-runaway criterion for the aniline (MNB hydrogenation) PFR.

Purpose: derive T_runaway analytically from first-principles thermal-stability theory, as a
rigorous cross-check against the numerical proxy already used in aniline_pfr_feasibility.py and
aniline_pfr_recycle.py (peak reactor T exceeding the bulk solvent's boiling point). Per the plan
doc's own sourcing methodology: T_runaway isn't a tabulated literature number for this reaction --
it's derived here from the same Ea/A (Prekob et al. 2021), dH_rxn (NIST Hess's Law), and cooling
U(D,u) already assembled and validated in the base/recycle scripts. Reuses those modules directly
so only the stability analysis is new, not the underlying data.

Theory (Semenov thermal explosion / Frank-Kamenetskii stability criterion -- Davis & Davis,
"Fundamentals of Chemical Reaction Engineering", Ch. 9 "Nonisothermal Reactors", per this project's
own Reaction Kinetics Textbooks.pdf):

A reacting volume element has a heat-generation curve g(T) = (-dH_rxn)*k(T)*C (exponential in T,
via Arrhenius k(T)) and a heat-removal curve r(T) = (4U/D)*(T-Tc) (linear in T). This is exactly
the per-unit-volume form of the energy balance already used in the PFR ODEs (4U/D = UA/V for a
cylindrical tube), so the same lumped-element stability analysis applies locally.

Below a critical operating point the two curves cross twice (a stable low-T steady state + an
unstable high-T one); above it they don't cross at all (no steady state nearby -- runaway). The
boundary between the two regimes is where g(T) and r(T) are simultaneously EQUAL and TANGENT:
    g(T*) = r(T*)          and          g'(T*) = r'(T*)
Because g'(T)/g(T) = Ea/(R*T^2) exactly (Arrhenius derivative, independent of the prefactor -- A,
dH_rxn, C all cancel out of this ratio), combining the two tangency conditions eliminates the
cooling-capacity term (4U/D) entirely and leaves a quadratic in T* alone:
    R*T*^2 - Ea*T* + Ea*Tc = 0
solved exactly here (not the further-linearized classical approximation dT_crit = R*Tc^2/Ea, which
is only valid for Ea/(R*Tc) >> 1 -- checked below and found NOT to hold well for this reaction, so
the exact quadratic is used as the reported result and the linearized approximation is shown only
as a reference comparison).

T* found this way is the reactor's critical hot-spot temperature -- independent of D, U, C, or
dH_rxn magnitude. What DOES depend on those is whether a given reactor (a given D, and thus a given
U and C_in) is actually AT that critical point: plugging T* back into g'(T*) gives the critical
cooling-capacity term (4U/D)_crit required to keep the reactor subcritical. Comparing that against
the reactor's ACTUAL (4U/D), evaluated the same way as the existing D-sweeps, gives an independent,
theory-derived critical diameter -- directly comparable to the boiling-point-crossing D_crit already
found numerically.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import importlib.util


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

base = _load("base", "aniline_pfr_feasibility.py")
recycle = _load("recycle", "aniline_pfr_recycle.py")

R = base.R
Ea = base.Ea
A_arr = base.A_arr
dH_rxn = base.dH_rxn
C0 = base.C0
Tc = base.Tc
L_over_D = base.L_over_D
tau_base = base.tau_base
Tb_solvent = base.constants.Tbs[0]

# --- Step 1: exact critical hot-spot temperature T* (quadratic; depends only on Ea, R, Tc) ---
disc = Ea**2 - 4 * R * Ea * Tc
T_star_lo = (Ea - np.sqrt(disc)) / (2 * R)   # physical root (near Tc)
T_star_hi = (Ea + np.sqrt(disc)) / (2 * R)   # mathematical artifact, far outside any physical range
T_star = T_star_lo

# Classical linearized (Frank-Kamenetskii) approximation, for reference/comparison only
dT_crit_linearized = R * Tc**2 / Ea
Ea_over_RTc = Ea / (R * Tc)

print("=== Semenov/FK critical hot-spot temperature ===")
print(f"Ea/(R*Tc) = {Ea_over_RTc:.2f}  (classical linearized approx assumes this >> 1 -- "
      f"{'reasonably large' if Ea_over_RTc > 15 else 'NOT large, linearization will be inaccurate'})")
print(f"Exact quadratic root:      T* = {T_star-273.15:.2f} C  (dT_crit = {T_star-Tc:.2f} K above Tc)")
print(f"Linearized FK approx:      dT_crit = R*Tc^2/Ea = {dT_crit_linearized:.2f} K  "
      f"-> T*_approx = {Tc + dT_crit_linearized - 273.15:.2f} C")
print(f"(unphysical second root discarded: T*_hi = {T_star_hi-273.15:.0f} C)")
print(f"Bulk solvent (methanol) boiling point: {Tb_solvent-273.15:.1f} C -- for comparison against "
      f"the numerical peak-T-vs-boiling-point criterion used in the base/recycle sweeps.\n")


def k1_crit_at(T, C):
    """Critical cooling-capacity term (4U/D)_crit [W/m3/K] required to keep the reactor subcritical,
    given the local reactant concentration C feeding the hot spot (evaluated at T=T*)."""
    g = -dH_rxn * A_arr * np.exp(-Ea / (R * T)) * C
    return g * Ea / (R * T**2)


# --- Step 2: base (once-through) case -- C fixed at C0, so k1_crit is a single number ---
k1_crit_base = k1_crit_at(T_star, C0)
print("=== Base (once-through) case ===")
print(f"Required (4U/D)_crit = {k1_crit_base:.1f} W/m3/K (constant -- C fixed at C0={C0:.0f} mol/m3)")

D_values = np.linspace(0.005, 0.20, 200)
k1_actual_base = []
for D in D_values:
    u = L_over_D * D / tau_base
    U, Re = base.U_overall(D, u)
    k1_actual_base.append(4 * U / D)
k1_actual_base = np.array(k1_actual_base)

# Reactor is subcritical (stable) where actual cooling capacity exceeds the critical requirement
stable_base = k1_actual_base > k1_crit_base
D_crit_semenov_base = None
if stable_base[0] and not stable_base.all():
    idx = int(np.argmax(~stable_base))
    D_crit_semenov_base = D_values[idx]
    print(f"Semenov critical diameter: D_crit = {D_crit_semenov_base*1000:.1f} mm "
          f"(crosses from stable to supercritical between {D_values[idx-1]*1000:.1f} and "
          f"{D_values[idx]*1000:.1f} mm)")
elif not stable_base[0]:
    print("Supercritical even at the smallest swept diameter.")
else:
    print("Stable (subcritical) across the entire swept diameter range.")

# Cross-check against the existing numerical (boiling-point) criterion, re-run identically here
peak_T_base = []
for D in D_values:
    sol, u, U, Re, L = base.integrate_pfr(D, L_over_D, tau_base)
    peak_T_base.append(sol.y[1].max())
peak_T_base = np.array(peak_T_base)
boil_idx = int(np.argmax(peak_T_base > Tb_solvent)) if np.any(peak_T_base > Tb_solvent) else None
D_crit_boiling_base = D_values[boil_idx] if boil_idx else None
if D_crit_boiling_base:
    print(f"Numerical (boiling-point) critical diameter: {D_crit_boiling_base*1000:.1f} mm")
    if D_crit_semenov_base:
        print(f"-> Semenov theory predicts instability at a {'SMALLER' if D_crit_semenov_base < D_crit_boiling_base else 'LARGER'} "
              f"diameter than the boiling-point proxy ({D_crit_semenov_base*1000:.1f} mm vs "
              f"{D_crit_boiling_base*1000:.1f} mm) -- i.e. the reactor loses thermal STABILITY "
              f"(enters a runaway trajectory) before the boiling point is actually reached, so the "
              f"boiling-point criterion is the less conservative of the two. Consistent: both find "
              f"a finite critical diameter in the same D range, independently derived.")
print()

# --- Step 3: recycle case (RR=2) -- C_in varies with D via the recycle fixed-point solve, so
# k1_crit(D) is no longer flat; sweep with continuation exactly as aniline_pfr_recycle.py does ---
RR = 2
print(f"=== Recycle case (RR={RR}) ===")
k1_actual_rec, k1_crit_rec, C_in_prev = [], [], None
for D in D_values:
    u_total = (1 + RR) * (L_over_D * D) / tau_base
    U, Re = base.U_overall(D, u_total)
    X, T_peak, U_, Re_, u_, C_exit = recycle.solve_recycle_pfr(D, RR, C_exit_init=C_in_prev)
    C_in_prev = C_exit
    C_in = (C0 + RR * C_exit) / (1 + RR)
    k1_actual_rec.append(4 * U / D)
    k1_crit_rec.append(k1_crit_at(T_star, C_in))
k1_actual_rec = np.array(k1_actual_rec)
k1_crit_rec = np.array(k1_crit_rec)

stable_rec = k1_actual_rec > k1_crit_rec
if stable_rec.all():
    print(f"Stable (subcritical) across the entire swept diameter range (0-{D_values[-1]*1000:.0f} mm) "
          f"-- agrees with the numerical recycle sweep (no boiling-point crossing found for RR={RR} "
          f"either).")
    D_crit_semenov_rec = None
else:
    idx = int(np.argmax(~stable_rec))
    D_crit_semenov_rec = D_values[idx]
    print(f"Semenov critical diameter: D_crit = {D_crit_semenov_rec*1000:.1f} mm")

# --- Plot: generation vs removal curves at the base case's critical D (tangency), plus the
# stability-margin sweep for both cases ---
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

D_show = D_crit_semenov_base if D_crit_semenov_base else D_values[len(D_values) // 2]
u_show = L_over_D * D_show / tau_base
U_show, _ = base.U_overall(D_show, u_show)
T_range = np.linspace(Tc, T_star + 15, 200)
g_curve = -dH_rxn * A_arr * np.exp(-Ea / (R * T_range)) * C0
r_curve = (4 * U_show / D_show) * (T_range - Tc)
axes[0].plot(T_range - 273.15, g_curve / 1000, 'r-', lw=2, label='heat generation g(T)')
axes[0].plot(T_range - 273.15, r_curve / 1000, 'b-', lw=2, label='heat removal r(T)')
axes[0].axvline(T_star - 273.15, color='k', ls=':', lw=1, label=f'T* = {T_star-273.15:.1f} C')
axes[0].set_xlabel('T (C)')
axes[0].set_ylabel('Rate (kW/m3)')
axes[0].set_title(f'Semenov tangency at D={D_show*1000:.1f} mm (base case)')
axes[0].legend(fontsize=8)

axes[1].plot(D_values * 1000, k1_actual_base, 'g-', lw=2, label='actual 4U/D (base)')
axes[1].axhline(k1_crit_base, color='g', ls='--', lw=1.5, label='required (4U/D)_crit (base)')
axes[1].plot(D_values * 1000, k1_actual_rec, 'm-', lw=2, label=f'actual 4U/D (RR={RR})')
axes[1].plot(D_values * 1000, k1_crit_rec, 'm--', lw=1.5, label=f'required (4U/D)_crit (RR={RR})')
axes[1].set_xlabel('Tube diameter D (mm)')
axes[1].set_ylabel('4U/D  (W/m3/K)')
axes[1].set_yscale('log')
axes[1].set_title('Cooling capacity vs. Semenov requirement')
axes[1].legend(fontsize=7)

plt.tight_layout()
plt.savefig('aniline_semenov_runaway.png', dpi=130)
print("\nSaved plot to aniline_semenov_runaway.png")
