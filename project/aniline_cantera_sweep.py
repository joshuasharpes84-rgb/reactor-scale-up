"""
Cantera cross-check of the full diameter sweep from aniline_pfr_feasibility.py -- does an
independent reactor-network solver reproduce the same critical-diameter finding (~55 mm)?
Same method as aniline_cantera_check.py (batch-reactor-as-PFR trick), same U(D) correlation
imported directly from the scipy model so only the reactor integration itself differs.
"""
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
import importlib.util

spec = importlib.util.spec_from_file_location("scipy_model", "aniline_pfr_feasibility.py")
scipy_model = importlib.util.module_from_spec(spec)
import sys
sys.argv = [sys.argv[0]]  # avoid the base module's own plt.show() blocking
import matplotlib
matplotlib.use('Agg')
spec.loader.exec_module(scipy_model)  # runs the scipy script once, gives us U_overall(), tau_base

L_over_D = scipy_model.L_over_D
tau = scipy_model.tau_base
T_feed = scipy_model.T_feed
Tc = scipy_model.Tc

D_values = np.linspace(0.005, 0.20, 25)  # coarser than the scipy sweep -- Cantera integration is slower
peak_T_ct = []
exit_X_ct = []

for D in D_values:
    L = L_over_D * D
    u = L / tau
    U, Re = scipy_model.U_overall(D, u)

    gas = ct.Solution('aniline_mech.yaml')
    gas.TPX = T_feed, ct.one_atm, {'CH3OH': 0.9802, 'C6H5NO2': 0.009901, 'H2': 0.009901}
    r = ct.ConstPressureReactor(gas, energy='on')
    r.volume = 1.0
    env = ct.Reservoir(ct.Solution('aniline_mech.yaml'))
    env.phase.TPX = Tc, ct.one_atm, {'CH3OH': 1.0}
    w = ct.Wall(r, env, A=4.0 / D, U=U)
    net = ct.ReactorNet([r])

    X0 = r.phase['C6H5NO2'].X[0]
    Tmax = T_feed
    n_steps = 150
    for t in np.linspace(0, tau, n_steps):
        net.advance(t)
        Tmax = max(Tmax, r.T)
    X_exit = r.phase['C6H5NO2'].X[0]
    peak_T_ct.append(Tmax)
    exit_X_ct.append(1 - X_exit / X0)
    print(f"D={D*1000:5.1f} mm  peakT={Tmax-273.15:6.1f} C  X_exit={100*(1-X_exit/X0):5.1f}%")

peak_T_ct = np.array(peak_T_ct)
exit_X_ct = np.array(exit_X_ct)
Tb = ct.Solution('aniline_mech.yaml').T  # placeholder, real boiling point set below
Tb_methanol_K = 337.8  # 64.65 C, from the scipy model's constants.Tbs[0] for methanol

runaway_idx = int(np.argmax(peak_T_ct > Tb_methanol_K)) if np.any(peak_T_ct > Tb_methanol_K) else None
print(f"\nMethanol boiling point: {Tb_methanol_K-273.15:.1f} C")
if runaway_idx:
    print(f"Cantera critical diameter ~ {D_values[runaway_idx]*1000:.1f} mm "
          f"(scipy model found ~54.6 mm)")
elif runaway_idx == 0:
    print("Cantera: runs away even at smallest D")
else:
    print("Cantera: no runaway found in swept range")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
axes[0].plot(D_values * 1000, peak_T_ct - 273.15, 'ro-', lw=2, label='Cantera')
axes[0].plot(scipy_model.D_values * 1000, scipy_model.peak_T - 273.15, 'r--', lw=1.5, alpha=0.6,
             label='scipy')
axes[0].axhline(Tb_methanol_K - 273.15, color='k', ls='--', label='MeOH boiling pt')
axes[0].set_xlabel('D (mm)'); axes[0].set_ylabel('Peak T (C)'); axes[0].legend()
axes[0].set_title('Peak T vs D: Cantera vs scipy')

axes[1].plot(D_values * 1000, exit_X_ct * 100, 'bo-', lw=2, label='Cantera')
axes[1].plot(scipy_model.D_values * 1000, scipy_model.exit_X * 100, 'b--', lw=1.5, alpha=0.6,
             label='scipy')
axes[1].set_xlabel('D (mm)'); axes[1].set_ylabel('Exit conversion (%)'); axes[1].legend()
axes[1].set_title('Conversion vs D: Cantera vs scipy')

plt.tight_layout()
plt.savefig('aniline_cantera_sweep_compare.png', dpi=130)
print("Saved comparison plot to aniline_cantera_sweep_compare.png")
