"""
Independent cross-check of aniline_pfr_feasibility.py using Cantera's own reactor-network solver,
instead of a hand-rolled scipy ODE.

Method: a steady 1D PFR is equivalent to a single fluid parcel evolving in time via t = z/u
(the standard Cantera PFR trick -- see Cantera's own PFR examples). Reactor jacket cooling is
reproduced via a Wall object with area/volume = 4/D (so U_wall*A_wall/V = 4*U/D, matching the
scipy model's (4U/D)*(T-Tc) cooling term exactly) connected to a fixed-temperature Reservoir.

Same reaction data as aniline_pfr_feasibility.py: A=1.9325 1/s, Ea=20604 J/mol (fitted from Prekob
et al. 2021 Table 2), dH_rxn=-552.9 kJ/mol (embedded via species formation enthalpies in
aniline_mech.yaml, sourced from the same NIST WebBook Hess's-Law derivation), same base-case
geometry (D=10mm, L/D=100), same coolant (15 C) and heat-transfer coefficient.
"""
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

# --- same base-case parameters as the scipy model ---
D = 0.01       # m
L_over_D = 100.0
tau = 4239.7   # s (from the scipy model's isothermal residence-time estimate)
T_feed = 313.15
Tc = 288.15
U = 68.0       # W/m2/K, from the scipy model's Dittus-Boelter+jacket calc at this D

L = L_over_D * D
u = L / tau

gas = ct.Solution('aniline_mech.yaml')
gas.TPX = T_feed, ct.one_atm, {'CH3OH': 0.9802, 'C6H5NO2': 0.009901, 'H2': 0.009901}

r = ct.ConstPressureReactor(gas, energy='on')
r.volume = 1.0  # m3 -- arbitrary; only the area/volume ratio below matters

env = ct.Reservoir(ct.Solution('aniline_mech.yaml'))
env.phase.TPX = Tc, ct.one_atm, {'CH3OH': 1.0}

wall_area = 4.0 / D  # m2, per m3 of reactor volume -- matches scipy model's (4U/D) term
w = ct.Wall(r, env, A=wall_area, U=U)

net = ct.ReactorNet([r])

t_end = tau
n_pts = 400
t_arr = np.linspace(0, t_end, n_pts)
T_arr = np.zeros(n_pts)
X_MNB_arr = np.zeros(n_pts)

X0 = r.phase['C6H5NO2'].X[0]

for i, t in enumerate(t_arr):
    net.advance(t)
    T_arr[i] = r.T
    X_MNB_arr[i] = r.phase['C6H5NO2'].X[0]

conversion = 1 - X_MNB_arr / X0
z_arr = u * t_arr

print(f"Cantera cross-check: D={D*1000:.1f} mm, L={L:.2f} m, u={u*100:.3f} cm/s, tau={tau:.0f} s")
print(f"Exit conversion = {conversion[-1]*100:.1f}%, peak T = {T_arr.max()-273.15:.1f} C "
      f"(feed {T_feed-273.15:.1f} C)")
print("scipy model gave: exit conversion = 79.3%, peak T = 40.0 C -- compare above")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
axes[0].plot(z_arr, T_arr - 273.15, 'r-', lw=2)
axes[0].set_xlabel('Axial position z (m)')
axes[0].set_ylabel('T (C)')
axes[0].set_title('Cantera reactor-network T(z)')

axes[1].plot(z_arr, conversion * 100, 'b-', lw=2)
axes[1].set_xlabel('Axial position z (m)')
axes[1].set_ylabel('Conversion (%)')
axes[1].set_title('Cantera reactor-network X(z)')

plt.tight_layout()
plt.savefig('aniline_cantera_check.png', dpi=130)
print("Saved plot to aniline_cantera_check.png")
