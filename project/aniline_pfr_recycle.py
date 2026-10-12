"""
Recycle-loop extension of aniline_pfr_feasibility.py: same physical reactor (D, L/D), but reactor
effluent is split into product + recycle, recycle passes through an (idealized, perfect) external
cooler back to T_feed, remixes with fresh feed, and re-enters the reactor. Tests the user's own
diagnosis: recycle boosts velocity/Re (fixing the laminar-flow process-side bottleneck found in
the base model) while also diluting reactor-inlet concentration and cutting per-pass residence
time (both reducing local heat generation per unit volume) -- attacking the runaway problem from
both the generation side and the removal side at once, unlike the jacket-only change tried before.

Topology: fresh feed Q0 (at C_A0, T_feed) + recycle Q_R (at C_A,exit, cooled back to T_feed) mix
-> reactor (D, L, same as the base model) -> exit stream splits into product (= Q0, steady state)
and recycle (= Q_R = RR*Q0). Reactor inlet concentration C_in = (C_A0 + RR*C_exit)/(1+RR) is
implicit (depends on the very C_exit it helps determine) -- solved by fixed-point iteration.

Reuses the base model's kinetics, thermo, dH_rxn, and U(D,u) correlation directly (imports the
module) so the only thing that differs is the recycle topology.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import importlib.util

spec = importlib.util.spec_from_file_location("base", "aniline_pfr_feasibility.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)  # runs the base script once; gives k_rate, U_overall, constants, etc.

C_A0 = base.C0
T_feed = base.T_feed
Tc = base.Tc
dH_rxn = base.dH_rxn
rho_feed = base.rho_feed
Cp_feed_mass = base.Cp_feed_mass
L_over_D = base.L_over_D
tau_oncethrough = base.tau_base  # per-pass residence time basis at RR=0 (same physical reactor size)
Tb_solvent = base.constants.Tbs[0]


def integrate_pass(D, u, C_in, U):
    """Integrate one reactor pass given inlet concentration/temperature and velocity."""
    L = L_over_D * D

    def odes(z, y):
        C, T = y
        k = base.k_rate(T)
        rate = k * C
        dCdz = -rate / u
        gen = -dH_rxn * rate
        remove = (4 * U / D) * (T - Tc)
        dTdz = (gen - remove) / (rho_feed * Cp_feed_mass * u)
        return [dCdz, dTdz]

    sol = solve_ivp(odes, [0, L], [C_in, T_feed], max_step=L / 500)
    return sol.y[0, -1], sol.y[1].max(), sol.y[1, -1]


def solve_recycle_pfr(D, RR, C_exit_init=None, tol=1e-8, max_iter=300, damping=0.5):
    """Damped fixed-point iteration (successive substitution damped 50% -- naive substitution was
    prone to oscillating rather than converging near the steep ignition/extinction region)."""
    u_total = (1 + RR) * (L_over_D * D) / tau_oncethrough
    U, Re = base.U_overall(D, u_total)

    C_exit_guess = C_A0 * 0.5 if C_exit_init is None else C_exit_init
    for _ in range(max_iter):
        C_in = (C_A0 + RR * C_exit_guess) / (1 + RR)
        C_exit_new, T_peak, T_exit = integrate_pass(D, u_total, C_in, U)
        step = C_exit_new - C_exit_guess
        if abs(step) < tol * C_A0:
            C_exit_guess = C_exit_new
            break
        C_exit_guess = C_exit_guess + damping * step
    else:
        print(f"    WARNING: fixed-point iteration did not converge at D={D*1000:.1f}mm, RR={RR}")
    X_overall = 1 - C_exit_guess / C_A0
    return X_overall, T_peak, U, Re, u_total, C_exit_guess


def continuation_sweep(D_values, RR, reverse=False):
    """Sweep D using the PREVIOUS D's converged C_exit as the next starting guess, to trace a
    single solution branch rather than risk landing on different branches at different D (which
    would show up as spurious discontinuities rather than genuine bistability)."""
    order = D_values[::-1] if reverse else D_values
    peak_Ts, X_overalls, C_exit_prev = [], [], None
    for D in order:
        X, Tp, U, Re, u, C_exit_prev = solve_recycle_pfr(D, RR, C_exit_init=C_exit_prev)
        peak_Ts.append(Tp)
        X_overalls.append(X)
    if reverse:
        peak_Ts, X_overalls = peak_Ts[::-1], X_overalls[::-1]
    return np.array(peak_Ts), np.array(X_overalls)


if __name__ == "__main__":
    print(f"Base (once-through) reactor sizing basis: tau={tau_oncethrough:.0f} s, L/D={L_over_D:.0f}")
    print(f"Solvent boiling point: {Tb_solvent - 273.15:.1f} C\n")

    D_values = np.linspace(0.005, 0.20, 40)
    RR_cases = [2, 5, 10]

    print("Continuation check: forward (increasing D) vs backward (decreasing D) sweeps.")
    print("If they agree, the earlier jumps were real; if they disagree/hysteresis-loop, it was")
    print("the naive fixed-point iteration landing on different branches inconsistently.\n")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for RR in RR_cases:
        peakT_fwd, X_fwd = continuation_sweep(D_values, RR, reverse=False)
        peakT_bwd, X_bwd = continuation_sweep(D_values, RR, reverse=True)

        max_diff = np.max(np.abs(peakT_fwd - peakT_bwd))
        print(f"RR={RR:3d}: max |peakT_fwd - peakT_bwd| = {max_diff:.2f} K "
              f"{'(AGREE -- smooth, no real bistability found)' if max_diff < 0.5 else '(DISAGREE -- real hysteresis/bistability)'}")

        runaway_idx = int(np.argmax(peakT_fwd > Tb_solvent)) if np.any(peakT_fwd > Tb_solvent) else None
        if runaway_idx:
            print(f"    critical D (forward) ~ {D_values[runaway_idx]*1000:.1f} mm")
        else:
            print(f"    no runaway up to D={D_values[-1]*1000:.0f} mm")

        axes[0].plot(D_values * 1000, peakT_fwd - 273.15, '-', lw=2, label=f'RR={RR} fwd')
        axes[0].plot(D_values * 1000, peakT_bwd - 273.15, '--', lw=1.5, label=f'RR={RR} bwd')
        axes[1].plot(D_values * 1000, X_fwd * 100, '-', lw=2, label=f'RR={RR} fwd')
        axes[1].plot(D_values * 1000, X_bwd * 100, '--', lw=1.5, label=f'RR={RR} bwd')

    axes[0].axhline(Tb_solvent - 273.15, color='k', ls=':', label='solvent b.p.')
    axes[0].set_xlabel('Reactor tube diameter D (mm)')
    axes[0].set_ylabel('Peak T (C)')
    axes[0].set_title('Peak T vs D: forward vs backward continuation')
    axes[0].legend(fontsize=7)

    axes[1].set_xlabel('Reactor tube diameter D (mm)')
    axes[1].set_ylabel('Overall (fresh-feed) conversion (%)')
    axes[1].set_title('Conversion vs D: forward vs backward continuation')
    axes[1].legend(fontsize=7)

    plt.tight_layout()
    plt.savefig('aniline_pfr_recycle_continuation.png', dpi=130)
    print("\nSaved plot to aniline_pfr_recycle_continuation.png")
