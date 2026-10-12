"""
Precompute static D-sweep data (RR=2 recycle-loop aniline reactor, per L/D) for the web
calculator, baked at build time since the deployed site is static HTML/JS (no Python backend
on Netlify) -- per "Reactor Scaleup Tool - Plan.md", Workflow step 5 and Input/Output Contract.

Physics: reuses the validated base (aniline_pfr_feasibility.py) and Semenov
(aniline_semenov_runaway.py) primitives directly -- kinetics, NIST-derived dH_rxn, U(D,u)
correlation, and the exact Semenov/FK critical temperature T* -- so the web tool shows the same
numbers already cross-checked against Cantera and confirmed in the plan doc, not a re-derivation.

RR=2 (recycle ratio) is fixed/baked in, matching the validated recycle-loop finding ("Recycle Loop
-- Fixes the mm-Scale Result", plan doc) that RR=0 gives an unrealistically small ~55mm critical
diameter for a "scale-up" tool. The sweep here re-implements the recycle fixed-point solve inline
(rather than importing aniline_pfr_recycle.py's module-level L_over_D=100 constant) so L/D itself
can be varied across the grid the web slider will expose -- the one thing that script held fixed.

Fresh-feed throughput (mdot) turns out to be RR-independent for a given (D, L/D): the recycle loop
only changes whether that throughput is thermally safe (by pushing the critical diameter outward),
not the production basis itself -- u_fresh = L_over_D*D/tau_base regardless of RR, since tau_base
is fixed by target conversion and feed kinetics alone.

Output: aniline_calculator_data.js, a single `const ANILINE_DATA = {...}` assignment, hand-copied
into portfolio-site/index.html's Interactive Calculations section (not auto-injected -- the plan
doc's own deployment note: build/verify everything in this folder first, hand off once finished).
"""
import numpy as np
from scipy.integrate import solve_ivp
import json
import importlib.util
import matplotlib
matplotlib.use('Agg')

spec = importlib.util.spec_from_file_location("base", "aniline_pfr_feasibility.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

R, Ea, A_arr, dH_rxn = base.R, base.Ea, base.A_arr, base.dH_rxn
C0, Tc, T_feed = base.C0, base.Tc, base.T_feed
rho_feed, Cp_feed_mass = base.rho_feed, base.Cp_feed_mass
tau_base = base.tau_base
Tb_solvent = base.constants.Tbs[0]
RR = 2.0

# Exact Semenov/FK critical temperature (depends only on Ea, R, Tc -- see aniline_semenov_runaway.py)
disc = Ea**2 - 4 * R * Ea * Tc
T_star = (Ea - np.sqrt(disc)) / (2 * R)


def k1_crit_at(C):
    g = -dH_rxn * A_arr * np.exp(-Ea / (R * T_star)) * C
    return g * Ea / (R * T_star**2)


def integrate_pass(D, L, u, C_in):
    U, _ = base.U_overall(D, u)

    def odes(z, y):
        C, T = y
        k = A_arr * np.exp(-Ea / (R * T))
        rate = k * C
        dCdz = -rate / u
        gen = -dH_rxn * rate
        remove = (4 * U / D) * (T - Tc)
        dTdz = (gen - remove) / (rho_feed * Cp_feed_mass * u)
        return [dCdz, dTdz]

    sol = solve_ivp(odes, [0, L], [C_in, T_feed], max_step=L / 500)
    return sol.y[0, -1], sol.y[1].max()


def solve_recycle(D, L, u_total, C_exit_guess, tol=1e-8, max_iter=300, damping=0.5):
    for _ in range(max_iter):
        C_in = (C0 + RR * C_exit_guess) / (1 + RR)
        C_exit_new, T_peak = integrate_pass(D, L, u_total, C_in)
        step = C_exit_new - C_exit_guess
        if abs(step) < tol * C0:
            C_exit_guess = C_exit_new
            break
        C_exit_guess += damping * step
    X = 1 - C_exit_guess / C0
    C_in = (C0 + RR * C_exit_guess) / (1 + RR)
    return X, T_peak, C_exit_guess, C_in


def find_crossing(D_mm_arr, arr, thresh, rising=True):
    """Linear-interpolated D where arr crosses thresh (rising = arr goes from below to above)."""
    cond = (arr > thresh) if rising else (arr < thresh)
    if not cond.any() or cond.all():
        return None
    idx = int(np.argmax(cond))
    if idx == 0:
        return float(D_mm_arr[0])
    x0, x1 = D_mm_arr[idx - 1], D_mm_arr[idx]
    y0, y1 = arr[idx - 1], arr[idx]
    frac = (thresh - y0) / (y1 - y0)
    return float(x0 + frac * (x1 - x0))


LOD_GRID = [25, 50, 75, 100, 125, 150, 175, 200]
D_MM = np.linspace(5, 200, 60)
D_M = D_MM / 1000.0

series = {}
for L_over_D in LOD_GRID:
    peakT_C, convX_pct, mdot_kgh, k1_actual, k1_crit = [], [], [], [], []
    C_exit_prev = C0 * 0.5
    for D in D_M:
        L = L_over_D * D
        u_fresh = L_over_D * D / tau_base
        u_total = (1 + RR) * u_fresh
        X, T_peak, C_exit_prev, C_in = solve_recycle(D, L, u_total, C_exit_prev)
        U, Re = base.U_overall(D, u_total)
        area = np.pi * D**2 / 4
        mdot = u_fresh * area * rho_feed * 3600
        peakT_C.append(T_peak - 273.15)
        convX_pct.append(X * 100)
        mdot_kgh.append(mdot)
        k1_actual.append(4 * U / D)
        k1_crit.append(k1_crit_at(C_in))

    peakT_C = np.array(peakT_C)
    k1_actual = np.array(k1_actual)
    k1_crit = np.array(k1_crit)
    margin = k1_actual - k1_crit

    D_crit_boil = find_crossing(D_MM, peakT_C, Tb_solvent - 273.15, rising=True)
    D_crit_sem = find_crossing(D_MM, -margin, 0.0, rising=True)

    series[str(L_over_D)] = {
        "peakT_C": [round(v, 2) for v in peakT_C.tolist()],
        "convX_pct": [round(v, 2) for v in convX_pct],
        "mdot_kgh": [round(v, 4) for v in mdot_kgh],
        "limits": [
            {"key": "semenov", "label": "Semenov/F-K stability limit", "kind": "thermal",
             "D_crit_mm": round(D_crit_sem, 1) if D_crit_sem else None,
             "T_ref_C": round(T_star - 273.15, 1)},
            {"key": "boiling", "label": "Boiling-point limit (hard physical ceiling)", "kind": "thermal",
             "D_crit_mm": round(D_crit_boil, 1) if D_crit_boil else None,
             "T_ref_C": round(Tb_solvent - 273.15, 1)},
        ],
    }
    print(f"L/D={L_over_D}: D_crit_boiling={D_crit_boil}, D_crit_semenov={D_crit_sem}")

data = {
    "reaction_name": "Aniline (MNB Hydrogenation, DuPont/KBR-type)",
    "equation": "C6H5NO2(l) + 3 H2(g) → C6H5NH2(l) + 2 H2O(l)",
    "D_mm": [round(v, 2) for v in D_MM.tolist()],
    "LOD_grid": LOD_GRID,
    "series": series,
    "reference": {
        "dH_rxn_kJmol": round(dH_rxn / 1000, 1),
        "Ea_kJmol": round(Ea / 1000, 2),
        "A_arr_1s": round(A_arr, 4),
        "C0_molm3": C0,
        "solvent": "methanol",
        "Tb_solvent_C": round(Tb_solvent - 273.15, 1),
        "Tc_C": round(Tc - 273.15, 1),
        "Tfeed_C": round(T_feed - 273.15, 1),
        "X_target_pct": 95,
        "RR": RR,
        "T_star_semenov_C": round(T_star - 273.15, 1),
        "throughputUnit": "kg/h MNB",
        "limitsExplanation": "<strong>Two independent runaway criteria, on purpose.</strong> The "
            "boiling-point limit is where the actual reactor ODEs (with axial reactant depletion) "
            "predict the solvent starts flashing — a hard physical ceiling. The Semenov/"
            "Frank-Kamenetskii limit is a more conservative, first-principles thermal-stability bound "
            "(exact quadratic solve of the classic tangency criterion) that freezes concentration at "
            "the inlet value, so it flags local loss of steady-state uniqueness before the reactor as "
            "a whole necessarily overheats. They disagree because real reactant burn-off downstream "
            "self-limits the exotherm. <strong>The recommendation above applies a 15% margin below "
            "whichever of these two is smaller (the binding limit)</strong> — a standard conservative "
            "design margin, not a further first-principles derivation; treat it as a starting point "
            "for detailed design, not a final answer.",
        "tableRows": [
            ["ΔH_rxn (Hess's Law, NIST)", f"{round(dH_rxn/1000,1)} kJ/mol"],
            ["Activation energy, Ea", f"{round(Ea/1000,2)} kJ/mol"],
            ["Pre-exponential factor, A", f"{round(A_arr,4)} 1/s"],
            ["Feed nitrobenzene conc., C0", f"{C0:.1f} mol/m³"],
            ["Solvent", f"methanol (b.p. {round(Tb_solvent-273.15,1)} °C)"],
            ["Coolant temperature, Tc", f"{round(Tc-273.15,1)} °C"],
            ["Feed temperature", f"{round(T_feed-273.15,1)} °C"],
            ["Target conversion", "95 %"],
            ["Recycle ratio, RR", f"{RR}× fresh feed"],
            ["Semenov critical temp., T*", f"{round(T_star-273.15,1)} °C"],
            ["Kinetics source", "Prekob et al., Nanomaterials 11(5):1172 (2021), Pd/GCF catalyst, dilute nitrobenzene-in-methanol, 283-323 K"],
            ["ΔH_rxn source", "Hess's Law from NIST WebBook liquid-phase heats of formation"],
        ],
        "gaps": [
            "Kinetics (Ea, A) are fitted from Prekob et al. 2021's dilute lab system (Pd/glassy-carbon-foam "
            "catalyst, nitrobenzene in methanol, 283–323 K) — the best verified public data for this reaction "
            "class, not DuPont's proprietary neat-feed industrial catalyst. Throughput numbers are "
            "lab-kinetics-scale, not representative of real DuPont capacity.",
            "Recycle loop (RR=2, fixed) uses an idealized, perfectly effective external cooler bringing "
            "recycle back to feed temperature — not a sized real heat exchanger.",
            "Properties (ρ, Cp, μ, k) are held constant at feed conditions along the reactor. An independent "
            "Cantera cross-check (full composition/T-dependent properties) found the same qualitative "
            "runaway behavior but a ~30% lower critical diameter — a legitimate modeling uncertainty band, "
            "not a bug.",
            "L/D is offered at 8 discrete precomputed points (25–200), not a continuous solve — each point "
            "is a fully converged reactor sweep, not an interpolation.",
            "T_runaway (Semenov) uses the classical inlet-concentration approximation — see the caution "
            "note above for why it disagrees with the boiling-point number.",
            "The recommended tube diameter applies a flat 15% margin below the binding stability limit — "
            "a reasonable conservative-practice starting point, not an optimized or code-derived design "
            "margin (no formal HAZOP/LOPA-style margin study behind this number).",
            "The reaction-scheme diagram is a simplified skeletal structure for illustration (ring "
            "hydrogens omitted, not to bond-angle/3D scale).",
        ],
        "footer": "Ideal 1D PFR, scipy solve_ivp · NIST WebBook Hess's-Law ΔH_rxn · Prekob et al. 2021 "
                  "kinetics · Davis &amp; Davis Ch. 9 Semenov/FK stability theory · For engineering "
                  "demonstration use only",
    },
}

with open("aniline_calculator_data.js", "w", encoding="utf-8") as f:
    f.write("// Auto-generated by aniline_webdata_precompute.py -- do not hand-edit.\n")
    f.write("const ANILINE_DATA = ")
    f.write(json.dumps(data, indent=None))
    f.write(";\n")

print("Wrote aniline_calculator_data.js")
