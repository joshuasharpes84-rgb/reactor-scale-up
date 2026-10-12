"""
Precompute static D-sweep data (Cumene packed-bed reactor, per L/D) for the web calculator, baked
at build time -- same rationale as aniline_webdata_precompute.py (no Python backend on Netlify).

L/D grid is denser than Aniline's 8-point grid (25-200) between 100-150, since
cumene_pfr_feasibility.py found the two binding criteria (Semenov thermal limit, rising with L/D;
Ergun pressure-drop limit, falling with L/D) cross somewhere in that range -- the binding diameter
has a real maximum there, and the web slider needs enough resolution to actually show it, not step
past it between two sparse grid points.

Output schema (generalized from Aniline's original hardcoded boiling/semenov fields, since Cumene
introduces a genuinely different limit TYPE -- pressure-drop, not another thermal criterion):
  series[lod] = {
    peakT_C, convX_pct, mdot_kgh, dP_kPa: [...],
    limits: [ {key, label, kind, D_crit_mm}, ... ]   -- generic, UI iterates this, no hardcoded names
  }
kind is "thermal" (has a meaningful T* to plot as a horizontal reference line) or "pressure" (no
temperature reference; only a vertical D marker makes sense).
"""
import numpy as np
import json
import importlib.util

spec = importlib.util.spec_from_file_location("base", "cumene_pfr_feasibility.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

R, Ea, A_arr, dH_rxn = base.R, base.Ea, base.A_arr, base.dH_rxn
C0, Tc, T_feed = base.C0, base.Tc, base.T_feed
rho_feed, Cp_feed_mass = base.rho_feed, base.Cp_feed_mass
tau_base = base.tau_base
T_star = base.T_star
dP_allowable = base.dP_allowable
P_op = base.P_op


def find_crossing(D_mm_arr, arr, thresh, rising=True):
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


LOD_GRID = [25, 50, 75, 100, 110, 120, 130, 140, 150, 160, 175, 200, 225, 250]
# D_MIN = 50mm: enforces D/dp >= 25 (dp=2mm), standard packed-bed design practice (D/dp >= 20) to
# avoid wall-channeling -- see cumene_pfr_feasibility.py's matching note for the full finding.
D_MM = np.linspace(50, 300, 60)
D_M = D_MM / 1000.0

k1_crit_base = base.k1_crit_at(C0)  # single-pass, no recycle -> C fixed at C0

series = {}
for L_over_D in LOD_GRID:
    peakT_C, convX_pct, mdot_kgh, dP_kPa, k1_actual = [], [], [], [], []
    for D in D_M:
        sol, u_s, U, Re_p, L, dP = base.integrate_pfr(D, L_over_D, tau_base)
        X = 1 - sol.y[0, -1] / C0
        T_peak = sol.y[1].max()
        area = np.pi * D**2 / 4
        mdot = u_s * area * rho_feed * 3600  # kg/h propylene-equivalent, fresh-feed basis
        peakT_C.append(T_peak - 273.15)
        convX_pct.append(X * 100)
        mdot_kgh.append(mdot)
        dP_kPa.append(dP / 1000)
        k1_actual.append(4 * U / D)

    peakT_C = np.array(peakT_C)
    dP_kPa = np.array(dP_kPa)
    margin = np.array(k1_actual) - k1_crit_base

    D_crit_pressure = find_crossing(D_MM, dP_kPa, dP_allowable / 1000, rising=True)
    D_crit_semenov = find_crossing(D_MM, -margin, 0.0, rising=True)

    series[str(L_over_D)] = {
        "peakT_C": [round(v, 2) for v in peakT_C.tolist()],
        "convX_pct": [round(v, 2) for v in convX_pct],
        "mdot_kgh": [round(v, 4) for v in mdot_kgh],
        "dP_kPa": [round(v, 2) for v in dP_kPa.tolist()],
        "limits": [
            {"key": "semenov", "label": "Semenov/F-K stability limit", "kind": "thermal",
             "D_crit_mm": round(D_crit_semenov, 1) if D_crit_semenov else None,
             "T_ref_C": round(T_star - 273.15, 1)},
            {"key": "pressure", "label": "Pressure-drop limit (10% of operating P)", "kind": "pressure",
             "D_crit_mm": round(D_crit_pressure, 1) if D_crit_pressure else None},
        ],
    }
    print(f"L/D={L_over_D}: D_crit_semenov={D_crit_semenov}, D_crit_pressure={D_crit_pressure}")

data = {
    "reaction_name": "Cumene (UOP Q-Max Benzene Alkylation)",
    "equation": "C6H6(l) + C3H6(dissolved) → C6H5CH(CH3)2(l)",
    "D_mm": [round(v, 2) for v in D_MM.tolist()],
    "LOD_grid": LOD_GRID,
    "series": series,
    "reference": {
        "dH_rxn_kJmol": round(dH_rxn / 1000, 1),
        "Ea_kJmol": round(Ea / 1000, 2),
        "A_arr_1s": A_arr,
        "C0_molm3": round(C0, 1),
        "Tfeed_C": round(T_feed - 273.15, 1),
        "Tc_C": round(Tc - 273.15, 1),
        "P_op_kPa": round(P_op / 1000, 0),
        "dP_allowable_kPa": round(dP_allowable / 1000, 0),
        "X_target_pct": 90,
        "T_star_semenov_C": round(T_star - 273.15, 1),
        "dP_allowable_kPa": round(dP_allowable / 1000, 0),
        "throughputUnit": "kg/h propylene",
        "limitsExplanation": "<strong>Two criteria of genuinely different character, on purpose.</strong> "
            "This reaction's exotherm (−110.6 kJ/mol) is much milder than Aniline's, and packed-bed "
            "wall heat transfer (Leva's correlation) is far more effective than an empty tube's — so "
            "the Semenov/F-K thermal-stability limit here rarely binds at practical diameters. What "
            "does bind is <strong>pressure drop</strong> (Ergun equation): as L/D rises, velocity "
            "rises, and ΔP grows steeply — while the thermal limit simultaneously grows too (higher "
            "velocity also means better cooling). The two move in <strong>opposite directions</strong> "
            "as L/D changes, so which one binds — and the diameter the recommendation lands on — "
            "genuinely depends on L/D, unlike Aniline's flat, diameter-invariant case. "
            "<strong>The recommendation above applies a 15% margin below whichever is smaller "
            "(the binding limit) at the selected L/D</strong> — a standard conservative design "
            "margin, not a further first-principles derivation.",
        "tableRows": [
            ["Reaction basis", "propylene-limiting, C6H6:C3H6 = 0.85:1 (molar)"],
            ["ΔH_rxn (Hess's Law, NIST)", f"{round(dH_rxn/1000,1)} kJ/mol"],
            ["Activation energy, Ea", f"{round(Ea/1000,2)} kJ/mol"],
            ["Pre-exponential factor, A", f"{A_arr:.0f} 1/s"],
            ["Feed propylene conc., C0", f"{round(C0,1)} mol/m³"],
            ["Feed temperature", f"{round(T_feed-273.15,1)} °C"],
            ["Operating pressure", f"{round(P_op/1000,0):.0f} kPa"],
            ["Coolant temperature, Tc", f"{round(Tc-273.15,1)} °C"],
            ["Target conversion", "90 %"],
            ["Semenov critical temp., T*", f"{round(T_star-273.15,1)} °C"],
            ["Allowable ΔP (10% of P_op)", f"{round(dP_allowable/1000,0):.0f} kPa"],
            ["Catalyst", "2mm beta-zeolite pellets, Si/Al=25"],
            ["Kinetics source", "Ben Hamida, Preprints.org 202601.1995.v1 (2026), Aspen HYSYS PFR sim, Q-Max"],
            ["ΔH_rxn source", "Hess's Law from NIST WebBook liquid-phase heats of formation"],
        ],
        "gaps": [
            "Kinetics (K1) are from an Aspen HYSYS PFR simulation study fitted to the real Q-Max process, "
            "treated here as an already bed-effective (lumped) rate constant -- no separate catalyst "
            "effectiveness factor is applied, to avoid double-counting mass-transfer limits already "
            "implicit in the fit.",
            "DIPB transalkylation side reaction has sourced kinetics (K2 = 450·exp(-55000/RT)) but is not "
            "yet modeled -- this version tracks only the main alkylation.",
            "Coolant temperature (Tc = Tfeed - 20K) and target conversion (90%) are reasonable engineering "
            "assumptions, not sourced from DuPont/UOP-specific design data.",
            "Allowable pressure drop (10% of operating pressure) is a commonly-cited packed-bed design "
            "heuristic, not verified against a specific page in this project's own reference library -- "
            "flag if it needs firming up before using this number for real design work.",
            "Benzene:propylene ratio (0.85:1) implies propylene in slight molar excess over benzene, taken "
            "at face value from the source -- unusual relative to typical industrial practice (excess "
            "benzene normally suppresses polyalkylation), not independently re-verified here.",
            "Properties (ρ, Cp, μ, k) held constant at feed conditions along the reactor -- same "
            "simplification as the Aniline case, not cross-checked against a composition-aware solver "
            "(e.g. Cantera) the way Aniline was.",
        ],
        "footer": "Ideal 1D packed-bed PFR, scipy solve_ivp · NIST WebBook Hess's-Law ΔH_rxn · "
                  "Ben Hamida 2026 kinetics · Leva packed-bed wall-U correlation · Ergun ΔP · "
                  "Davis & Davis Ch. 9 Semenov/FK stability theory · For engineering demonstration use only",
    },
}

with open("cumene_calculator_data.js", "w", encoding="utf-8") as f:
    f.write("// Auto-generated by cumene_webdata_precompute.py -- do not hand-edit.\n")
    f.write("const CUMENE_DATA = ")
    f.write(json.dumps(data, indent=None))
    f.write(";\n")

print("Wrote cumene_calculator_data.js")
