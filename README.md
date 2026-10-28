# Jacketed / Packed-Bed Reactor Scale-Up Calculator

How large can a tubular reactor running an exothermic reaction get before it stops being safe or practical?
Pick a reaction and an L/D ratio. The tool finds the tube diameter at which each physical limit binds, then
recommends a design point (diameter and throughput) 15 % below the first limit reached.

> **Make it yours, free:** download this repo (**Code → Download ZIP**, or `git clone`), open the folder in
> VS Code or any editor, and extend it for your own case with your AI assistant. This README is written for
> both of you. Taking these tools further is a free service of
> [Sharpes Process Engineering](https://jsharpes-portfolio.netlify.app); the site's Ask box can point you to the right starting file.
>
> **License:** [MIT](LICENSE) for the code and original text: use, modify and share freely, keeping the
> copyright notice. Provided as is, without warranty of any kind. It is an engineering study aid, not a
> substitute for a design checked and sealed by a licensed professional engineer.

**The main finding:** different reactions hit different limits.

- **Aniline (empty jacketed tube)** is limited by heat removal. Heat generation scales with volume (∝ D²L)
  and jacket removal with surface (∝ DL). Flow stays laminar (Re ≈ 4–375), and laminar Nu = 3.66 does not
  depend on velocity, so the critical diameter **does not change with L/D**. Only throughput moves.
- **Cumene (packed bed)** is usually limited by pressure drop. Packed-bed wall heat transfer (Leva,
  Nu_w ∝ Re_p^0.8) improves with velocity, so the thermal limit *rises* with L/D while the Ergun ΔP limit
  *falls*. They cross near L/D ≈ 140, which gives a real optimum.

- **Live tool:** https://jsharpes-portfolio.netlify.app/#rsu-tool
- **Run it locally:** download this repo and open `index.html`. Reaction data is precomputed into
  `aniline_calculator_data.js` and `cumene_calculator_data.js`, so the page runs no Python.
- **Author:** Joshua Sharpes, EIT. A screening study for teaching scale-up physics, not a reactor design basis.

> **For AI assistants:** `project/Reactor Scaleup Tool - Plan.md` is the detailed source of truth: goals,
> reaction data with citations, every confirmed result, and known caveats. The web page renders each
> reaction generically from its data file's `limits`, `tableRows`, `gaps`, `footer`, and
> `limitsExplanation` fields. To explain a result, read those fields for the selected reaction. The Python
> scripts in `project/` produce the data files and the validation plots.

## Method

1. **Mass balance fixes volume.** Kinetics plus target conversion set the required reactor volume,
   independent of tube geometry.
2. **L/D plus volume fix D and L.** Reactor type (empty tube or packed bed) only decides which
   correlations give U and ΔP.
3. **Energy balance along the PFR** (steady state, jacketed, scipy ODE) gives peak temperature vs. diameter.
4. **Limits, each reported separately and never blended:**
   - temperature ceiling crossing (for example, solvent boiling point, for aniline)
   - **Semenov / Frank-Kamenetskii** thermal-runaway onset
   - **Ergun** pressure drop vs. an allowable ΔP (10 % of operating pressure, a common heuristic flagged as unverified)
5. **Recommended design point:** 15 % margin below the smallest (binding) critical diameter, with the
   throughput at that diameter.

## Reaction data

| | Aniline (DuPont/KBR nitrobenzene hydrogenation) | Cumene (UOP Q-Max benzene alkylation) |
|---|---|---|
| ΔH_rxn | −552.9 kJ/mol (Hess's law, NIST WebBook formation enthalpies) | −110.6 kJ/mol (same method) |
| Kinetics | first order in nitrobenzene; Ea ≈ 20.6 kJ/mol, A = 1.93 s⁻¹ (fit to the source's data, Pd catalyst) | Ea = 52.6 kJ/mol, A = 6510 s⁻¹. DIPB side reaction Ea = 55.0 kJ/mol (sourced, not modeled) |
| Reactor | empty tube, laminar, jacketed | beta-zeolite packed bed, 2 mm pellets, 170–180 °C, ~3.5 MPa |
| Binding limit | thermal (L/D-invariant) | usually pressure drop; crossover near L/D ≈ 140 |

Reactions were shortlisted from Meyers, *Handbook of Petrochemical Production Processes*. Kinetics were
checked against the primary articles (see the Plan doc for the exact papers). Runaway temperatures are
derived with Semenov/F-K from the assembled kinetics, because tabulated values are generally proprietary.
The Semenov/F-K basis is Davis & Davis, *Fundamentals of Chemical Reaction Engineering*, Ch. 9.

## Validation and cross-checks

- **Cantera cross-check (aniline):** an independent Cantera reactor-network model reproduces the same trend.
  Its critical diameter is about 37.5 mm vs. 54.6 mm from the scipy model. The 30 % gap comes from
  Cantera's composition- and temperature-dependent properties (see `aniline_cantera_*.py` and the plots).
- **Semenov onset (aniline):** D ≈ 28.5 mm once-through, which is more conservative than the
  boiling-point crossing.
- **Recycle (aniline):** recycle at RR = 2 removes the boiling-point limit up to 200 mm. The Semenov onset
  moves to about 43 mm. Forward and backward continuation agree (`aniline_pfr_recycle*.py`).

## Limits

- Two reactions so far. EDC, methanol synthesis, and ethylbenzene are shortlisted but not built.
- Allowable ΔP = 10 % of operating pressure is a heuristic, not a cited value.
- Steady-state, plug-flow, 1-D model: no radial gradients, no catalyst deactivation.

## Files

```
index.html                    the tool, as on the live site
aniline_calculator_data.js    precomputed web data (generated by project/*_webdata_precompute.py)
cumene_calculator_data.js
project/
  Reactor Scaleup Tool - Plan.md          full build history, sources, results
  aniline_pfr_feasibility.py              base PFR energy-balance sweep
  aniline_pfr_recycle.py (+ continuation) recycle extension
  aniline_semenov_runaway.py              Semenov/F-K stability analysis
  aniline_cantera_check.py / _sweep.py    independent Cantera cross-check (+ aniline_mech.yaml)
  cumene_pfr_feasibility.py               packed-bed model (Leva U, Ergun ΔP, Semenov T*)
  *_webdata_precompute.py                 write the two data .js files
  reactor_scaleup_calculator_block.html   embeddable block; test_harness.html runs it on its own
  *.png                                   result plots
```

Python needs `numpy`, `scipy`, `matplotlib`, `thermo`, and (for the cross-check) `cantera`.
