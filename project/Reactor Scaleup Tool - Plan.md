# Reactor Scale-Up Tool — Plan

## Goal

An interactive Netlify calculator (portfolio piece, `jsharpes-portfolio.netlify.app`) demonstrating
the generic process of scaling a **jacketed PFR running an exothermic reaction** from lab scale to
pilot scale. Originated as a spin-off of `Design\Zeta Energy`'s axial-dispersion/nonideal-PFR
modeling anchor (built while that project was still active) but applied to common, named
petrochemical reactions rather than Zeta's own SPAN chemistry — moved into this standalone folder
2026-08-17 since it was never Zeta-specific to begin with.

## Workflow

1. User selects a reaction from a dropdown — common reactions drawn from
   `Reference Materials\Handbook for Petrochemical Production Processes.pdf` (Meyers).
2. For the selected reaction: identify reacting species, reaction properties (ΔH_rxn, T_runaway,
   kinetics), and typical operating conditions. Requires supporting literature — see sourcing notes
   below.
3. Reaction properties + operating conditions give the utility type for cooling, completing the
   energy balance.
4. User inputs L/D ratio of the reactor; tool computes maximum scale-up size (as product mass flow
   rate) before heat-removal capacity (UA/V, shrinking with diameter) can no longer keep up with the
   heat-generation rate (scales with V) — this ratio is the actual engineering core of the tool.
5. Fluid properties pulled in via Python (`thermo` library, PR/SRK EOS per this machine's standard
   pattern) — baked into a static data table at build time, since the deployed site is static
   HTML/JS/CSS (no Python backend on Netlify).
6. Other visualization ideas (graphs, trends) TBD as the tool takes shape.

## Deployment

Build everything in this folder first (`Reactor Scaleup Tool\`). Only hand the finished
HTML/JS/CSS/data files to `Interview Prep\portfolio-site` for deploy, to avoid burning Netlify
deploy tokens on iteration.

## Reaction Property Sourcing

Per this project's own `CLAUDE.md` convention (Confirmed vs. Inference) and the Heat Exchangers
project's source-hierarchy discipline (own library → shared library → sibling SME folder → web,
last resort and always verified against the primary source before use):

- **ΔH_rxn**: prefer deriving via Hess's Law from standard heats of formation (ΔHf°) pulled from a
  primary source (e.g. NIST WebBook) rather than trusting a secondary "heat of reaction" claim
  wholesale — verifiable and reproducible.
- **T_runaway**: not generally available as a tabulated per-reaction literature number outside
  proprietary/instrument-derived data (ARC, DSC, VSP2 calorimetry — not owned in this library).
  Better approach: *derive* a critical/runaway temperature per reaction from Semenov /
  Frank-Kamenetskii thermal-stability criteria, using assembled kinetics (Ea, A), ΔH_rxn, and the
  cooling UA — this is the rigorous version of what the tool is already computing, not a separate
  lookup.
- Methodology reference already in this project's own library: `Reference Library\Reactor Design &
  Scale-Up\Reaction Kinetics Textbooks.pdf` (Davis & Davis, Ch. 9 "Nonisothermal Reactors — energy
  balances, stability/sensitivity") — the standard basis for the Semenov-type approach above.
- Perry's (Sec. 4 Thermodynamics, Sec. 7 Reaction Kinetics) and Elliot & Lira (energy balances for
  reacting systems, Ch. 3; reaction equilibria, Ch. 17) are the fallback deep-reference sources,
  both already in the shared `Textbooks (Level 1)` library.
- Any new source beyond what's already owned (e.g. Bretherick's *Handbook of Reactive Chemical
  Hazards*, a DIERS/CCPS runaway-reaction dataset) should be flagged to the user for a vetting/
  acquisition decision before being pulled in unilaterally.

## Reaction Shortlist (confirmed 2026-08-13)

Screened all 18 Parts of the Meyers handbook for genuinely exothermic, PFR/fixed-bed-type
reactions with public (not just proprietary-vendor) kinetics. Ruled out: styrene (EB
dehydrogenation, endothermic), ethylene cracking (endothermic pyrolysis), propane dehydrogenation
(endothermic), and the polymer parts (PE/PP/PS/PET — polymerization systems, don't reduce to a
single-reaction jacketed-PFR energy balance cleanly).

**Building first, in this order:**

1. **Aniline** (Part 2, `Handbook...pdf` p.2.3–2.8) — DuPont/KBR process.
   Reaction: C₆H₅NO₂ (MNB) + 3H₂ → C₆H₅NH₂ (aniline) + 2H₂O, catalytic, liquid-phase.
   **Confirmed from the handbook: runs in a liquid-phase plug-flow hydrogenation reactor** — direct
   match for this tool's jacketed-PFR framing, no analogy needed. Also the flagship reaction
   because nitrobenzene hydrogenation has a well-known real-world thermal-runaway case history in
   process safety literature — good narrative fit. No numeric ΔH_rxn stated in the handbook itself;
   to be derived via Hess's Law from NIST WebBook heats of formation (ΔHf°) for MNB, H₂, aniline,
   H₂O rather than trusted from a secondary source.

2. **Cumene** (Part 4, `Handbook...pdf` p.4.3–4.13) — **two competing technologies in the
   handbook, use Q-Max (UOP, Ch. 4.2), not CDCumene (ABB Lummus/CDTECH, Ch. 4.1):**
   - CDCumene is a **catalytic distillation column**, not a PFR — heat of reaction is removed by
     vaporizing benzene in-situ, and the chapter explicitly states this makes "the risk of a
     runaway reaction... virtually nonexistent." Real and interesting (same self-regulating-cooling
     idea as the boiling-jacket point in this project's own findings doc), but the wrong reactor
     type for this tool — don't build the cumene case on it.
   - **Q-Max** is a conventional Friedel-Crafts alkylation: C₆H₆ (benzene) + C₃H₆ (propylene) ⇌
     C₆H₅CH(CH₃)₂ (cumene), plus a DIPB transalkylation side reaction
     (C₆H₄[CH(CH₃)₂]₂ + C₆H₆ → 2 cumene). This is the conventional fixed-bed case — use this one.
     No numeric ΔH_rxn stated in the handbook; same NIST/Hess's-Law derivation plan as Aniline.

**Queued for after the pipeline is proven on the above two:** EDC direct chlorination (Part 18),
methanol synthesis (Part 7), possibly ethylbenzene (Part 5, same alkylation family as cumene —
mostly useful to confirm the tool generalizes, likely redundant otherwise).

## Reaction Data (confirmed 2026-08-14)

### Aniline — DuPont/KBR MNB hydrogenation

C₆H₅NO₂(l) + 3H₂(g) → C₆H₅NH₂(l) + 2H₂O(l)

- **ΔH_rxn = −552.9 kJ/mol** (nitrobenzene basis). Derived via Hess's Law from NIST WebBook liquid-
  phase ΔfH° values (internally consistent source pairs, not mixed across studies): nitrobenzene(l)
  12.5 kJ/mol [Pedley 1986 reanalysis of Lebedeva et al. 1971], aniline(l) 31.3 kJ/mol [Hatton,
  Hildenbrand et al. 1962], H₂O(l) −285.83 kJ/mol [CODATA, Cox/Wagman et al. 1984], H₂(g) = 0
  (reference state). Consistent with the commonly-cited ~−550 kJ/mol figure for this reaction in
  process safety literature.
- **Kinetics**: first-order in nitrobenzene, **Ea ≈ 20–23 kJ/mol** (Pd catalyst; 20.18–22.9 kJ/mol
  depending on support). Source: Prekob et al., *Nanomaterials* 11(5):1172 (2021), open access via
  PMC8146346 — verified by fetching the primary article directly, not a search-summary claim.
  Conditions studied: 283–323 K, 20–50 bar, nitrobenzene in methanol solvent, Pd/glassy-carbon-foam
  catalysts. **Caveat**: not the exact DuPont proprietary catalyst or solvent-free industrial
  conditions — the best verified public literature value for the reaction class, not a vendor-exact
  number. Discarded several unverified search-summary Ea claims (35.81, 40.18, 35±1 kJ/mol from
  other studies) since they couldn't be confirmed against a fetchable primary source.

### Cumene — UOP Q-Max benzene alkylation

C₆H₆(l) + C₃H₆(g) ⇌ C₆H₅CH(CH₃)₂(l), plus DIPB side/transalkylation chemistry

- **ΔH_rxn = −110.6 kJ/mol**. Derived via Hess's Law from NIST WebBook values: benzene(l) 49.0 kJ/mol
  and cumene(l) −41.2 kJ/mol [Roux/Temprado 2008 and Prosen et al. 1945 respectively], propylene(g)
  20.41 kJ/mol [Furuyama/Golden 1969] — gas-phase for propylene since its critical temperature
  (91.8°C) is below reactor operating temperature, so it can't exist as a distinct liquid phase there.
  Consistent with the commonly-cited ~−113 kJ/mol figure.
- **Kinetics**: pseudo-first-order in propylene, rate = K₁·C_propylene, **K₁ = 6510·exp(−52564/RT)**
  → **Ea = 52.6 kJ/mol**, A = 6510 s⁻¹ (E in kJ/kmol, rate in kmol/m³·s, per the source's own units
  statement). Also gives the DIPB side reaction: K₂ = 450·exp(−55000/RT), Ea = 55.0 kJ/mol.
  Source: Ben Hamida, *Preprints.org* preprint 202601.1995.v1 (2026), read in full (not just
  search-summarized) — an Aspen HYSYS PFR simulation study citing Hazmi 2022 (*JCPE*) and
  Norouzi & Fatemi 2012 (*Chem. Eng. Commun.*) as the original kinetics sources. This paper's own
  reference list cites Al-Ghamdi et al.'s "Production of Cumene via Alkylation of Benzene and
  Propylene By Honeywell UOP QMax Technology" directly, confirming alignment with the Q-Max
  technology already selected in this plan, not a generic/mismatched zeolite process.
  Reactor conditions given: 170–180°C, 3540–3600 kPa, beta-zeolite (Si/Al=25), 2mm pellets,
  liquid-phase fixed-bed PFR, benzene:propylene molar ratio ~0.75–0.95:1 — usable as the tool's
  default base-case operating point.
  Cross-check: pseudo-first-order-in-olefin is also what Ercan et al. 1998 (cited in this project's
  own `Reaction Kinetics Textbooks.pdf`, p.235–236) found for the closely related zeolite-Y-catalyzed
  alkylation of ethylbenzene with propylene — same reaction family, same rate-law form, independent
  source.

## Feasibility Proof — Aniline (confirmed 2026-08-14)

Built and ran `aniline_pfr_feasibility.py` — a full 1D ideal-PFR mass+energy balance (scipy
`solve_ivp`), proving the tool's core physics end-to-end on real, verified data before building
the web UI or adding more reactions. Result: **the concept works.**

- Model: constant-property PFR, kinetics `k(T)=A*exp(-Ea/RT)` fitted directly from Prekob et al.'s
  4 real (T,k) data points (A=1.93 s⁻¹, Ea=20.6 kJ/mol — matches their own reported value),
  ΔH_rxn=−552.9 kJ/mol, Dittus-Boelter process-side h combined with an assumed 2000 W/m²K
  jacket-side h, sensible-heat coolant at 15°C.
- **Key finding, and an important self-consistency fix**: a first attempt modeled a *neat*
  nitrobenzene feed (matching the real DuPont process's likely concentration) and got immediate,
  nonphysical runaway even at the smallest lab diameter — undiluted nitrobenzene's exotherm at
  these rate constants vastly exceeds what any simple single-pass jacketed tube can remove. This
  is itself a real, honest result (it's exactly why real industrial hydrogenations this exothermic
  run diluted or with product recycle, never neat single-pass) but it mixed Prekob's dilute-system
  kinetics with a concentration Prekob never tested — not self-consistent. Fixed by modeling
  Prekob's *actual* system instead: dilute nitrobenzene (250 mol/m³) in methanol solvent, using
  methanol's own bulk ρ/Cp/μ/k from `thermo`.
- **Result on the corrected, self-consistent basis**: at L/D=100 and a 95%-conversion-targeted
  residence time (~71 min), peak reactor temperature stays flat at feed temperature for small
  diameters, then rises smoothly as D grows (cooling surface/volume ∝ 1/D degrading), crossing
  methanol's boiling point (64.5°C) at **D≈55 mm**. Max safe MNB throughput at that point:
  **~8.4 kg/h**, vs. ~400+ kg/h if cooling were unconstrained at D=200mm — a clean, physically
  grounded runaway boundary, and a real demonstration of the D³-throughput-vs-1/D-cooling tension
  that's the whole point of this tool.
- **What this proves**: the full pipeline (reaction data → kinetics → real thermo properties →
  PFR energy balance → L/D-driven scale-up limit → visualization) works end-to-end with real,
  citable numbers, not placeholders. Throughput numbers themselves are lab-kinetics-scale, not
  representative of DuPont's real industrial capacity — expected, since the only public kinetics
  available is a dilute lab system, not DuPont's proprietary neat-feed catalyst.
- Script and its output plot (`aniline_pfr_feasibility.png`) live in this folder.

## Input/Output Contract (confirmed 2026-08-14)

Deliberately minimal, matching the style of the existing portfolio calculators (compressor-train
sizing deck: a couple of free inputs, everything else cascades from fixed data + internal checks
against a physical limit) — the visitor isn't meant to spend time tuning this.

**Only two user inputs:**
1. **Reaction** (dropdown) — loads a fixed per-reaction dataset (see below).
2. **L/D ratio** (slider, bounded per-reaction to a sensible range so no nonsensical combination
   is reachable).

**Everything else is fixed, reaction-coupled data, not exposed as a slider:** feed concentration/
dilution, coolant type (sensible-liquid vs. boiling/evaporative), coolant temperature, target
conversion, feed temperature, and (confirmed 2026-08-14) **reactor type** (empty-tube/homogeneous
vs. packed-bed). Selecting a reaction *is* selecting that whole bundle — e.g. Aniline's entry locks
in the validated dilute-methanol basis, sensible-liquid jacket at 15°C, 95% target conversion,
40°C feed, empty-tube heat transfer (Dittus-Boelter/Gnielinski), exactly as proven out in
`aniline_pfr_feasibility.py`. Cumene will instead carry packed-bed reactor type (Ergun pressure
drop, a packed-bed wall-U correlation, catalyst effectiveness factor η₀ — matching the real Q-Max
bed specs already in hand: 2mm beta-zeolite pellets, Si/Al=25, bed density 0.70 g/cm³).

**Architectural principle (confirmed 2026-08-14): reactor type affects heat transfer (and other
hardware-specific physics — pressure drop, effectiveness factor), never the mass balance.** The
kinetics + target conversion determine the required reactor *volume* (residence time × flow) on
their own, independent of tube shape or reactor type — chemistry doesn't know whether it's sitting
in a packed bed or an empty tube. L/D + that volume then fix the geometry (D, L). Reactor type only
enters afterward, in which correlation family computes U (and, for packed beds, ΔP and η₀) at that
geometry. This is already how the code is structured (`U_overall()` is the only reactor-type-
dependent piece; the mass balance and volume/residence-time sizing never reference it) — stated
here explicitly so it doesn't drift when Cumene's packed-bed correlations get added.

**Outputs**, per reaction + L/D combination:
- Headline number: max safe throughput (kg/h)
- Critical diameter D_max
- Plot: peak reactor T vs. D (with the physical runaway limit marked)
- Plot: exit conversion vs. D
- Plot: throughput vs. D (with the cutoff marked)
- Reference panel echoing the reaction's fixed data (ΔH_rxn, Ea, A, MW, Tb, etc.)

## Independent Cross-Check via Cantera (confirmed 2026-08-14)

Built a second, independent model of the same aniline base case and full diameter sweep using
Cantera's actual reactor-network solver (not another hand-rolled ODE), to validate the scipy
model's physics rather than just trusting one implementation.

- **Method**: steady PFR modeled as a fluid parcel evolving in time (t=z/u, the standard Cantera
  PFR trick), using a `ConstPressureReactor` + `Wall`-to-`Reservoir` for jacket cooling. Reaction
  and thermo defined in `aniline_mech.yaml`, a custom `ideal-condensed` phase (5 species: methanol
  solvent + nitrobenzene, H2, aniline, water) with the *same* A/Ea (fitted from Prekob et al.) and
  the *same* NIST-derived species formation enthalpies as the scipy model — so only the reactor
  integration itself differs, not the underlying data.
- **Real gotcha, worth remembering**: Cantera's `ideal-condensed` phase rate law uses mole
  fraction (X), not molar concentration, as its activity basis — the forward rate *constant*
  matched the scipy model exactly (confirming the Arrhenius parameters parsed correctly), but the
  rate of progress was ~100x too low until `A` was rescaled by the mixture's total molar density
  (`A_yaml = A_scipy * C_total`). First rescale attempt used C_total in mol/m³ and overshot by
  1000x — Cantera's native molar unit is kmol, not mol — caught by comparing
  `net_rates_of_progress` against the scipy model's expected rate directly rather than trusting
  the fix on the first try.
- **Result**: Cantera reproduces the same qualitative finding as the scipy model — peak T rises
  smoothly with D, crossing methanol's boiling point, conversion asymptoting toward 100% — with a
  critical diameter of **~37.5 mm**, vs. the scipy model's **~54.6 mm**. Same shape, ~30% lower
  threshold. The gap has a clean, understood cause: Cantera tracks real composition/T-dependent
  enthalpy as the reaction proceeds; the scipy model held ρ and Cp fixed at feed conditions
  throughout. Cantera is arguably the more accurate of the two on this point, not just a duplicate.
- **What this proves**: the reactor-network physics (hot-spot growth, a genuine finite critical
  diameter, the D³-throughput-vs-1/D-cooling tension) is real and robust across two independent
  numerical implementations, not an artifact of one ODE setup. The ~30% quantitative gap between
  constant-property and composition-aware treatments is itself a legitimate, reportable modeling
  uncertainty band, not a bug to hide.
- Files: `aniline_mech.yaml`, `aniline_cantera_check.py` (base-case check),
  `aniline_cantera_sweep.py` (full sweep + comparison plot), `aniline_cantera_sweep_compare.png`.

## Recycle Loop — Fixes the mm-Scale Result (confirmed 2026-08-14)

The base (once-through) aniline model's critical diameter (~55 mm, ~8.4 kg/h) is unrealistically
small for a "scale-up" tool. Diagnosed why: at that scale the reactor is deeply laminar (Re~4-125
across the whole runaway range), so the process-side film coefficient is pinned at the laminar
constant (Nu=3.66, h~k/D) — cooling capacity was degrading as ~1/D² (laminar h *and* the usual
4/D surface-to-volume falloff), not just 1/D. Tried switching to a boiling/evaporative jacket
first (motivated by the confirmed real DuPont steam-export fact) — barely moved the critical
diameter (55.0mm), because the jacket side was never the bottleneck; U was ~70 W/m²K either way
since the process side dominated the resistance.

**Real fix: a recycle loop.** Reactor effluent splits into product + recycle; recycle passes
through an idealized external cooler back to T_feed, remixes with fresh feed
(C_in=(C_A0+RR·C_exit)/(1+RR), solved by fixed-point iteration since C_in depends on C_exit),
re-enters the same physical reactor (D, L/D unchanged) at total flow (1+RR)× higher. This attacks
both sides at once: recycle dilution lowers reactor-inlet concentration (directly cutting local
generation), and the higher velocity breaks the reactor out of laminar flow, which disproportionately
helps *removal* over generation (U grows with Re on top of the shared 1/u term both sides get).
**Result: even RR=2 (recycle = 2× fresh feed) pushes the runaway limit past the entire swept
200mm range**, vs. 55mm at RR=0 — confirmed correct, not a fluke, by extending the same
kinetics/thermo/U(D,u) machinery from the base model, changing only the recycle topology
(`aniline_pfr_recycle.py`).

**Caught and ruled out a possible false finding along the way**: the first sweep showed sharp
discontinuities in peak T vs. D (e.g. a clean drop from ~56°C to ~45°C at RR=2 between D=135-140mm).
Two explanations were live: genuine bistability (real, well-documented for autothermal/recycle
reactors) vs. a naive-solver artifact (successive substitution landing on different fixed points at
different D). Resolved by continuation — re-solving each D using the previous D's converged
solution as the starting guess, swept both forward (increasing D) and backward (decreasing D): the
two directions agreed to within 0.00 K everywhere, ruling out bistability. **The jumps are real,
but they're a correlation artifact, not physics or a solver bug**: they land exactly where Re
crosses 2300, the hard if/else threshold in `U_overall()` between the laminar constant (Nu=3.66)
and the turbulent Dittus-Boelter correlation. Real transitional flow is smooth; this correlation
isn't. Doesn't change the safety conclusion (no jump crosses the boiling-point line), but the exact
curve shape near those transition points needs a smooth correlation (e.g. Gnielinski) before it's
quantitatively trustworthy — open item.

**Closed (2026-08-14)**: replaced the hard if/else Nu switch in `U_overall()` with a smooth cubic
blend (smoothstep) between the laminar constant (Nu=3.66) and the Gnielinski correlation
(smoother and more accurate near the transition than Dittus-Boelter) across Re 2300-4000. Base
model unaffected (never leaves deep laminar flow before hitting its own limit). Re-ran the
recycle continuation sweep: the cliff-edge jumps are gone, replaced by smooth, physically sensible
bell-shaped peaks (T rises as D grows, then tapers back down as the reactor moves further into
turbulent flow and cooling keeps winning) — still zero crossings of the boiling-point line for
RR=2/5/10 across the full 0-200mm sweep, and forward/backward continuation still agree exactly
(bistability finding unaffected by the correlation fix, as expected since they're independent
issues). `U_overall()` is shared by both the base and recycle scripts (recycle imports it from
the base module), so this fix applies everywhere in one place.

Files: `aniline_pfr_recycle.py` (recycle model + first sweep), continuation check added in the
same file (`continuation_sweep()`), `aniline_pfr_recycle_continuation.png` (forward/backward
overlay, confirms no bistability).

## Semenov/Frank-Kamenetskii Thermal-Runaway Analysis — Aniline (confirmed 2026-08-15)

Derived T_runaway analytically (Davis & Davis Ch. 9 methodology, per the sourcing plan above)
rather than trusting a tabulated number, and cross-checked it against the numerical boiling-point
criterion already used in the base/recycle sweeps. `aniline_semenov_runaway.py` imports the base
and recycle modules directly (same pattern as `aniline_pfr_recycle.py`) so the underlying kinetics/
thermo/U(D,u) data is identical — only the stability analysis is new.

- **Method**: the PFR's per-unit-volume energy balance already has the classic Semenov form —
  heat generation g(T)=(-ΔH_rxn)·k(T)·C (exponential in T) vs. heat removal r(T)=(4U/D)(T-Tc)
  (linear in T), since 4U/D = UA/V for a cylindrical tube. The stability boundary is where the two
  curves are simultaneously equal AND tangent (g=r, g'=r'). Because g'(T)/g(T) = Ea/(RT²) exactly
  (Arrhenius derivative, independent of the A/ΔH_rxn/C prefactor), combining the two tangency
  conditions cancels the cooling term entirely and collapses to a quadratic in T* alone:
  R·T*² − Ea·T* + Ea·Tc = 0. Solved exactly (not the further-linearized classical approximation
  ΔT_crit=R·Tc²/Ea, which requires Ea/(R·Tc) >> 1 — checked and found NOT to hold well here,
  Ea/(R·Tc)=8.6, giving a 25% gap between the exact and linearized answers).
- **Result**: T* = 59.7°C (ΔT_crit = 44.7 K above the 15°C coolant), independent of D, U, C, or
  ΔH_rxn magnitude — those parameters instead set whether a given D is actually AT that critical
  point, via the required (4U/D)_crit = g'(T*), compared against the reactor's actual (4U/D) at
  each D.
- **Base (once-through) case**: Semenov instability onset at **D≈28.5 mm** — smaller (more
  conservative) than the 54.6 mm boiling-point crossing already found.
- **Recycle case (RR=2)**: Semenov instability onset at **D≈43.2 mm** — inside the 0–200 mm range
  the boiling-point sweep called fully safe (no crossing found up to 200 mm).
- **Important caveat, not a contradiction of the recycle fix**: the classical Semenov criterion
  freezes reactant concentration at its inlet value (a "no burn-off yet" worst case) — it flags
  *local* loss of steady-state uniqueness near the reactor entrance. The actual PFR ODEs let
  concentration deplete downstream, self-limiting the exotherm before temperature can actually
  diverge, so a reactor can be Semenov-unstable at the inlet in the formal multiplicity sense while
  still having a bounded peak temperature overall. This is a well-documented real feature of
  tubular reactors (why Barkelew/Bilous-Amundson tubular-reactor runaway criteria differ from the
  zero-order batch/CSTR Semenov result) — reported here as an honest second, more conservative
  criterion, not used to overturn the recycle-loop finding above, and not swept under the rug
  either. **Practical read: the recycle loop's boiling-point safety margin is real, but it has less
  margin than the once-through comparison alone suggested — D≈43mm is where a rigorous stability
  analysis says caution should start, well before the ~200mm boiling limit.**
- Files: `aniline_semenov_runaway.py`, `aniline_semenov_runaway.png` (tangency plot at the base
  case's critical D + stability-margin sweep for both cases).

## Open Questions (STALE as of 2026-08-17 — see Cumene section below for current state)

- ~~T_runaway for cumene still needs to be *derived*~~ — done, see below.
- ~~Final page layout~~ — done: web UI built and functional (both reactions), see below.
- ~~Whether to report boiling-point, Semenov, or both~~ — resolved: both shown, whichever binds
  first (smallest D_crit) drives a single "recommended diameter" with a 15% margin.

## Cumene — Packed-Bed Reactor (built 2026-08-17)

**User feedback that triggered this work**: after using the Aniline calculator, the user found it
"dull" — L/D only ever changed throughput, never the recommended diameter itself, because
Aniline's empty-tube heat transfer uses Nu_lam=3.66, a constant independent of velocity across the
entire practical range (deeply laminar, Re~4-375 even with the recycle loop's 3x boost). Traced and
confirmed: not a bug, just Aniline's real physics — reaching turbulent flow near the critical
diameter would need L/D in the thousands (a tube hundreds of meters long), not a real design point.
User's explicit instruction: flag it if Cumene heads toward the same dead end rather than
completing a "useless" build — check before investing further effort.

**Scoping check (before full build)**: packed-bed heat transfer (Leva's correlation, Nu_w =
2.03·Re_p^0.8·exp(-6dp/D)) is genuinely Re-dependent across its whole valid range, unlike the
laminar empty-tube constant — confirmed by direct numerical check. But the *governing* constraint
shifted: Cumene's exotherm (-110.6 kJ/mol) is much milder than Aniline's (-552.9 kJ/mol), and
packed beds cool far more effectively (U in the thousands of W/m²K vs. Aniline's tens-hundreds) —
so thermal runaway is rarely reached at practical diameters. What *does* bind is pressure drop
(Ergun equation), which rises by >1000x over the same diameter range. **The two criteria move in
opposite directions with L/D** — Semenov thermal limit rises with L/D (more velocity → better
cooling), pressure-drop limit falls with L/D (more velocity → more ΔP) — crossing near L/D≈140.
This gives the tool a genuine optimization surface (the binding/recommended diameter has a real
maximum near the crossover) that Aniline never had. User decision: show both criteria, recommend
whichever binds first (smaller D_crit) — same pattern as Aniline's dual-criteria approach,
generalized.

**Two real modeling issues caught and fixed during the build, not swept under the rug**:
1. Leva's correlation is only valid for Re_p > ~20; extrapolating below that produced a spurious
   ~400°C "runaway" at the smallest diameters (D=10mm, L/D=25, Re_p=4.98) that vanished the instant
   Re_p crossed 20 — a correlation-extrapolation artifact, not real physics. Fixed by flooring the
   value used *in* the correlation at Re_p=20 (the correlation's own boundary), while still
   reporting the true Re_p for diagnostics.
2. At D comparable to the 2mm pellet size (D/dp < ~10), Leva's own exp(-6dp/D) wall-proximity term
   correctly predicts severe channeling heat-transfer penalties — real physics, not an artifact,
   but not a regime any real packed-bed design would use (standard practice: D/dp ≥ 20 to avoid
   wall-channeling/flow maldistribution). Fixed by raising the sweep's minimum diameter to 50mm
   (D/dp=25), not by suppressing the correlation's own valid behavior.

**Kinetics/reaction data**: K1 = 6510·exp(-52564/RT) (Ben Hamida 2026 preprint, an Aspen HYSYS PFR
simulation study of the real Q-Max process), treated as an already bed-effective (lumped)
pseudo-homogeneous rate constant — no separate catalyst effectiveness factor layered on top, to
avoid double-counting whatever mass-transfer limitation the fit already implies. DIPB
transalkylation side reaction has sourced kinetics but is not yet modeled (main alkylation only).
Benzene:propylene = 0.85:1 (mid of source's 0.75-0.95:1), taken at face value. Unsourced,
reasonable-judgment assumptions: target conversion 90%, Tc = Tfeed-20K, allowable ΔP = 10% of
operating pressure (common packed-bed heuristic, not verified against a specific page in this
project's own library — flag if it needs firming up).

**Files**: `cumene_pfr_feasibility.py` (base packed-bed model, Ergun ΔP, Leva U, exact Semenov T*),
`cumene_webdata_precompute.py` (generates `cumene_calculator_data.js`, 14-point L/D grid with extra
resolution 100-160 around the crossover, D swept 50-300mm).

## Web UI — Generalized to Multi-Reaction (2026-08-17)

Data schema generalized from Aniline's original hardcoded boiling/semenov fields to a generic
`limits: [{key, label, kind, D_crit_mm, T_ref_C?}]` list per L/D, since Cumene introduced a limit
*type* (pressure-drop) Aniline never had. `reference` now also carries `tableRows`, `gaps`,
`footer`, `limitsExplanation`, and `throughputUnit` as precomputed strings/arrays, so the JS
renders generically off whatever a reaction's data provides — no reaction-specific branches.
Aniline's own data regenerated under the new schema (numbers unchanged, only the JSON shape
changed). `reactor_scaleup_calculator_block.html` and `test_harness.html` fully rewritten:
reaction dropdown now has both options enabled, two reaction-scheme SVGs (aniline/cumene) toggle
by selection, a 4th "pressure drop vs. D" chart shows only when a reaction's data includes
`dP_kPa` (Cumene only). L/D slider range now rebuilds per-reaction (Aniline: 8 points 25-200;
Cumene: 14 points 25-250). Not yet re-confirmed by user in-browser after this rewrite.
