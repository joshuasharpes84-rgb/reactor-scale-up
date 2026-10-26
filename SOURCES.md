# Sources — Jacketed / Packed-Bed Reactor Scale-Up

Where each type of information in this project (feasibility models, Cantera cross-check, web
calculator) comes from. Topic-level, not number-by-number: the aim is to show nothing was invented,
and to point a reader to the right book and section if they want to dig. Full sourcing history:
`Reactor Scaleup Tool - Plan.md`.

**Tiers**
- **Library** — reference textbook in `Downloads\ChemE\Textbooks (Level 1)\` (used across projects)
- **Project PDF** — paper or excerpt held in this folder's `Reference Materials\`
- **Web** — public page or open-access paper, URL or identifier given
- **Print only** — book cited from a paper copy; no digital copy for page lookup
- **Own work** — derived, fitted or chosen in this project; the basis is stated

Copyrighted PDFs are cited here but never published.

| Topic | Source | Tier | Location |
|---|---|---|---|
| Reaction shortlist and industrial process descriptions (DuPont/KBR aniline, UOP Q-Max cumene) | Meyers, *Handbook of Petrochemical Production Processes* | Library + Project PDF | `Textbooks (Level 1)\Process-Specific References\`; copy in `Reference Materials\` |
| Heats of formation for ΔH_rxn by Hess's law (aniline −552.9, cumene −110.6 kJ/mol) | NIST Chemistry WebBook (Pedley 1986; Hatton/Hildenbrand 1962; CODATA 1984; Roux/Temprado 2008; Prosen 1945) | Web | https://webbook.nist.gov |
| Aniline kinetics (first order in nitrobenzene, Pd catalyst, Ea 20–23 kJ/mol) | Prekob et al., *Nanomaterials* 11(5):1172 (2021), open access, PMC8146346 | Web | primary article, fetched and read |
| Cumene kinetics (pseudo-first order in propylene, Ea 52.6 kJ/mol; DIPB side reaction 55.0 kJ/mol), Q-Max operating conditions | Ben Hamida, *Preprints.org* 202601.1995.v1 (2026), citing Hazmi 2022 and Norouzi & Fatemi 2012 | Project PDF | `Reference Materials\Ben Hamida 2026 ...pdf` |
| Independent check on the cumene rate-law form | Ercan et al. (1998), zeolite-catalyzed alkylation, as cited in Davis & Davis pp. 235–236 | Library | `Textbooks (Level 1)\Core Textbooks\Reaction Kinetics Textbooks.pdf` |
| Thermal runaway criterion (Semenov / Frank-Kamenetskii) | Davis & Davis, *Fundamentals of Chemical Reaction Engineering*, Ch. 9 | Library | as above |
| Laminar tube heat transfer (Nu = 3.66) | Standard fully developed laminar-flow result | Standard correlation | e.g. Incropera, `Textbooks (Level 1)\Core Textbooks\Heat and Mass Transfer Textbook.pdf` |
| Packed-bed wall heat transfer | Leva's correlation, Nu_w = 2.03·Re_p^0.8·exp(−6d_p/D), valid Re_p > ~20 | Standard correlation | page not yet recorded |
| Packed-bed pressure drop | Ergun equation | Standard correlation | page not yet recorded |
| Fluid properties (density, c_p, viscosity, conductivity) | `thermo` Python library | Own work (generated) | feasibility scripts |
| Independent reactor-network cross-check | Cantera, with the same kinetics and formation enthalpies | Own work | `aniline_cantera_*.py`, `aniline_mech.yaml` |
| Allowable ΔP = 10 % of operating pressure; target conversion 90 %; coolant T = feed − 20 K; 15 % design margin | Common heuristics / judgment, not tied to a page | Own work | flagged in tool |
| Packed-bed minimum D/d_p ≥ 20 | Standard practice to avoid wall channeling | Own work | sweep starts at D = 50 mm (D/d_p = 25) |
