# Engineering Decisions — Jacketed / Packed-Bed Reactor Scale-Up

The calls that shaped the scale-up tool: what was decided, what was rejected, and why. Sources for
each topic are in `SOURCES.md`; the full build history is `Reactor Scaleup Tool - Plan.md`.

*Built by Joshua Sharpes, EIT, with AI assistance (Claude) for drafting, coding and literature
checks. Scope, method selection, validation criteria and the corrections below are his.*

## The central finding

**There is no single scale-up limit.** The diameter at which a tubular reactor stops scaling depends
on which physical constraint binds first, and that differs by reaction:
- **Aniline (empty tube, jacketed):** thermal. Heat generation grows with volume, removal with jacket
  area, and the reactor stays deeply laminar, so the critical diameter doesn't change with L/D.
- **Cumene (packed bed):** usually pressure drop. The thermal limit rises with L/D while the
  pressure-drop limit falls, and they cross near L/D ≈ 140, giving a real optimum.

The tool shows every applicable limit and recommends from whichever binds first, with a 15 % margin.
Limits are never blended into one number.

## Data sourcing discipline

**Derive rather than trust secondary numbers.** Heats of reaction are derived by Hess's law from
NIST formation enthalpies. Runaway temperatures are derived from the kinetics (Semenov /
Frank-Kamenetskii), because tabulated values are generally proprietary.

**Primary articles only.** Kinetics were taken from papers read in full. Several activation
energies found only in search summaries (35–40 kJ/mol for aniline) were discarded because they
couldn't be traced to a primary source.

**Stated gap to industry.** The best public aniline kinetics are a dilute lab system, not DuPont's
proprietary catalyst, so throughputs are lab-kinetics scale. This is stated, not hidden.

## Decisions along the way

**Aniline's L/D insensitivity is physics, not a bug.** Joshua found the first calculator "dull":
L/D moved throughput but never the recommended diameter. It was traced to the laminar heat-transfer
constant; reaching turbulent flow would need a tube hundreds of meters long.

**Scope check before building the second case.** Joshua's instruction was to stop if cumene hit the
same dead end. It was checked first: packed-bed heat transfer does depend on velocity, so cumene adds
a genuinely different behavior and was worth building.

**Fix the real bottleneck.** The once-through aniline limit (~55 mm) was too small to be a useful
scale-up answer. A boiling jacket was tried and barely helped, because the process side dominates
heat transfer. A recycle loop through an external cooler was the effective fix (limit moves past
200 mm at recycle ratio 2).

**Correlation limits respected, not exploited.** Two artifacts were caught during the cumene build:
- Leva's correlation extrapolated below its valid range produced a false runaway at the smallest
  diameters; it is now held at its own boundary while the true value is still reported.
- Diameters near the pellet size show real but unrealistic wall channeling; the sweep starts at
  standard practice (D/d_p ≥ 20) instead of suppressing the correlation.

**No double counting.** The cumene rate constant is already a bed-effective lumped value, so no
separate catalyst effectiveness factor is applied on top.

## Validation

- An independent Cantera reactor model reproduces the aniline trend. Its critical diameter is about
  37.5 mm vs 54.6 mm from the main model; the ~30 % gap comes from Cantera's temperature- and
  composition-dependent properties. Reported as a range.
- The cumene rate-law form matches an independent study of the same reaction family.
- Recycle results confirmed by forward and backward continuation.
