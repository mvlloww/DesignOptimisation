# ME52037 Design Optimisation – Coursework Code

MATLAB code for the two coursework tasks. Tested in MATLAB R2025b.

## Task 1 – Wind turbine blade structure

The task: find the least ABS material for a 3D-printed NACA 0012 blade (chord 0.15 m, span 1.5 m) under the worst-case load, with a yield stress of 30 MPa and a safety factor of 1.5. Bending is solved by the supplied `beamBending.p`.

| File | Purpose |
|---|---|
| `RunBeamBending.m` | Analyse one design: thin skin, solid spar, hollow box spar or fully solid blade (set with `structure`). |
| `ParetoSkin.m` | Volume vs tip-deflection Pareto front for the skin (ε-constraint, `fmincon`): linear taper, then printable stepped skin in whole 0.4 mm perimeters. |
| `ParetoSpar.m` | Fair comparison of solid spar, hollow box spar and skin fronts (multistart `fmincon`). |
| `report/` | LaTeX write-up of how the code was developed. |

Approach: the blade is statically determinate, so a single `beamBending` call gives the bending moment for every design. Tip deflection is then computed with the unit-load method, which matches the solver to within 0.02%.

## Task 2 – Black-box optimisation

The task: maximise MCSI (0–10) over three knobs: x1 advert 0–200 £k, x2 targeting 0–60 £k, x3 skew 0–1. The budget is 50 calls.

| Folder | Purpose |
|---|---|
| `bb session toolboxes/` | **Main version.** One-file workflow using the Statistics, Optimization and (optionally) Global Optimization toolboxes. |
| `bb session no toolboxes/` | One-file version that needs no toolboxes. |
| `blackbox_toolkit/` | Earlier multi-file, toolbox-free version (see its own README). |

Workflow in `bb_session.m`:
1. **Fixed block (30 calls):** 25 maximin LHS training points, screened on minimum distance and L2-star discrepancy, plus 5 hold-out points.
2. **Model choice (every call):** polynomials, five RBF kernels and Kriging are compared by LOOCV.
3. **Adaptive calls:** alternate between exploit and explore (top-25% region or anywhere), then local refinement once the decision checks pass.
4. **Confirm calls:** ±5% sensitivity steps, then the final surrogate optimum is measured.
5. **Recommendation:** the best *measured* point.

```matlab
bb_session('plan')       % build the sampling plan (before the session)
bb_session               % run / resume; type each MCSI reading at the prompt
bb_session('status')     % model table, recommendation, report figures
bb_session('selftest')   % check toolboxes and fast-LOOCV formulas
```

Studies, all on separate files so the session log is never touched:

| Command | What it does |
|---|---|
| `bb_session('rehearse','hard')` | Dry run on a fake box. **This one deletes `bb_log.csv`.** |
| `bb_session('study')` | Sweeps the initial sample size `nTrain` on easy and deceptive fake boxes. |
| `bb_session('plancompare')` | Compares the chosen LHS plan with Halton and random plans. |
| `bb_session('maxcheck')` | Checks the surrogate maximisation against grid size and seed count. |
| `bb_session('truthplot')` | Plots the true fake response against the final surrogate. |

Generated outputs (`.csv`, `.mat`, `.png`) are not tracked.
