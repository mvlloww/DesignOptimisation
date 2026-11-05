# Black-box session toolkit (Optimisation Skills 2)

Toolbox-free MATLAB (no Statistics or Optimization Toolbox needed). Tested end to end in Octave 8.4.

## Before the session
1. `rehearse` — full dry run against `fake_blackbox`. Run it on the laptop you'll take into the room.
2. Delete `bb_log.csv`, then run `make_sampling_plan` once. Check the printed phi_q, L2-star and D* values and keep them for the report.
3. Print `bb_plan.csv` as a paper backup, with blank columns for MCSI.
4. On the day, check the dial resolution and update `cfg.res` in `bb_config.m` before generating the plan if it differs.

## In the session (one person dials, one types, one reads back)
```
bb_next          % prints the point to dial and the predicted MCSI
bb_record(7.42)  % log the reading (written to CSV straight away)
```
| Calls | Phase | What to do |
|---|---|---|
| 1–25 | LHS training | `bb_next` / `bb_record` loop |
| 26 | repeat of run 1 | if not identical: box is noisy, set `cfg.rbfLambda = [1e-6 1e-3 1e-2 1e-1]` |
| 27–31 | hold-out | predicted value shown before you query (genuine out-of-sample test) |
| 32 | repeat of run 13 | checks drift as well as noise |
| — | **decision point** | `bb_status` → train-vs-hold-out table and model choice |
| 33–47 | infill | `bb_next('auto')` alternates exploit and explore, and auto-switches to explore once converged |
| 48–50 | confirm | `bb_next('confirm')`: optimum, then ±5% steps per axis (local sensitivity) |

Use `bb_next('exploit')` or `bb_next('explore')` to override. `bb_status(true)` makes the report figures.

## Files
| File | Purpose |
|---|---|
| `bb_config.m` | bounds, dial resolution, call budget, model grid |
| `make_sampling_plan.m` | best of 2000 LHS ranked on phi_q (maximin) and L2-star; exact D*; greedy maximin hold-out |
| `l2star.m`, `dstar_exact.m`, `phi_q.m` | sampling-plan metrics |
| `fit_surrogate.m`, `predict_surrogate.m` | poly (order 1–3) and RBF (Gaussian, MQ, IMQ, cubic, TPS, with linear tail) on [0,1]-normalised x and z-scored y |
| `loocv.m`, `compare_models.m` | LOOCV NRMSE and Q², plus hold-out NRMSE |
| `maximise_surrogate.m` | grid → multi-start fminsearch |
| `bb_next.m`, `bb_record.m`, `bb_status.m` | session interface |
| `bb_plots.m` | LOOCV parity plot and 2D slices through the optimum |
| `fake_blackbox.m`, `rehearse.m` | rehearsal only |
