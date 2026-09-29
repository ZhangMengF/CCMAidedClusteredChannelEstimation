# CCM variance-mismatch compensation experiment

## Protocol

- 100 paired trials, SNR 5 dB, Kp=32, N=32, K=256. Same seed (260925),
  finite-subpath channels, noise and standardized prior errors as the earlier
  CCM-prior mismatch experiment. No blockage.
- Scan a=[0.5,1,2,4], the actual/assumed error **SD** ratio; variance ratio=a².
  Estimator-assumed delay/angle SDs remain 3 ns and 3 degrees. Vary one actual
  SD at a time. Retain the original zero-error first-delay anchor.
- Proposed refits all cluster powers and both spread vectors at each point
  using the current single pilot. The control freezes each trial's fitted
  powers/spreads from a=1, but uses the CURRENT a-dependent prior centers.
- Oracle uses actual uncertainty statistics, true powers and truncated
  profiles, but the same noisy centers, not perfect centers.
- Original accelerated AO/SQUAREM/L-BFGS limits 3/3/2 and FFT-PCG are unchanged.
  No estimator bounds, regularizers, clipping, or outlier removal were added.
- Right axes: equally weighted mean across all 23 clusters and 100 trials,
  corresponding delay (ns) or angle (degrees), with 95% bootstrap intervals.
  Power-weighted means are saved as additional diagnostics, not substituted.

## NMSE results

All entries below are in dB. Gain is Frozen minus Proposed (positive is better).

| Scanned error | a | Proposed | Oracle | Frozen at a=1 | Refit gain | 95% paired gain interval |
|---|---:|---:|---:|---:|---:|---|
| Delay | 0.5 | -6.440 | -6.854 | -6.337 | 0.103 | [0.046, 0.177] |
| Delay | 1 | -6.253 | -6.653 | -6.253 | 0 | [0, 0] |
| Delay | 2 | -5.791 | -6.216 | -5.917 | -0.126 | [-0.207, -0.050] |
| Delay | 4 | -4.923 | -5.575 | -5.099 | -0.176 | [-0.295, -0.059] |
| Angle | 0.5 | -6.348 | -6.732 | -6.287 | 0.062 | [0.023, 0.108] |
| Angle | 1 | -6.253 | -6.653 | -6.253 | 0 | [0, 0] |
| Angle | 2 | -6.087 | -6.409 | -6.110 | -0.023 | [-0.090, 0.064] |
| Angle | 4 | -4.827 | -5.762 | -5.386 | -0.560 | [-0.961, -0.228] |

These results do NOT establish beneficial spread compensation for underestimated
prior uncertainty (a>1). The angular a=2 difference is inconclusive; at a=4,
refitting is significantly worse than freezing in both scans. This control
tests JOINT power/spread adaptation, not the causal effect of spreads alone.

## Extreme spread estimates: retain, do not hide

The requested arithmetic means are dominated by extreme estimates:

| a | Mean delay spread in delay scan (ns) | Mean angular spread in angle scan (deg) |
|---:|---:|---:|
| 0.5 | 9.441 | 9.899 |
| 1 | 32.833 | 10.504 |
| 2 | 23.103 | 14.241 |
| 4 | 1.364e28 | 2.036e34 |

This is not a plotting error. `inspect_CCMVarianceCompensation` replays four
trials without changing the optimizer and reproduces their fitted parameters.
Examples: trial 42 (delay a=4) reaches a delay factor of 3.137e22 seconds;
trial 64 (angle a=4) reaches an angular factor of 8.174e35 radians. Their
recorded AO objective decreases throughout. Trial 97 already has an extreme
delay factor at a=1 (50.85 microseconds).

Code inspection shows unconstrained log-spread optimization without an upper
bound. Objective descent therefore does not ensure physically meaningful
spread values. Large spreads can approach limiting covariance shapes, so
extreme parameter magnitudes need not imply an equally extreme channel NMSE.
Parameter bounds/regularization would be a separate algorithm change; none
was introduced in this experiment.

For context only, medians across the same cluster/trial entries are
9.356, 9.559, 10.133, 10.979 ns in the delay scan and
9.332, 9.574, 9.948, 10.423 degrees in the angle scan. These weak trends do not
override the arithmetic means or establish an NMSE benefit. The requested
plot retains the actual means and all trials.

## Files and reproduction

- `variance_compensation.png` and visible editable `.fig`: requested dual-axis plot.
- `variance_compensation.mat`: all fitted parameters, anchor parameters,
  per-trial errors/energies, configurations, source snapshots and bootstrap draws.
- `variance_compensation.csv`: unrounded central statistics.
- `extreme_spread_diagnostics.mat` and `diagnostic.log`: original AO traces
  for trial/axis/a cases [97,1,1], [42,1,4], [60,2,4], [64,2,4].
- Run `exp_CCMVarianceCompensation(100,folder,8)`, then
  `inspect_CCMVarianceCompensation(folder)`; redraw without simulation using
  `plot_CCMVarianceCompensation(folder)`. Default worker count is 4; 0 is serial.
- `parallel_run.log` records the complete 100-trial execution with 8 process
  workers and 2 computational threads each. `run.log` is an interrupted serial
  attempt, not the reported result. Smoke subfolders contain 2-trial checks only.

Validation: serial/parallel smoke outputs are exactly equal; a=1 Proposed and
Frozen errors are exactly equal; both matched-axis anchors agree. Full-run
Proposed/Oracle per-trial errors reproduce the previous mismatch study with
relative differences below 8e-10 (roundoff from scaled prior reconstruction).
All stored errors and fitted values are finite, despite the extreme magnitudes.
No paper or response-letter files were changed, and no PDF figure was exported.
