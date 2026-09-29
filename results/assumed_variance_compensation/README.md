# Quick fixed-observation variance-compensation study

20 paired trials, 15 dB, Kp=64, seed 260927. This is exploratory evidence,
not a replacement for the paper's original Proposed results.

## Controlled changes

- True channel, pilot observations, and noisy CCM centers remain identical
  across all points. Actual delay/angle error SDs stay at 3 ns / 3 degrees.
- a=actual/assumed SD, so a=[0.25,0.5,1,2,4] corresponds to assumed SDs
  [12,6,3,1.5,0.75] ns or degrees. Variance ratio is a². Change one axis at a time.
- A matched-uncertainty joint fit provides one common power vector per trial.
  All scan fits hold this power fixed and optimize both spreads, starting
  from the same reference factors with the same AO/SQUAREM/L-BFGS budget 3/3/2.
  The matched spread-only fit defines the frozen-profile control.
- Experimental change: log-spread optimization rejects steps outside fixed
  [0.05,4] reference-factor bounds (delay [0.55,44] ns, angle [0.5,40] degrees).
  No post-hoc clipping or trial exclusion. Base's default remains unbounded.
- Original finite-subpath generator, diagonal likelihood and FFT-PCG are retained.
  Oracle knows the actual statistics and true truncated profiles; it uses the
  same noisy centers and is computed once per trial, hence exactly invariant.

## Findings

Oracle NMSE is -13.68485 dB throughout. At a=1 both fitted/frozen controls
give -12.90787 dB. Refit gains below are Frozen minus Spread-refit in dB:

| a | Delay-scan gain | Angular-scan gain | Angular gain 95% paired interval | Mean angular spread (deg) |
|---:|---:|---:|---|---:|
| 0.25 | -0.070 | 0.268 | [0.217,0.325] | 7.112 |
| 0.5 | -0.015 | 0.136 | [0.089,0.184] | 8.161 |
| 1 | 0 | 0 | [0,0] | 9.329 |
| 2 | -0.030 | 0.003 | [-0.029,0.055] | 9.701 |
| 4 | -0.055 | 0.007 | [-0.042,0.082] | 9.779 |

In this small study, overestimating angular-prior variance causes the fitted
angular spread to shrink, accompanied by a measurable benefit over frozen
spreads. This supports PARTIAL compensation in that regime. Underestimation
has no statistically clear benefit. Delay spreads shift too, but delay-scan
gains are nonpositive; at a=0.25 the interval is [-0.155,-0.004] dB, while
the other nontrivial delay intervals include zero. Do not claim universal
variance-mismatch robustness or that parameter shifts always help NMSE.

## Validation and reproduction

- 20-trial first two results exactly reproduce the 2-trial smoke run.
- Matched frozen/refit errors are identical; Oracle is identical at all points.
- Powers remain exactly fixed throughout all spread-only AO histories.
- Recorded AO stage objectives are non-increasing; all fitted values are finite.
- No final parameter lies within 1% of an upper or lower bound. Feasible line
  search may still have affected intermediate steps; this does not establish
  equivalence to the unconstrained algorithm.
- Observed delay range: 0.565 to 34.105 ns; angular range: 0.638 to 32.177 deg.
- Existing `test_ResponseExperiments` passed after the optional-bound extension.
- Run `exp_AssumedVarianceCompensation(20,folder,8)` and redraw using
  `plot_AssumedVarianceCompensation(folder)`. Default trial count is 20.
- MAT includes raw errors, full fitted vectors, powers, AO histories, captured
  optimizer output, seeds/configuration and source snapshots. CSV is unrounded;
  PNG and visible FIG are saved, without PDF images.
- No main.tex or response-letter change was made.
