# Stability check: 100 trials, SNR = 15 dB, Kp = 64

Independent result directory; do not replace the response-letter figure automatically.

Run from `OpenSourceCodes_new`:

```matlab
exp_AssumedVarianceCompensation(100, ...
    fullfile(pwd,'results','assumed_variance_compensation_15dB_Kp64_mc100'),4,15,64);
```

The only simulation change relative to `results/assumed_variance_compensation`
is increasing the Monte Carlo count from 20 to 100. Seed 260927 and all channel,
calibration, fixed-power spread-refitting, spread-bound, iteration, and final-solver
settings remain unchanged. The common seed should retain the original first 20
realizations; verify their numerical results before interpreting differences.

Report paired bootstrap 95% confidence intervals for `Frozen NMSE - Refit NMSE`
(positive means improvement). Do not infer stability from visual smoothness alone.
MAT/CSV/PNG/FIG are saved here; no standalone PDF figure is exported.

## Completed comparison

- The first 20 trials reproduce the original energy and squared-error arrays exactly
  (maximum absolute difference zero). All 100-trial error entries are finite.
- Angular mismatch, a = 0.25: gain 0.2683 -> 0.2667 dB;
  100-trial paired 95% CI [0.2456, 0.2892] dB.
- Angular mismatch, a = 0.5: gain 0.1361 -> 0.1355 dB;
  100-trial paired 95% CI [0.1159, 0.1564] dB.
- Delay mismatch remains without an NMSE benefit: gains range from -0.0446
  to -0.0111 dB away from a = 1. Some small degradations are statistically resolved.
- Angular a = 2 and 4 no longer show tiny positive point gains: -0.0153 and
  -0.0166 dB respectively. Do not claim uniform improvement at all mismatch levels.
- Absolute NMSE values shift by up to 0.2433 dB between 20 and 100 trials;
  the main compensation trend is stable, not every numerical value.
