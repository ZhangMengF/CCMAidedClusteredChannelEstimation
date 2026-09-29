# Assumed CCM-error SD mismatch: SNR = 5 dB, Kp = 32

Independent result; not used by the response letter yet.

From `OpenSourceCodes_new`, reproduce with:

```matlab
exp_AssumedVarianceCompensation(20, ...
    fullfile(pwd,'results','assumed_variance_compensation_5dB_Kp32'),4,5,32);
```

Only SNR and pilot-subcarrier count differ from the 15 dB / Kp = 64 experiment.
The 20 trials, RNG seed 260927, channel model, actual CCM-error SDs, ratio grid,
power calibration, fixed-power spread refitting, spread bounds, iteration budgets,
and channel estimator are unchanged. Each ratio uses identical channel/pilot data
within a trial; changing Kp changes the random draw dimensions, so the two operating
points should not be treated as identical channel/noise realizations.

MAT and CSV contain the numerical results. PNG and editable FIG use
`plot_AssumedVarianceCompensation.m`; no PDF figure is exported.
The original result directory and response-letter figure are not replaced.
