# Completed response-letter experiments

All four studies are complete: **100 paired Monte Carlo trials per
configuration**, MATLAB R2025a Update 1, eight computational threads, SNR
5 dB, N=32, K=256. All ordinary fits use the current safeguarded optimizer
with AO/SQUAREM/L-BFGS maxima 3/3/2 and FFT-PCG final solving.
No main.tex or original-release file was modified.

## Results

### LoS: response_los.mat / response_los.csv

Existing CDL-D-based branch; a random-phase specular component carries 88.8%
of the power. Same finite-subpath generator and legacy angle convention.

| Kp | Pilot-only DFT | CPM | Fixed-profile | Proposed | Oracle |
|---|---:|---:|---:|---:|---:|
| 32 | -3.75 | -11.07 | -12.14 | -12.45 | -13.45 |
| 64 | -4.60 | -12.00 | -13.94 | -15.26 | -15.64 |

Values are NMSE in dB. Proposed gains over Fixed-profile are 0.31 and
1.32 dB, with paired 95% bootstrap intervals [0.09,0.52] and [1.17,1.47].
Plot: [los_scan/scan_Kp.png](los_scan/scan_Kp.png).

### Prior uncertainty: prior_mismatch.mat / prior_mismatch.csv

Estimator-assumed SDs stay at 3 ns / 3 degrees. Actual SD is multiplied
by the listed factor, one dimension at a time. Channels/noise and
standardized prior errors are shared across the scan. Oracle knows the
actual uncertainty and profile statistics, not the realized true centers.

| Scanned SD | Actual / assumed | Fixed-profile | Proposed | Oracle | Proposed gain over Fixed |
|---|---:|---:|---:|---:|---:|
| Delay | 0.5 | -5.70 | -6.44 | -6.85 | 0.74 |
| Delay | 1 | -5.60 | -6.25 | -6.65 | 0.66 |
| Delay | 2 | -5.30 | -5.79 | -6.22 | 0.49 |
| Delay | 4 | -4.76 | -4.92 | -5.57 | 0.16 |
| Angle | 0.5 | -5.63 | -6.35 | -6.73 | 0.72 |
| Angle | 1 | -5.60 | -6.25 | -6.65 | 0.66 |
| Angle | 2 | -5.47 | -6.09 | -6.41 | 0.62 |
| Angle | 4 | -4.99 | -4.83 | -5.76 | -0.16 |

All values/gains are dB. At fourfold angular error the mean ranking reverses;
the gain interval [-0.58,0.21] includes zero. Do not claim superiority there,
or a statistically established inferiority. Pilot-only is exactly invariant
at -1.3658036348 dB and both matched anchors are identical.
Plot: [prior_mismatch.png](prior_mismatch.png).

### Block ablation: parameter_ablation.mat / parameter_ablation.csv

| Kp | Fixed-profile | Power only | Spread only | Proposed | Oracle |
|---|---:|---:|---:|---:|---:|
| 32 | -5.69 | -6.29 | -5.57 | -6.45 | -6.82 |
| 64 | -8.35 | -9.16 | -8.63 | -9.41 | -9.65 |

Joint fitting gains 0.76/1.05 dB over Fixed-profile and 0.15/0.25 dB over
Power-only. Their paired gain intervals are [0.70,0.82]/[1.02,1.09] and
[0.12,0.20]/[0.23,0.27], respectively. Spread-only is 0.12 dB worse at
Kp=32 (gain interval [-0.20,-0.04]); this reversal is retained.
The frozen blocks were verified exactly unchanged in every trial.
Plot: [parameter_ablation.png](parameter_ablation.png).

### Accuracy: parameter_accuracy.mat / parameter_accuracy.csv

**Different channel law:** full untruncated model-matched complex-Gaussian
covariance, not finite rays or a diagonal-covariance generator. Fixed prior
means with uncertainty already integrated in the covariance. Independent
training snapshots and an independent test channel/pilot. S excludes the
test pilot. These are not the manuscript's same-snapshot online results.

| Kp | Fitting snapshots S | Fitted relative power error | Fitted relative covariance error | Initialized test NMSE | Proposed test NMSE | Oracle test NMSE |
|---|---:|---:|---:|---:|---:|---:|
| 32 | 1 | 0.333 | 0.327 | -5.45 | -5.64 | -5.88 |
| 32 | 8 | 0.162 | 0.148 | -5.48 | -5.78 | -5.88 |
| 64 | 1 | 0.288 | 0.272 | -7.93 | -8.15 | -8.34 |
| 64 | 8 | 0.143 | 0.106 | -7.95 | -8.29 | -8.34 |

Errors are mean norm ratios, not squared errors or dB; test NMSE is in dB.
One-snapshot fitting worsens individual-spread log-RMSE relative to the
already close reference initialization, despite substantially improving
covariance accuracy. Eight snapshots improve all reported fitted mean
errors. This is not an identifiability, consistency, or global-optimality proof.

Plots: [parameter_accuracy.png](parameter_accuracy.png),
[parameter_accuracy_nmse.png](parameter_accuracy_nmse.png).

## Reproduction and validation

- Run run_ResponseExperiments from the package, or call the four experiment
  entries separately. The package README documents their interfaces.
- Three independent data-only plotting scripts redraw the new result types.
  The LoS study reuses the original Kp scan plot.
- MAT files include raw per-trial errors and energies, configurations, seeds,
  source snapshots and 2000 paired bootstrap draws. Ablation/accuracy files
  additionally retain fitted parameters and AO histories. CSV values are
  unrounded. PNG and visible editable FIG are supplied; no PDF images.
- Three run logs document the completed full runs. Runs were parallelized by
  study, not combined into a new timing benchmark.
- Existing numerical/optimizer tests and test_ResponseExperiments passed.
  FFT matvec relative error was 5.75e-16; FFT-PCG/Cholesky estimate error
  3.36e-9 in the numerical check. The small Gaussian sampling check used 5000
  draws and had relative covariance error 0.0247.
- The pre-extension matched two-trial baseline was exactly unchanged. All
  four full studies reproduce their two-trial smoke prefixes exactly.
  The smoke subfolder contains only validation data, not reported results.
- Existing paper data predate the current safeguarded implementation and
  use different seeds; they were not overwritten or silently pooled with
  these new studies.
- Formal interpretations and channel/algorithm differences were added to
  the response letter. Its internal manuscript-edit checklist remains pending.
