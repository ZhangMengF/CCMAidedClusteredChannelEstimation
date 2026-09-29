# Saved-result index

See [../README.md](../README.md) for executable commands. Retain MAT/CSV evidence;
do not infer the selected dataset from the newest filename or a default wrapper.
Plotting from saved data is separate from rerunning simulations.

## Selected manuscript / response evidence

| Target | Data (relative to this directory) | Image / table |
|---|---|---|
| Main Fig.1(a–c) | `combined_accelerated/` | `combined_accelerated/Result.png` |
| Initialization | `initialization.mat`, `initialization.csv` | Redraw with `PlotResults`; current-style audit image: `release_audit/plots/initialization.png` |
| Coupling table, R1.3 | `copula.mat`, `copula.csv` | CSV includes Oracle and method gaps |
| Acceleration table | See source breakdown in release index | Do not use root `acceleration.csv` alone |
| R1.1 LoS | `response_supplement_20260925/response_los.mat` and `.csv` | Table values |
| R1.5 Pilot-only | `pilot_only.mat`, `pilot_only.csv` | `pilot_only.png` |
| R2.4 prior-SD mismatch | `assumed_variance_compensation_15dB_Kp64_mc100/` | `assumed_variance_compensation.png` within that directory |
| R2.5 accuracy | `parameter_accuracy_truncated/` | `parameter_accuracy.png` and `.csv` within that directory |
| R3.4 ablation | `response_supplement_20260925/parameter_ablation.mat` and `.csv` | Table values |

## Historical and diagnostic material retained

- Root `Result.png` / `Combined_Plots.fig`: earlier mixed-source Fig.1, superseded
  by `combined_accelerated/`. Root SNR/Kp scans are earlier independent runs.
- `response_supplement_20260925/` also contains superseded prior-mismatch and
  older accuracy outputs. Only LoS and ablation are selected from that folder.
- `assumed_variance_compensation/`, `assumed_variance_compensation_5dB_Kp32/`,
  `ccm_variance_compensation/`: earlier/sample-size/configuration controls, not
  the selected 100-trial R2.4 dataset.
- `initialization_style_preview/`, root `initialization.png`: older rendering
  variants. The numerical source remains root `initialization.mat`.
- `acceleration_20trial_partial.mat`, acceleration logs, regression tests and
  `release_audit/`: provenance and verification, not additional paper curves.
- `acceleration_kp8_64_timing50_20260929/`: separate timing recheck; partial rows
  must not be treated as complete or silently substituted into the table.
- Historical PDF images are retained unchanged; current scripts do not create
  additional PDFs. Logs are retained because they document experiment provenance.

The earlier narrative is retained as [HISTORICAL_NOTES.md](HISTORICAL_NOTES.md).
No numerical dataset was deleted, moved or overwritten during documentation cleanup.
