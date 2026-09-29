# Submission data and reproduction map

This is the release-facing index. Exploratory files/results are retained for
provenance, not silently deleted or substituted for the manuscript results.
The manuscript and response-letter TEX files are not written by any entry below.

## Figure 1: unified accelerated run

Run `exp_CombinedResults` (Parallel Computing Toolbox; eight workers by default).
For a serial run use `exp_CombinedResults(folder,0)`. It calls the same
`Base.CalcuNMSEsByMonteCarlo` implementation as the existing three scan scripts:

| Panel | Grid | Independent trials per point | RNG seeds |
|---|---|---:|---|
| a | SNR = -5,0,5,10,15 dB; Kp=32 | 100 | 101:105 |
| b | Kp = 8,16,32,64,128; SNR=5 dB | 50 | 101:105 |
| c | S=0,1,3,7,15; blockage=.5,.3,0; Kp=64, SNR=5 dB | 100 | 1000*block-index + 100*S-index |

S=0 selects single-current-snapshot online fitting; S>0 fits using historical
snapshots. In panel c, S=0 is displayed as the horizontal online reference,
not connected to the historical curve. All fitting uses `beam_delay_diag`,
all final solves use `fft_pcg`, and AO/EM/L-BFGS limits remain 3/3/2.
All five estimators in a/b share trial data. Panel c reuses the existing
independent historical-channel generation and blockage convention unchanged.

Default output is `results/combined_accelerated/`. Per-configuration MAT
checkpoints contain raw errors/energies, MC counts and seeds; `source.mat`
records Base and the experiment entry. Restarting the same command resumes
completed configurations only when the source matches. Workers use two threads;
this is not a runtime-benchmark experiment.

Run `plot_CombinedResults(folder)` to redraw without simulation. It uses
`figure_assets/Fig1_layout.fig`, a frozen copy of the user-approved three-panel
layout, and changes only curve YData. It checks the original XData, axes positions,
limits, ticks and fonts. Legends, titles, colors, styles and line widths are retained.
The only user-approved layout exception is panel c's upper Y limit: -5.5 to
-5 dB, to include the new -5.0236 dB point; its tick spacing is unchanged.
The old `Fig1_original.fig` is retained as historical provenance, not as a numerical
source for the new result.

## Other manuscript and response results

| Target | Saved numerical source | Experiment entry / plotting entry |
|---|---|---|
| Main Fig.2 initialization | `results/initialization.mat` | `exp_InitializationConvergence(20,folder)` / `PlotResults(folder,folder,{'initialization'})` |
| Main Table I acceleration | Separate NMSE/timing sources listed below | `exp_AccelerationNMSERecheck`, `exp_Acceleration` (timing must run alone) |
| Main Table II / R1.3 coupling | `results/copula.mat` | `exp_CopulaCoupling(100,folder)` |
| R1.1 LoS table | `results/response_supplement_20260925/response_los.mat` | `exp_ResponseLoS(100,folder)` |
| R1.5 pilot-only figure | `results/pilot_only.mat` | `exp_PilotOnlyComparison([100,50],folder)` / `PlotResults(folder,folder,{'pilot_only'})` |
| R2.4 assumed-error mismatch | `results/assumed_variance_compensation_15dB_Kp64_mc100/assumed_variance_compensation.mat` | `exp_AssumedVarianceCompensation(100,folder,4,15,64)` / `plot_AssumedVarianceCompensation(folder)` |
| R2.5 accuracy figure/table | `results/parameter_accuracy_truncated/parameter_accuracy.mat` | `exp_ParameterAccuracy(20,folder,true)` / `plot_ParameterAccuracy(folder)` |
| R3.4 power/spread ablation | `results/response_supplement_20260925/parameter_ablation.mat` | `exp_ParameterAblation(100,folder)` |

Use fresh folders for verification. Do not overwrite historical data to hide
differences. Saved-data agreement, complete numerical reruns, representative-trial
reruns, and exact graphical reproduction are separate verification levels.

## Current acceleration-table sources

The current manuscript was re-read during release preparation. Its numerical
entries now agree, after rounding, with the following sources; the earlier
audit's mixed-NMSE warning describes an older manuscript state.

| Entries | Source | Sampling |
|---|---|---|
| NMSE and NMSE difference, Kp=8,16,32,64 | `results/acceleration_nmse_recheck_20260929/acceleration_nmse_recheck.csv` | 100 paired trials each |
| NMSE, difference and speedup, Kp=128 | `results/acceleration_kp128_timing_20260929/acceleration.csv` | 100 paired trials |
| Speedup, Kp=8,16,32,64 | Root `results/acceleration.csv` | 20,20,20,5 paired trials respectively |

Reproduce the three source groups independently, in fresh directories:

```matlab
exp_AccelerationNMSERecheck(100,fullfile(out,'acceleration_nmse'),8,[16,64,8,32]);
% Timing: run alone; do not overlap with the parallel NMSE run.
exp_Acceleration(100,fullfile(out,'acceleration_128'),false,128,26092900);
exp_Acceleration([20,20,20,5],fullfile(out,'acceleration_8_64'),false,[8,16,32,64]);
```

Here `out` is the fresh root from README. Runtime speedups are not expected to
match exactly across machines or loads. The old small-sample timing sources
predate the optimizer safety update; keep their provenance distinct from fresh
current-code measurements. Concurrent-workload intervals of the saved Kp=128
run are not isolated timing measurements. A separate 50-trial timing recheck
directory is retained, but is not the source of the current table values.

## Verification scope and retained audit caveats (2026-09-29)

- Current table values match the multi-source breakdown above, not one output
  file. No manuscript TEX was modified during release preparation.
- Initialization/coupling/old acceleration were saved before the 2026-09-25
  optimizer safety fix. The current-code full initialization replay and
  per-configuration first-trial coupling replay exactly match saved raw errors.
  This does not validate old runtime measurements for the updated implementation.
- The LoS and ablation Base snapshots differ from current Base only by optional
  spread-bound handling, disabled for those experiments. The accuracy and
  100-trial assumed-variance study contain an exact current Base snapshot.
- R1.2's illustrative values (-5.79/-6.49 dB for Fixed-profile/Proposed)
  agree with the unified new Fig.1 at 5 dB, Kp=32 after rounding.
- The legacy README's equal-weight pilot-only selection rule picks linear;
  the paper chooses DFT because it is best at 7/10 tested points. DFT is explicitly
  selected by the new Fig.1 experiment; no per-point method cherry-picking is used.

Actual rerun checks are saved under `results/release_audit/`. The earlier audit
is recorded in `results/release_audit/REPORT.md`; its table-status snapshot
predates the current table reconciliation documented above.
