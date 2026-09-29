# CCM-aided channel estimation

MATLAB code and saved results for the revised manuscript **and its response
letter**. Beam-delay diagonal parameter fitting and FFT-PCG final LMMSE solves
are enabled by default. Use the selected entries below for the final results;
historical controls are retained for provenance, not silently substituted.

## Requirements and first check

- Tested with MATLAB R2025a Update 1 on Linux; other releases are not certified.
- Parallel Computing Toolbox is needed for the default parallel Fig.1,
  assumed-variance and acceleration-NMSE-recheck entries. Fig.1 and
  assumed-variance entries accept `workers=0` for serial execution.
- No external dataset, Sionna, Python, GPU, Optimization Toolbox or Signal
  Processing Toolbox is required by these experiments.
- Open MATLAB in this directory. `which Base` must point to this package,
  not another version on your MATLAB path.

```matlab
test_Numerics;
test_OptimizerSafety;
test_ResponseExperiments;
```

These are small implementation checks, not full Monte Carlo experiments.
Start with saved-data plotting if you only want to inspect results. Full runs,
especially the dense timing baseline, can be slow.

## Reproduce manuscript and response results

Run entries individually from this directory. Use a fresh output root to avoid
overwriting supplied evidence. Each experiment creates its output folder.

```matlab
out = fullfile(pwd,'results','my_reproduction');

% Main Fig.1: all panels, including accelerated online/offline fitting.
exp_CombinedResults(fullfile(out,'combined'),8);
% Use 0 instead of 8 for serial execution. Plotting is included.

% Initialization: 20 trials, 0:30 recorded, 0:15 displayed.
exp_InitializationConvergence(20,fullfile(out,'initialization'));
PlotResults(fullfile(out,'initialization'),fullfile(out,'initialization'),{'initialization'});

% Main coupling table / Response R1.3.
exp_CopulaCoupling(100,fullfile(out,'copula'));

% Response R1.1: LoS table.
exp_ResponseLoS(100,fullfile(out,'los'));

% Response R1.5: six pilot-only methods, 100/50 trials per SNR/Kp point.
exp_PilotOnlyComparison([100,50],fullfile(out,'pilot_only'));
PlotResults(fullfile(out,'pilot_only'),fullfile(out,'pilot_only'),{'pilot_only'});

% Response R2.4: 100 trials, 15 dB, Kp=64; close any existing pool first.
exp_AssumedVarianceCompensation(100,fullfile(out,'assumed_variance'),4,15,64);
plot_AssumedVarianceCompensation(fullfile(out,'assumed_variance'));

% Response R2.5: finite-subpath model, 20 trials, 1/8 fitting snapshots.
exp_ParameterAccuracy(20,fullfile(out,'accuracy'),true);
plot_ParameterAccuracy(fullfile(out,'accuracy'));

% Response R3.4: power/spread ablation table.
exp_ParameterAblation(100,fullfile(out,'ablation'));
```

The acceleration table uses separate NMSE and timing datasets. Its exact
sources and commands are in [RELEASE_REPRODUCIBILITY.md](RELEASE_REPRODUCIBILITY.md).
Default `exp_Acceleration` alone reproduces the older dataset, not every current
table entry. Run timing experiments alone, without competing simulation jobs.

**Legacy wrappers are not exact final-paper recipes:** `run_all` runs four
earlier studies and can overwrite root results. `run_ResponseExperiments` also
runs the superseded prior-mismatch control and uses 100 accuracy trials, whereas
the response uses 20. Both remain unchanged for backward compatibility.

## Redraw saved results without simulation

```matlab
plot_CombinedResults(fullfile(pwd,'results','combined_accelerated'));
PlotResults('results','results/redraw',{'initialization','pilot_only'});
plot_AssumedVarianceCompensation('results/assumed_variance_compensation_15dB_Kp64_mc100');
plot_ParameterAccuracy('results/parameter_accuracy_truncated');
```

Plot functions do not regenerate channels. Dedicated plotters write into their
input directory; copy that result directory first to preserve its PNG/FIG.
Fig.1 needs `figure_assets/Fig1_layout.fig`; retain this asset when distributing.
New figures are PNG, with editable FIG where supported; no extra PDF export.

## Reproducibility and model conventions

- Fixed Twister seeds and paired channel/noise/prior draws are used. Monte Carlo
  counts are independent trials, not historical pilot-snapshot counts.
- NMSE pools squared errors and channel energies before conversion to dB.
  MAT files retain configurations and raw statistics; CSV values are unrounded.
- Reproduction requires unchanged code, seeds, parameters and trial counts.
  MATLAB versions/thread counts may cause small floating-point differences;
  runtime and speedup depend on hardware/load, not just seeds.
- Resume only compatible checkpoints. Fig.1 checks its source snapshot; use a
  fresh folder after source changes or to force an independent complete rerun.
- Ordinary experiments use CDL macro parameters and finite weighted subpaths,
  not a full TR 38.901 implementation. The copula control instead samples its
  marginals with equal-power rays. Oracle integrates the true profile and CCM
  uncertainty, not the realized subpath covariance.
- R2.4 is a controlled **spread-refitting study with frozen powers and explicit
  spread bounds**, not an unchanged joint-fit estimator. See
  [IMPLEMENTATION_NOTES.md](IMPLEMENTATION_NOTES.md) for full details.
- Fig.1 uses DFT at every Pilot-only point. DFT wins 7/10 comparison conditions;
  the separate equal-weight linear-domain diagnostic favors linear interpolation.
  No per-point method switching is used.

## Files and verification status

| Item | Purpose |
|---|---|
| `Base.m` | Shared channel, estimator and numerical routines |
| `exp_*.m` | Independent experiments; selected entries listed above |
| `PlotResults.m`, `plot_*.m` | Saved-data plotting |
| `test_*.m` | Small numerical regression checks |
| `figure_assets/` | Frozen Fig.1 layout and historical provenance |
| [Release index](RELEASE_REPRODUCIBILITY.md) | Exact manuscript/response-to-data mapping |
| [Results index](results/README.md) | Current versus historical outputs |
| [Audit report](results/release_audit/REPORT.md) | Previous reproduction checks and their scope |
