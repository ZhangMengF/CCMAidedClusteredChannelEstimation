# CCM-Aided Clustered Channel Estimation

MATLAB implementation of channel cluster map (CCM)-aided channel estimation
for clustered SIMO-OFDM channels. The estimator combines uncertain cluster
delay/angle priors with pilot-based estimation of cluster powers and
intra-cluster spreads.

The package provides:

- CDL-parameter-based channel generation and Monte Carlo evaluation.
- Alternating power/spread estimation using EM/SQUAREM and L-BFGS.
- Beam-delay diagonal likelihood approximation for fast parameter fitting.
- FFT-PCG LMMSE channel estimation.
- Oracle, CPM-aided, fixed-profile, and pilot-only comparison methods.
- Experiments on initialization, historical observations, delay-angle coupling,
  prior uncertainty, and parameter-estimation accuracy.

Both computational accelerations are enabled by default.

## Requirements

Tested with **MATLAB R2025a Update 1 on Linux**.

Parallel Computing Toolbox is required for the default parallel experiments.
The combined comparison and assumed-variance experiments accept a worker count
of `0` for serial execution. No external datasets, Python, Sionna, GPU,
Optimization Toolbox, or Signal Processing Toolbox are required.

Open MATLAB in this directory. Check that `which Base` points to this
repository rather than another copy on the MATLAB path.

## Quick start

Run the numerical checks:

```matlab
test_Numerics;
test_OptimizerSafety;
test_ResponseExperiments;
```

For a short initialization experiment:

```matlab
out = fullfile(pwd,'results','quick_start');
exp_InitializationConvergence(2,out);
PlotResults(out,out,{'initialization'});
```

This uses two independent trials for a quick check. Use 20 trials for the
supplied initialization configuration.

## Experiments

Run entries individually from the repository directory. Choose a fresh output
folder to avoid overwriting saved results.

| Experiment | Entry | Example configuration |
|---|---|---|
| NMSE versus SNR, pilot count, and history length | `exp_CombinedResults(folder,8)` | 100/50/100 trials per configuration; eight workers |
| Initialization convergence | `exp_InitializationConvergence(20,folder)` | 5 dB, 64 pilots, three initializations |
| Parameter-estimation acceleration | `exp_Acceleration([20,20,20,5,5],folder)` | Dense versus accelerated fitting |
| NMSE-only acceleration comparison | `exp_AccelerationNMSERecheck(100,folder,8)` | 100 paired trials per pilot count |
| Delay-angle coupling | `exp_CopulaCoupling(100,folder)` | Gaussian-copula dependence, fixed marginal profiles |
| Pilot-only interpolation | `exp_PilotOnlyComparison([100,50],folder)` | Six LS-based methods; SNR and pilot-count scans |
| LoS channel estimation | `exp_ResponseLoS(100,folder)` | CDL-D-based channel, 32/64 pilots |
| Assumed prior-error mismatch | `exp_AssumedVarianceCompensation(100,folder,4,15,64)` | 100 trials, four workers, 15 dB, 64 pilots |
| Parameter and covariance accuracy | `exp_ParameterAccuracy(20,folder,true)` | Finite subpaths, 1/8 fitting snapshots |
| Power/spread ablation | `exp_ParameterAblation(100,folder)` | Fixed, power-only, spread-only, and joint fitting |

For example, run the combined comparison with:

```matlab
folder = fullfile(pwd,'results','my_comparison');
exp_CombinedResults(folder,8);  % use 0 for serial execution
```

It saves per-configuration checkpoints and generates the three-panel figure.
Use a fresh folder after source changes; compatible completed checkpoints are
reused on restart.

The individual scan scripts `exp_SNRVsMonteCarloNMSE`,
`exp_KpVsMonteCarloNMSE`, and `exp_HistypSampNumScan` are also available.
In the first two scripts, set `estor.pilot_only_method` to `dft`, `linear`,
`quadratic`, `fir`, `pchip`, or `spline`. DFT is the default.

**Timing experiments:** run `exp_Acceleration` alone, without competing jobs.
Runtime depends on hardware and system load. The parallel NMSE-only recheck
does not measure acceleration speedup. Close an existing parallel pool before
starting `exp_AssumedVarianceCompensation`, which creates its own pool.

Batch helpers `run_all` and `run_ResponseExperiments` run predefined subsets,
not every experiment. Inspect their settings before use; their default outputs
may overwrite existing result files.

## Plot saved data

Simulation and plotting are separate. To redraw the supplied results:

```matlab
plot_CombinedResults(fullfile(pwd,'results','combined_accelerated'));
PlotResults('results','results/redraw',{'initialization','pilot_only'});
plot_AssumedVarianceCompensation('results/assumed_variance_compensation_15dB_Kp64_mc100');
plot_ParameterAccuracy('results/parameter_accuracy_truncated');
```

Dedicated plotters save into their input directory; copy the result directory
first if you want to preserve its existing figures. The combined plot uses
`figure_assets/Fig1_layout.fig`, which must remain in the package.

Outputs are MAT data, CSV summaries, PNG images, and editable FIG files where
supported. Plotting does not regenerate channels. MATLAB rendering can vary
slightly across environments even when numerical data agree.

## Reproducibility

- Experiments explicitly set Twister random seeds. Compared methods use paired
  channel, noise, and prior realizations.
- Monte Carlo counts denote independent trials, not historical pilot snapshots.
- NMSE is computed from pooled squared errors and channel energies, then
  converted to dB. CSV summaries are unrounded.
- Preserve code, parameters, seeds, and trial counts when reproducing results.
  MATLAB versions and thread counts can introduce small floating-point differences.
- Standard fitting uses AO/EM/L-BFGS iteration limits of 3/3/2. The initialization
  study records 30 AO iterations and plots iterations 0 through 15.
- The supplied combined result is
  [results/combined_accelerated/Result.png](results/combined_accelerated/Result.png).
  Other result variants are retained; see the [data index](results/README.md)
  for their configurations and provenance.

## Model conventions

The channel generator uses CDL macro parameters with finite weighted subpaths;
it is not a full implementation of 3GPP TR 38.901. The copula experiment instead
uses equal-power rays sampled from fixed delay and angular marginals.

Oracle uses true cluster powers/spreads and the continuous truncated profile
while retaining CCM prior uncertainty, rather than using realized subpath
covariances.

The assumed-variance experiment isolates spread adaptation: it freezes calibrated
powers and constrains fitted spreads. It is a controlled variant of joint
power/spread fitting. The ordinary estimator does not enable those extra bounds.

## Code organization

| Files | Purpose |
|---|---|
| `Base.m` | Shared channel generation, estimators, and numerical routines |
| `exp_*.m` | Experiment entry points |
| `PlotResults.m`, `plot_*.m` | Plotting from saved data |
| `test_*.m` | Numerical and regression checks |
| `SaveResult.m`, `ResponseStatistics.m` | Result export and statistical summaries |
| `GaussianCovarianceTools.m` | Structured covariance diagnostics and Gaussian controls |
| `figure_assets/` | Reusable plotting assets |
| `results/` | Saved configurations, numerical results, and figures |

Detailed algorithm conventions and additional controls are documented in
[IMPLEMENTATION_NOTES.md](IMPLEMENTATION_NOTES.md).
