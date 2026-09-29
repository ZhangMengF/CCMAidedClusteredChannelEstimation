# Implementation notes and development history

This document preserves detailed implementation notes and historical workflows.
For final manuscript/response commands and selected datasets, start with
[README.md](README.md) and [RELEASE_REPRODUCIBILITY.md](RELEASE_REPRODUCIBILITY.md).
Historical configurations below are not all used in the submitted revision.

**Submission release index:** see [RELEASE_REPRODUCIBILITY.md](RELEASE_REPRODUCIBILITY.md)
for the exact manuscript/response data mapping, unified accelerated Fig.1 entry,
and historical-data audit caveats. `exp_CombinedResults` runs the new
Fig.1; `plot_CombinedResults` changes only data in the approved frozen layout.
The older experiment descriptions below also include exploratory controls and
must not all be interpreted as the final submitted configuration.

This self-contained MATLAB package starts from the six original release
files. It adds the paper's beam-delay diagonal likelihood approximation,
FFT-PCG final LMMSE solves, and four reproducible experiments. Both
accelerations are **enabled by default**. Tested with MATLAB R2025a Update 1;
no Python, external datasets, optimization toolbox, or signal-processing
toolbox is required.

## Run

Open MATLAB in this directory and run `run_all`. This runs numerical checks
and the four added experiments, then writes figures and tables to `results/`.
For independent runs:

| MATLAB entry | Paper target | Default configuration |
|---|---|---|
| `exp_Acceleration` | `tab:Acceleration` | 5 dB; Kp = 8,16,32,64,128; respectively 20,20,20,5,5 paired trials |
| `exp_InitializationConvergence` | `fig:InitializationConvergence` | 5 dB; Kp = 64; 20 paired trials; 30 AO iterations |
| `exp_CopulaCoupling` | `tab:CopulaCouplingGap` | 5 dB; Kp = 32; r = 0:0.2:0.8; 100 paired trials |
| `exp_PilotOnlyComparison` | Pilot-only selection evidence | SNR = -5:5:15 (100 trials); Kp = 8,16,32,64,128 (50 trials) |

Each entry accepts an optional trial count and output directory, e.g.
`exp_InitializationConvergence(1,'results/smoke')`; the pilot-only entry takes
two counts, e.g. `exp_PilotOnlyComparison([2,2],'results/smoke')`.
`PlotResults` redraws the four experiments from saved MAT files, without
rerunning simulations. To preview only convergence in a separate directory, use
`PlotResults('results','results/initialization_style_preview',{'initialization'})`.
The initialization plot displays iterations 0--15, with both independent
vertical axes using the full plotting height; colors and markers distinguish
initializations, while solid/dashed lines distinguish NLL/NMSE. Raw
30-iteration data are retained. Figures are exported as PNG, not PDF.
The acceleration entry accepts a scalar count or one count per Kp and saves
after every trial. To continue its checkpoint, use
`exp_Acceleration([20,20,20,5,5],'results',true)`. `completed_trials` records
available draws; each summary uses only the first `mc_num(j)` draws, even if
more were retained from an earlier run. Resume only with the same model/code.
Optional fourth and fifth arguments select pilot counts and a new seed bank:
`exp_Acceleration(100,'results/acceleration_kp128_timing_20260929',false,128,26092900)`.
This runs Kp=128 alone, sequentially with eight computational threads, and
records both paired NMSE and parameter-estimation timing. The explicit seed
rule is `seed_base + 1000*Kp + trial`, matching the parallel NMSE recheck.
Use the same arguments with `true` to resume. Omitting these options preserves
the original grid and random draws. Run timing experiments without competing
simulation jobs; use a detached tmux session for long runs.

The original three entries remain available: `exp_SNRVsMonteCarloNMSE`,
`exp_KpVsMonteCarloNMSE`, and `exp_HistypSampNumScan`. The first two now include
Fixed-profile and Pilot-Only curves. Set `estor.pilot_only_method` near the top
of either entry to `dft` (default), `linear`, `quadratic`, `fir`, `pchip`, or
`spline`; keywords are case-insensitive. The figure legend and saved estimator
configuration identify the selected method. This user-selected DFT baseline
does not change the six-method comparison or its reported aggregate metric.
Their scan helpers save MAT,
CSV, FIG, and PNG files to this package's `results/`, never to the caller's
working directory. Original dense algorithms remain available by setting
`estor.parameter_backend = 'dense'` and `estor.final_solver = 'dense_chol'`.

## File organization

`Base.m` contains the shared computational functions, including trial generation,
pilot-only interpolation, and safeguarded line search. Experiment entries and
tests access these through the returned `funcs` structure. The `exp_*.m` entries,
two scan helpers, `PlotResults.m`, `SaveResult.m`, and tests remain separate.

## Statistical and timing conventions

For an independent parallel NMSE-only replication of the acceleration table,
run `exp_AccelerationNMSERecheck`. It uses 100 paired trials per Kp, a new fixed
seed bank, and eight `parfor` workers with four computational threads each.
Kp=16 and 64 run first. Per-trial checkpoints and paired 95% percentile
bootstrap intervals (2,000 resamples) are saved in
`results/acceleration_nmse_recheck_20260929/`; rerunning resumes saved trials.
The channel model, initialization, iteration limits and FFT-PCG final solver
match `exp_Acceleration`. Parallel wall times are **not** speedup measurements;
the original serial timing experiment and existing paper data are unchanged.

- Every experiment explicitly fixes the Twister seed. Methods compared at a
  point share channel, noise, and prior realizations. Coupling points share
  the latent random draws. Trial counts are independent Monte Carlo trials,
  not historical pilot snapshots.
- NMSE is the sum of squared estimation errors divided by the sum of channel
  energies, converted to dB only after aggregation. MAT files contain these
  per-trial energies, configurations, seeds, and MATLAB environment details.
  CSV values are unrounded. No smoothing or trial selection is applied.
- Initialization scales are 0.5, 1, 1.5 times both reference spreads; each
  scale computes its **own** moment-matched power initialization. Iteration 0
  is the initialized estimate. NLL is averaged and divided by N*Kp. Both
  NLL and NMSE are measured after each complete AO iteration.
- Ordinary runs retain AO/EM/L-BFGS maxima 3/3/2. The convergence experiment
  disables outer early stopping and records iterations 0 through 30.
  As in the original code, the EM limit counts SQUAREM iterations, not the
  number of underlying EM fixed-point evaluations.
- The acceleration table compares dense versus diagonal-likelihood parameter
  estimation, both followed by FFT-PCG. The timer includes initialization and
  required precomputation, excludes channel generation, final solves and
  plotting. Both paths are warmed up at each Kp; their execution order is
  alternated. Run this experiment alone. It fixes eight computational threads
  and restores the previous setting afterward. Runtime is hardware/load
  dependent, not bitwise reproducible.
- `acceleration.mat` separately records FFT-PCG and Cholesky final-solve times
  and errors at identical estimated parameters. It checks relative estimate
  error < 1e-6 and pooled NMSE difference < 0.001 dB. FFT-PCG tolerance is
  1e-8, maximum 5000 iterations; nonconvergence raises an error rather than
  silently substituting another solver.

## Model and baseline conventions

- Original CDL macro parameters and subpath generator are retained. The
  legacy receive-angle vectors contain the standard's **AOD column**, retained
  by agreement for the reciprocal-uplink angle convention; they are not the
  table's AOA column. This is a CDL-parameter-based simulation with the paper's
  intra-cluster profiles, not a full implementation of TR 38.901.
- Fixed-profile uses normalized exp(-cluster delay / DS) power and the same
  reference spread in every diffuse cluster. True DS and cluster delays are
  privileged inputs. The original empirical baseline used true per-cluster
  spreads; that is intentionally corrected here to match the revised paper.
- The copula experiment alone samples truncated exponential delay and
  truncated Laplacian angular marginals by inverse CDF, with equal-power rays.
  The latent Gaussian correlation couples delay to **absolute angle**; the
  angle sign is an independent fair draw. It is not the Pearson correlation
  between signed angle and delay. The ordinary original generator samples
  uniform offsets and weights ray powers instead; finite-ray realizations
  are not identical even at r=0.
- Copula Oracle uses the population joint PDAP and the same noisy CCM centers
  and uncertainty as the other estimators, not the actual realized rays or
  perfect centers. Conditional quadrature uses 128 nodes and 12 angular-prior
  nodes; its first realization at each nonzero r is checked with 256 nodes
  (relative estimate difference < 1e-3). At r=0 it is checked against the
  existing truncated separable Oracle (difference < 1e-6).
- Pilot-only methods use the same normalized pilot LS samples and periodic
  boundary convention, without CCM, DS, CP-based truncation, or history.
  The six candidates are linear, quadratic, fixed Hamming-windowed sinc FIR,
  PCHIP, spline, and causal IDFT/zero-padding/DFT interpolation. The FIR has
  8D+1 taps for interpolation ratio D, cutoff pi/D, zero-phase circular
  convolution and unity polyphase DC gains. DFT appends zeros to the causal
  IDFT sequence; it is **not** centered Fourier interpolation (`interpft`).
- The released Fig.1 explicitly uses DFT at every point, consistent with the
  response letter's observation that it wins at 7/10 points. The saved equal-weight
  mean **linear** NMSE is a separate diagnostic (it favors linear interpolation),
  not the final selection rule. No pointwise method switching is used. These
  comparisons do not establish universal superiority. The classical reference is Coleri et al.,
  IEEE Trans. Broadcasting, 2002, DOI: 10.1109/TBC.2002.804034; the fixed FIR
  parameters above specify this implementation, not an exact reproduction of
  that paper's filter design.

## Minimal validation

`test_OptimizerSafety` additionally checks non-descent directions, exhausted
bracketing/zoom, nonfinite trials, objective/parameter history alignment, and
16 dense/diagonal AO traces with one/four snapshots and four initializations.
It is also called by `run_all` and does not overwrite experiment results.

The 2026-09-25 optimizer fix keeps strong-Wolfe search as the primary path.
`SafeLineSearch` resets invalid/non-descent directions, accepts only finite
sufficient-decrease steps, and uses a verified Armijo fallback on search
exhaustion. If no such step exists, the spread update stops at its old point
with `Base:LineSearchFailed`, rather than accepting a failed trial. This is
not a guarantee of global optimality or monotone channel NMSE. EM diagnostic
histories now pair each evaluated objective with the powers at that point
(including rejected extrapolation trials); the EM update itself is unchanged.
The channel generator and continuous-PDAP Oracle are unchanged. Existing
historical outputs predate this fix; subsequent replay checks are recorded in
`results/release_audit/REPORT.md`. Use a fresh output directory for post-fix
experiments, not a pre-fix acceleration checkpoint.

`test_Numerics` checks seeded draw reproducibility, the exact FFT Toeplitz
operator against a dense matrix, FFT-PCG against Cholesky, interpolation at
pilot positions/constant channels/full pilot grids, online/history calls,
and copula marginal support and uniform quantiles. Original public Monte
Carlo outputs 1--6 retain their ordering; outputs 7--10 add parameter timing,
final-solve timing, pilot-only NMSE, and per-trial energy statistics.

No manuscript files are changed by any entry. This package does not depend
on, copy, or modify the exploratory result caches or original release folder.

## Response-letter supplementary experiments (2026-09-25)

Run `run_ResponseExperiments` for the four supplementary studies (100 trials
per configuration), or `run_ResponseExperiments(2,'results/response_smoke')`
for a short check. The original `run_all` is unchanged. New results go to
`results/response_supplement_20260925/`, never over the earlier paper results.
Each experiment below accepts `(mc_num,folder)`; each plotting entry accepts
only `folder` and never runs an experiment.

| Experiment entry | Configuration | Data-only plotting entry |
|---|---|---|
| `exp_ResponseLoS` | CDL-D branch; 5 dB; Kp=32,64 | Reuses the existing Kp scan plot in `los_scan/` |
| `exp_CCMPriorMismatch` | CDL-A; 5 dB; Kp=32; actual/assumed error SD=0.5,1,2,4, delay and angle separately | `plot_CCMPriorMismatch` |
| `exp_ParameterAblation` | CDL-A; 5 dB; Kp=32,64; fixed/power-only/spread-only/joint/Oracle | `plot_ParameterAblation` |
| `exp_ParameterAccuracy` | Default: manuscript finite truncated subpaths; RMS-width accuracy; 5 dB; Kp=32,64; 1/8 training snapshots and independent test channel. Third argument `false` selects the older Gaussian control. | `plot_ParameterAccuracy` |

The updated 20-trial accuracy study is stored separately in
`results/parameter_accuracy_truncated/`; its README specifies the RMS-width
targets and reproduction command. Historical Gaussian-control results below
remain unchanged and should not be confused with this updated experiment.

Ordinary AO/SQUAREM/L-BFGS limits remain 3/3/2, with the existing stopping
rules and safeguarded line search. Both accelerations remain enabled. No
hyperparameter is selected using the supplementary outcomes.

Implementation differences and interpretation:

- `estor.sigma_tau_hat` and `estor.sigma_theta_hat` now actually control
  CPM/Fixed-profile/Proposed's assumed errors. Previously these stored fields
  were unused. Channel/prior draws and Oracle still use `chann.sigma_*`.
  Equal settings reproduce the pre-change default numerically. The mismatch
  scan fixes assumed SDs at 3 ns/3 degrees, shares channel/noise/standardized
  prior-error draws at all scales, and retains the first-delay anchor.
  Pilot-only must be identical across all these scales; Oracle knows actual
  error statistics and true profiles, not the true realized cluster centers.
- Optional `TraceAO` fields `fit_powers` and `fit_spreads` default to true.
  Setting either false skips that update without changing the other update.
  Histories mark skipped stages `PowerFixed`/`SpreadFixed`. Power-only and
  joint fitting use the same moment-matched initial powers and 1x spreads;
  spread-only keeps normalized exp(-delay/DS) powers fixed. No fitted powers
  are renormalized after optimization. These are controlled ablations, not
  new implementations of the full estimator.
- LoS reuses the existing CDL-D macro-parameter branch (14 components,
  including the specular component), DS=100 ns, reference spreads 5 ns and
  3 degrees, and a random uniform specular phase. The specular covariance
  has zero intra-component spreads. This is not a deterministic known-mean
  Rician model or a complete TR 38.901 implementation. Existing angle-column
  and first-delay-anchor conventions are retained.
- The older Gaussian-control accuracy study (third argument `false`) changes the channel law: it draws zero-mean
  complex-Gaussian channels from the full, untruncated analytic covariance,
  including the 3 ns/3-degree uncertainty integrated around fixed nominal
  means. It does not perturb those means again or use the diagonal surrogate
  to generate channels. True powers/spreads are the existing CDL-A-based
  vectors; no blockages. Training snapshots are independent conditional on
  one fixed covariance; the one-snapshot case is the first column of the
  eight-snapshot case. The test channel and its noise are independent of
  training and shared across the compared fits. Consequently these held-out
  NMSE values are not the paper's same-snapshot online NMSE.
- `GaussianCovarianceTools` samples through separable covariance roots and
  evaluates the full noise-free covariance error using exact Hermitian
  Toeplitz lag multiplicities, without forming an NK-by-NK matrix. Only
  negative eigenvalues at floating-point roundoff are clipped for sampling;
  this is not a PSD projection added to the estimator. Tests compare the
  norm calculation with dense matrices and check empirical covariance.
- Accuracy metrics are norm(p_hat-p)/norm(p), RMS natural log of the
  delay-spread ratio, RMS natural log of the angular-spread ratio, and
  norm(C_hat-C,'fro')/norm(C,'fro'). Initialization and fitted metrics are
  both saved; weak clusters are not silently omitted from the spread errors.
  Parameter accuracy, covariance accuracy and channel NMSE are different
  criteria; none proves identifiability or global optimality.

All methods at a configuration use paired trials. MAT files contain raw
errors/energies, configurations, seeds, relevant source snapshots, and (for
ablation/accuracy) estimated parameters and AO histories. `ResponseStatistics`
uses 2000 paired bootstrap resamples with a separate random stream to supply
95% percentile intervals; NMSE is the ratio of pooled energies, not mean dB.
MAT fields `bootstrap` permit paired method-difference intervals; the ablation
also saves `gain_db`/`gain_ci`. Accuracy metrics use across-trial means.
CSV files report unrounded central values. Plot error bars are 95% bootstrap
intervals, not standard deviations. Figures are PNG plus visible editable FIG,
never additional PDF exports. `test_ResponseExperiments` checks unchanged
default AO, frozen blocks, actual/assumed prior separation, Gaussian sampling,
structured covariance norms, and isolation of the bootstrap random stream.

### CCM-variance compensation control

Run `exp_CCMVarianceCompensation` for 100 paired trials at 5 dB and Kp=32;
redraw with `plot_CCMVarianceCompensation`. Outputs are saved separately in
`results/ccm_variance_compensation`, without replacing the earlier mismatch
results. The a-grid remains [0.5,1,2,4], where a is the actual/assumed prior
error **standard-deviation ratio**, not the variance ratio (which is a^2).
The optional third argument sets process workers (default 4; 0 runs serially).
All channel/noise/prior draws occur before parallel fitting, so scheduling
does not change the samples. Each fitting worker uses two computational threads.

For each trial, fit Proposed once at a=1. Freeze that trial's powers and both
spread vectors in the control, but use the current a-dependent noisy CCM
centers at every point. Proposed instead refits all parameters from the same
current pilot. Both retain assumed prior SDs of 3 ns and 3 degrees; Oracle
uses actual prior SDs and the true truncated profile, not exact centers.
Channels, noise and standardized prior errors are paired across the scan.
The original finite-subpath generator, accelerated AO and FFT-PCG are unchanged.

Left axes compare Proposed, Oracle and the frozen-a=1 control. Right axes
show the fitted delay spread (ns) or angular spread (degrees), averaged
equally across clusters and trials, with 95% bootstrap intervals. Full fitted
vectors and power-weighted means are also saved. A shift in mean spread alone
does not prove beneficial compensation: examine paired NMSE gains and their
intervals as well. The frozen control is a counterfactual diagnostic using a
matched prior for calibration, not a deployable method under unknown mismatch.

### Fixed-observation, spread-only compensation pilot study

`exp_AssumedVarianceCompensation(20,folder,8)` runs the quick 20-trial study
at 15 dB and Kp=64; `plot_AssumedVarianceCompensation(folder)` redraws it.
Default folder: `results/assumed_variance_compensation`. Unlike the previous
scan, the true channel, pilot, CCM centers and actual error SDs (3 ns / 3 deg)
stay fixed. Scan assumed SD = actual SD / a for a=[0.25,0.5,1,2,4], one
dimension at a time. Oracle is exactly invariant across the scan.

At a=1, a joint fit calibrates one shared power vector per trial. Hold this
power fixed, then fit both spreads at each a with identical 1x initialization
and the same 3/3/2 iteration limits. The a=1 spread-only fit supplies the
frozen-profile control. Thus the curves compare spread adaptation, not the
unmodified joint Proposed estimator. No channel-generation change is made.

This experiment opts into `ao_opts.spread_bounds`, a 2-by-2 matrix with
delay/angle rows and lower/upper columns. Bounds are fixed in advance at
[0.05,4] times the reference factors: delay [0.55,44] ns, angle [0.5,40] deg.
Existing L-BFGS line search rejects out-of-box candidate points; final values
are not clipped and no trials are excluded. This is an experimental algorithm
constraint, not a 3GPP prescription. Without the option, Base retains its
original behavior. Saved histories and full parameter arrays allow inspection
of objective descent, frozen powers and near-boundary estimates.
