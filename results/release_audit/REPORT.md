# Submission data audit — 2026-09-29

## Scope

Neither manuscript nor response-letter TEX was changed. Original submission
and original open-source code remain read-only. Base.m is unchanged; new work
adds a paired, checkpointed Fig.1 execution wrapper and reproducibility checks,
not a different channel generator or estimator.

## Fig.1

All 25 configurations completed in `../combined_accelerated/`:
100 trials per SNR point, 50 per Kp point, and 100 per historical-snapshot /
blockage configuration. Both online and offline Proposed fitting use the
beam-delay diagonal acceleration; final LMMSE solves use FFT-PCG. Pilot-only
is DFT. Raw errors, energies, seeds, configuration and source snapshot are saved.

The frozen plotting template retains fonts, line styles, legend ordering and
axes placement. The sole approved layout exception extends panel c's upper
Y limit from -5.5 to -5 dB: its new first offline point is -5.023587 dB.
The original -0.5 dB tick spacing is retained.
The checked PNG has replaced `0817/Result.png`; the previous PNG is retained
as `figure_assets/Fig1_before_accelerated.png`. No TEX or manuscript PDF was
rewritten. The updated image will appear on the next manuscript compilation.

At Kp=16 Proposed is 0.00274 dB worse than Fixed-profile in this paired run;
this tiny reversal is retained, not retuned or hidden. R1.2's rounded illustrative
values at 5 dB, Kp=32 remain correct (-5.79 / -6.49 dB).

## Reproduction checks

`check_SubmissionReplay` saves `replay_checks.mat`. Raw squared-error differences
were exactly zero in every replay listed below (same seeds):

| Result | Replay scope | Outcome |
|---|---|---|
| Initialization | Complete 20 trials, all 31 recorded iterates and three initializations | Exact raw-error and NLL match |
| Pilot-only comparison | Complete 100/50 trials per SNR/Kp point, six methods | Exact raw-error match |
| LoS | First trial at each Kp | Exact raw-error match |
| Power/spread ablation | First trial at each Kp, all methods | Exact raw-error match |
| Statistical accuracy | First trial per Kp, all snapshot conditions | Exact raw-error match |
| Copula coupling | First trial per correlation value | Exact raw-error match |

Representative-trial checks are not full Monte Carlo reruns. The LoS, ablation,
accuracy and coupling table values were independently checked against their
complete stored datasets and agree after rounding. The 100-trial assumed-error
mismatch figure points to its correct dataset; its Base snapshot equals current
Base, and the reported trends/numbers agree with that dataset. That expensive
100-trial experiment was not rerun during this audit.

`test_Numerics` and `test_OptimizerSafety` passed. FFT matvec relative error was
5.75e-16 and FFT-PCG solve relative error was 3.36e-9 in the numerical test.

## Remaining discrepancies — not silently corrected

1. **Main Table I is not internally aligned with one dataset.** Its absolute
   accelerated NMSE and speedups match the old acceleration dataset, but the
   Kp=16 difference 0.15 dB comes from the newer 100-trial recheck, whose absolute
   NMSE differs. Old Kp=16 difference is 0.2923 dB. At Kp=64 the stated 0.02 dB
   matches neither old -0.0111 dB nor new 0.03925 dB. New 100-trial Kp=8/32
   differences are 0.20052/0.02968 dB; old values are 0.08713/0.02338 dB.
   Reconcile the whole table after the separately running Kp=128 job completes.
2. **Runtime provenance:** a pre-existing Kp=128 timing experiment remained
   running while this audit ran. We did not stop or alter it. Timing samples
   overlapping concurrent workloads are not isolated benchmark measurements.
3. **Pilot-only graphic:** stored data and full replay agree exactly, but a fresh
   `PlotResults` export has different native-legend spacing/crop from the approved
   PNG (2699x998 versus 2728x1020). The response's approved PNG was not replaced.
   Numerical reproducibility is confirmed; pixel-identical layout is not.
4. The initialization figure *does* redraw pixel-identically to the main-paper
   PNG using current `PlotResults` and saved data, despite different PNG metadata.
5. DFT wins 7/10 pilot-only test conditions; the equal-weight linear-domain
   diagnostic score instead favors linear interpolation. These are different
   selection criteria; Fig.1 explicitly uses DFT, not pointwise selection.

## Protected TEX hashes

- `0817/main.tex`: `8d8401808dcf99658d0e1257074cf9fa69d6ac24209f735c71f744dca1f147fe`
- `0817/ResponseLetter/response_to_reviewers.tex`: `3545fe1146f58e2c00df4de83f5862145fdfe57d05da319593feac4f92589846`

See `../../RELEASE_REPRODUCIBILITY.md` for exact experiment and plotting entries.
