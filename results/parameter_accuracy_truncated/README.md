# Truncated-profile parameter accuracy

Reproduce from `OpenSourceCodes_new` in MATLAB:

```matlab
maxNumCompThreads(4);
exp_ParameterAccuracy(20);
plot_ParameterAccuracy(fullfile(pwd,'results','parameter_accuracy_truncated'));
```

20 paired trials, Base defaults (5 dB), Kp=32/64, one/eight independent fitting snapshots and an independent test snapshot. Seeds are 280925+(1:20). The first fitting snapshot is shared; initialized/fitted/Oracle methods share the test data. No blockage.

The generator is Base.GenerateChannelObservation: finite uniformly sampled subpath locations, truncated exponential/Laplacian power weights and independent random phases. Each snapshot draws Gaussian center perturbations around fixed CCM nominal centers (first delay held fixed), following Base's historical-mode convention. This replaces the earlier Gaussian covariance sampler, not the estimation algorithm. AO/EM/L-BFGS limits and accelerated fitting/FFT-PCG remain unchanged; no new spread bounds.

Spread accuracy compares centered RMS widths of the generating truncated population profiles with those of the fitted infinite-tail profiles. It does not compare the original pre-truncation scale parameters or the random empirical width of one finite realization. The targets are Base.c_ds_eq_vec and Base.c_asa_eq_vec; fitted RMS widths equal their respective fitted scale parameters. The metric is the root mean squared natural log fitted/target ratio over all 23 clusters (delay and angle separately). Center uncertainty is excluded from intra-cluster widths. Power and covariance errors remain relative Euclidean/Frobenius norms. Oracle and covariance-error reference use the analytic truncated profile, as in the manuscript; finite normalized subpath sampling approximates this population spectrum.

Raw trial metrics, estimates, seeds, source, histories, bootstrap intervals and NMSE are saved in MAT/CSV; PNG and visible FIG are generated separately. Old results are preserved. For the previous Gaussian control use an explicit separate output folder and `exp_ParameterAccuracy(n,folder,false)`.
