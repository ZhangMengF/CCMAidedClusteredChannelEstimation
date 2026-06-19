# CCM-Aided Clustered Channel Estimation

MATLAB implementation of the simulations for robust channel estimation using Channel Cluster Map (CCM) information.

## Requirements

- MATLAB R2025b
- No additional toolboxes are required.

## Files

- `Base.m`: System/channel configuration, channel generation, covariance construction, and estimation algorithms.
- `exp_KpVsMonteCarloNMSE.m`: NMSE versus pilot-subcarrier number.
- `exp_SNRVsMonteCarloNMSE.m`: NMSE versus SNR.
- `MonteCarloNMSE_ScanParas.m`: Common parameter-scan routine.
- `exp_HistypSampNumScan.m`: NMSE versus historical pilot-snapshot number under blockage.
- `HisthpNumVsNMSE_ParaScan.m`: Historical-snapshot scan routine.

## Usage

Place all six files in the same directory, set this directory as the MATLAB current folder, and run one of the following experiment scripts:

```matlab
exp_KpVsMonteCarloNMSE
exp_SNRVsMonteCarloNMSE
exp_HistypSampNumScan
The scripts use the CDL-A channel model and save the generated figures as .png and .fig files.
Simulation parameters, Monte Carlo trial counts, scan values, and enabled estimators can be configured directly in the corresponding experiment script.
In the historical-snapshot experiment, S = 0 denotes online estimation using the current pilot, whereas S > 0 denotes offline estimation using historical pilot snapshots.
```
