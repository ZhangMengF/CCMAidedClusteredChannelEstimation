%EXP_SNRVSMONTECARLONMSE Scan SNR and plot Monte Carlo NMSE.
% Compare CPM-aided, Fixed-profile, Proposed, Oracle and Pilot-Only for CDL-A.
% Figures and raw statistics are saved by MonteCarloNMSE_ScanParas.

clear; clc; close all;

[sys, chann, simu, estor, funcs] = Base(false);

%% Experiment configuration
scan_name = 'SNR_dB';
scan_values = [-5, 0, 5, 10, 15];
estor.pilot_only_method = 'dft'; % dft | linear | quadratic | fir | pchip | spline

simu.MCNum = 100;                     % Monte Carlo trials per scan point.
S_hist = 0;                           % Zero selects online estimation.
MCNMSEsSwitch = [false,false,true,true,true,true,true]; % CPM, fixed, proposed, oracle, pilot-only.

MonteCarloNMSE_ScanParas(sys, chann, simu, estor, funcs, S_hist, MCNMSEsSwitch, scan_name, scan_values);
