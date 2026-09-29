%EXP_KPVSMONTECARLONMSE Scan the pilot-subcarrier count and plot NMSE.
% Compare CPM-aided, Fixed-profile, Proposed, Oracle and Pilot-Only for CDL-A.
% Figures and raw statistics are saved by MonteCarloNMSE_ScanParas.

clear; clc; close all;

[sys, chann, simu, estor, funcs] = Base(false);

%% Experiment configuration
scan_name = 'Kp';
scan_values = [8, 16, 32, 64, 128];
estor.pilot_only_method = 'dft'; % dft | linear | quadratic | fir | pchip | spline

simu.MCNum = 50;                      % Monte Carlo trials per scan point.
S_hist = 0;                           % Zero selects online estimation.
MCNMSEsSwitch = [false,false,true,true,true,true,true]; % CPM, fixed, proposed, oracle, pilot-only.

MonteCarloNMSE_ScanParas(sys, chann, simu, estor, funcs, S_hist, MCNMSEsSwitch, scan_name, scan_values);
