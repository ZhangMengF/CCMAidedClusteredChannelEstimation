%EXP_KPVSMONTECARLONMSE Scan the pilot-subcarrier count and plot NMSE.
% The script compares the CPM-aided, proposed, and oracle estimators for
% CDL-A channels. Figures are saved by MonteCarloNMSE_ScanParas.

clear; clc; close all;

[sys, chann, simu, estor, funcs] = Base(false);

%% Experiment configuration
scan_name = 'Kp';
scan_values = [8, 16, 32, 64, 128];

simu.MCNum = 50;                      % Monte Carlo trials per scan point.
S_hist = 0;                           % Zero selects online estimation.
MCNMSEsSwitch = [false,false,true,false,true,true]; % Enable CPM, proposed, and oracle.

MonteCarloNMSE_ScanParas(sys, chann, simu, estor, funcs, S_hist, MCNMSEsSwitch, scan_name, scan_values);
