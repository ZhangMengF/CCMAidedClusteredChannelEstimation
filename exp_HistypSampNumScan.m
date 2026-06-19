%EXP_HISTYPSAMPNUMSCAN1 Scan historical-snapshot count under blockage.
% The script evaluates the proposed estimator for several historical sample
% counts and blockage probabilities. The scan function saves the figures.

clear; clc; close all;

[sys, chann, simu, estor, funcs] = Base(false);
sys.Kp = 64;

S_values = [0, 1, 3, 7, 15];      % Historical pilot snapshot counts.
simu.MCNum    = 100;               % Monte Carlo trials per scan point.
MCNMSEsSwitch = [false; false; false; false; true; false];  % Proposed estimator only.

% Alternative scan configurations:
% 'SNR_dB',    [-5, 5, 15];
% 'BlockProb', [0, 0.3, 0.5];
scan_name   = 'BlockProb';
scan_values = [0.5, 0.3, 0];
HisthpNumVsNMSE_ParaScan(sys, chann, simu, estor, funcs, S_values, MCNMSEsSwitch, scan_name, scan_values);
