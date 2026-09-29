%RUN_ALL Reproduce the four added experiments (no manuscript writes).
% Run this directory alone during the timing experiment.
clear; clc;
maxNumCompThreads(8);
test_Numerics;
test_OptimizerSafety;
exp_PilotOnlyComparison;
exp_InitializationConvergence;
exp_CopulaCoupling;
exp_Acceleration;
PlotResults;
