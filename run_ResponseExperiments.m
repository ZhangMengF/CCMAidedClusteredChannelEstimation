function run_ResponseExperiments(mc_num,folder)
%RUN_RESPONSEEXPERIMENTS Four standalone response studies; no manuscript writes.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','response_supplement_20260925'); end
old_threads=maxNumCompThreads(8);
cleanup=onCleanup(@() maxNumCompThreads(old_threads));
test_Numerics; test_OptimizerSafety; test_ResponseExperiments;
exp_ResponseLoS(mc_num,folder);
exp_CCMPriorMismatch(mc_num,folder); plot_CCMPriorMismatch(folder);
exp_ParameterAblation(mc_num,folder); plot_ParameterAblation(folder);
exp_ParameterAccuracy(mc_num,folder); plot_ParameterAccuracy(folder);
end
