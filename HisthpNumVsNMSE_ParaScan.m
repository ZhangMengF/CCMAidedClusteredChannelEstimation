function HisthpNumVsNMSE_ParaScan(sys, chann, simu, estor, funcs, S_values, MCNMSEsSwitch, scan_name, scan_values)
%HISTHPNUMVSNMSE_PARASCAN Scan historical-snapshot count and one parameter.
% Inputs:
%   sys           - OFDM and antenna-array parameters.
%   chann         - Channel, prior, spread, SNR, and blockage parameters.
%   simu          - Monte Carlo simulation settings.
%   estor         - Estimator optimization settings.
%   funcs         - Function handles exported by Base.
%   S_values      - Historical pilot snapshot counts.
%   MCNMSEsSwitch - Six-element logical vector selecting estimators.
%   scan_name     - Parameter name: 'SNR_dB', 'Kp', or 'BlockProb'.
%   scan_values   - Values defining the plotted curves.
% Output:
%   None. The function displays and saves PNG and FIG result files.

num_S  = length(S_values);
num_J  = length(scan_values);
NMSE_esti_all = zeros(num_S, num_J);

start = tic;
for j = 1:num_J
    % Apply the parameter value associated with the current curve.
    [sys, chann] = apply_scan_param(sys, chann, scan_name, scan_values(j));

    % Construct a uniform pilot pattern for the current Kp.
    pilot_sc  = 1:sys.K/sys.Kp:sys.K;
    pilot_ant = (1:sys.Np).';

    for i = 1:num_S
        S = S_values(i);

        rng(1000*j + 100*i);  % Use a deterministic seed for reproducibility.
        [~, ~, ~, ~, NMSE_esti_all(i, j), ~] = ...
            funcs.CalcuNMSEsByMonteCarlo(sys, chann, simu, estor, pilot_sc, pilot_ant, S, MCNMSEsSwitch);

        fprintf('\n====== [%s] %s=%g, S=%d (%d/%d, %d/%d) finished, 历时%.1fs ======\n', ...
            scan_name, scan_name, scan_values(j), S, j, num_J, i, num_S, toc(start));
    end
end

% Plot one curve for each scanned parameter value.
markers = {'o-', 's-', '^-', 'v-', 'd-', 'p-', 'h-'};
f = figure('Position', [100 100 800 500]);
hold on;
for j = 1:num_J
    mk = markers{mod(j-1, length(markers)) + 1};
    plot(S_values, 10*log10(NMSE_esti_all(:, j)), mk, ...
        'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', format_legend(scan_name, scan_values(j)));
end
xlabel('Historical Pilot Snapshot Count (S)');
ylabel('NMSE (dB)');
grid on;
legend('Location', 'northeast');

% Build a title containing the simulation parameters.
if chann.is_LoS, chann_name = 'CDL-D(LoS)'; else, chann_name = 'CDL-A(NLoS)'; end
pilot_str = 'PilotPatt=Uniform';
param_str1 = sprintf(['N=%d,K=%d,Np=%d,Kp=%d,{\\Delta}f=%dkHz,', ...
    '%s,SNR=%ddB,\\sigma_{\\tau}=%.1fns,c_{ds}=%.1fns(mean),\\sigma_{\\theta}=%.1fdeg,c_{asa}=%.1fdeg(mean),SubpathNum=%d,SpreadTruncFactor=%d,%s,MCNum=%d,BlockProb=%.1f'], ...
    sys.N,sys.K,sys.Np,sys.Kp,sys.Deltaf/1e3,chann_name,chann.SNR_dB, ...
    chann.sigma_tau*1e9,chann.c_ds_ref*1e9,chann.sigma_theta*180/pi,chann.c_asa_ref*180/pi,chann.M_subpath,chann.SpreadTruncFactor,pilot_str,simu.MCNum,chann.ClusterBlockProb);
param_str2 = sprintf('AO\\_MaxIter=%d,EM\\_MaxIter=%d,LBFGS\\_MaxIter=%d',estor.ao_opts.ao_max_iter,estor.ao_opts.em_max_iter,estor.ao_opts.lbfgs_max_iter);

title({['\fontsize{8}\rm ', param_str1], ...
    ['\fontsize{8}\rm ', param_str2], ...
    ['\fontsize{8}\bf AO Estimated NMSE vs S  (scan: ' scan_name ')']});

% Save the figure in the current working directory.
currentTime = char(datetime('now', 'Format', 'MMdd_HHmm'));
filename_png = sprintf('HistSampScan_%s_%s.png', scan_name, currentTime);
filename_fig = sprintf('HistSampScan_%s_%s.fig', scan_name, currentTime);
saveas(f, fullfile(pwd, filename_png));
saveas(f, fullfile(pwd, filename_fig));
fprintf('\nFigure saved: %s, %s\n', filename_png, filename_fig);

end


%% Local helper functions

function [sys, chann] = apply_scan_param(sys, chann, scan_name, val)
%APPLY_SCAN_PARAM Update one simulation parameter.
% Inputs:
%   sys, chann - Current system and channel structures.
%   scan_name  - Parameter name: 'SNR_dB', 'Kp', or 'BlockProb'.
%   val        - New parameter value.
% Outputs:
%   sys, chann - Updated system and channel structures.

switch scan_name
    case 'SNR_dB'
        chann.SNR_dB = val;
    case 'Kp'
        sys.Kp = val;
    case 'BlockProb'
        chann.ClusterBlockProb = val;
    otherwise
        error('Unsupported scan_name: %s. Use ''SNR_dB'', ''Kp'', or ''BlockProb''.', scan_name);
end
end

function s = format_legend(scan_name, val)
%FORMAT_LEGEND Create a legend label for one scanned parameter value.
% Inputs:
%   scan_name - Parameter name.
%   val       - Parameter value.
% Output:
%   s         - Formatted legend string.

switch scan_name
    case 'SNR_dB'
        s = sprintf('SNR = %g dB', val);
    case 'Kp'
        s = sprintf('K_p = %d', val);
    case 'BlockProb'
        s = sprintf('Proposed, p_{block} = %.1f', val);
    otherwise
        s = sprintf('%s = %g', scan_name, val);
end
end
