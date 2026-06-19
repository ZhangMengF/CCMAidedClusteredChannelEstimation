function MonteCarloNMSE_ScanParas(sys, chann, simu, estor, funcs, S_hist, MCNMSEsSwitch, scan_name, scan_values)
%MONTECARLONMSE_SCANPARAS Scan one parameter and plot estimator NMSE.
% Inputs:
%   sys           - OFDM and antenna-array parameters.
%   chann         - Channel, prior, spread, SNR, and blockage parameters.
%   simu          - Monte Carlo simulation settings.
%   estor         - Estimator optimization settings.
%   funcs         - Function handles exported by Base.
%   S_hist        - Historical snapshot count; zero selects online mode.
%   MCNMSEsSwitch - Six-element logical vector selecting plotted estimators.
%   scan_name     - Parameter name: 'SNR_dB' or 'Kp'.
%   scan_values   - Values of the selected parameter.
% Output:
%   None. The function displays and saves PNG and FIG result files.

scan_len = length(scan_values);
NMSE_Oracle_equa_all      = zeros(scan_len, 1);
NMSE_sigmaEqu0_cEqu0_all  = zeros(scan_len, 1);
NMSE_cEqu0_all            = zeros(scan_len, 1);
NMSE_emp_all              = zeros(scan_len, 1);
NMSE_esti_all             = zeros(scan_len, 1);
NMSE_Oracle_all           = zeros(scan_len, 1);

for i = 1:scan_len

    worker_tic = tic;
    chann_woker = chann;
    sys_worker = sys;
    if strcmp(scan_name, 'SNR_dB')
        chann_woker.SNR_dB = scan_values(i);
    elseif strcmp(scan_name, 'Kp')
        sys_worker.Kp = scan_values(i);
    else
        error('Unsupported scan_name: %s', scan_name);
    end

    pilot_sc = 1:sys_worker.K/sys_worker.Kp:sys_worker.K;
    pilot_ant = (1:sys_worker.Np).';

    rng(100 + i);  % Use a deterministic seed for reproducibility.
    [NMSE_Oracle_equa_all(i), NMSE_sigmaEqu0_cEqu0_all(i), NMSE_cEqu0_all(i), ...
        NMSE_emp_all(i), NMSE_esti_all(i),NMSE_Oracle_all(i)] = ...
        funcs.CalcuNMSEsByMonteCarlo(sys_worker, chann_woker, simu, estor, pilot_sc, pilot_ant, S_hist, MCNMSEsSwitch);

    fprintf('\n====== %s = %g (%d/%d) finished, consumed %.1fs ======\n\n', scan_name, scan_values(i), i, scan_len, toc(worker_tic));
end

%% Plot results
f = figure('Position', [100 100 800 500]);
hold on;
if MCNMSEsSwitch(1)
    plot(scan_values, 10*log10(NMSE_Oracle_equa_all), 'o-', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Oracle (True Power + Eq. Spread)');
end
if MCNMSEsSwitch(2)
    plot(scan_values, 10*log10(NMSE_sigmaEqu0_cEqu0_all), 'd-', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', '\sigma=0, c=0 (No Prior Uncert.)');
end
if MCNMSEsSwitch(3)
    plot(scan_values, 10*log10(NMSE_cEqu0_all), '^-', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'CPM-aided');
end
if MCNMSEsSwitch(4)
    plot(scan_values, 10*log10(NMSE_emp_all), 'v-', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Empirical Power and Spread');
end
if MCNMSEsSwitch(5)
    plot(scan_values, 10*log10(NMSE_esti_all), 's-', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Proposed');
end
if MCNMSEsSwitch(6)
    plot(scan_values, 10*log10(NMSE_Oracle_all), 'o-', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Oracle');
end

if strcmp(scan_name, 'SNR_dB')
    xlabel('SNR (dB)');
elseif strcmp(scan_name, 'Kp')
    xlabel('Pilot Subcarrier Number (Kp)');
else
    xlabel(scan_name);
end
ylabel('NMSE (dB)');
grid on;
legend('Location', 'northeast');

% Build a title containing the fixed simulation parameters.
if chann.is_LoS
    chann_name = 'CDL-D(LoS)';
else
    chann_name = 'CDL-A(NLoS)';
end
pilot_str = 'PilotPatt=Uniform';

param_str1 = sprintf(['N=%d,K=%d,Np=%d,Kp=%d,{\\Delta}f=%dkHz,', ...
    '%s,SNR=%ddB,\\sigma_{\\tau}=%.1fns,c_{ds}(mean)=%.1fns,\\sigma_{\\theta}=%.1fdeg,c_{asa}(mean)=%.1fdeg,SubpathNum=%d,SpreadTruncFactor=%d,%s,MCNum=%d,BlockProb=%.1f'], ...
    sys.N,sys.K,sys.Np,sys.Kp,sys.Deltaf/1e3,chann_name,chann.SNR_dB, ...
    chann.sigma_tau*1e9,chann.c_ds_ref*1e9,chann.sigma_theta*180/pi,chann.c_asa_ref*180/pi,chann.M_subpath,chann.SpreadTruncFactor,pilot_str,simu.MCNum,chann.ClusterBlockProb);

param_str2 = sprintf('AO\\_MaxIter=%d,EM\\_MaxIter=%d,LBFGS\\_MaxIter=%d',estor.ao_opts.ao_max_iter,estor.ao_opts.em_max_iter,estor.ao_opts.lbfgs_max_iter);

title({['\fontsize{8}\rm ', param_str1], ...
    ['\fontsize{8}\rm', param_str2], ...
    ['\fontsize{8}\bf MonteCarlo NMSE vs ', scan_name]});

% Save the figure beside this function file.
currentTime = char(datetime('now', 'Format', 'MMdd_HHmm'));
filename_png = sprintf('MonteCarloNMSE_Scan_%s_%s.png', scan_name, currentTime);
filename_fig = sprintf('MonteCarloNMSE_Scan_%s_%s.fig', scan_name, currentTime);
saveas(f, fullfile(fileparts(mfilename('fullpath')), filename_png));
saveas(f, fullfile(fileparts(mfilename('fullpath')), filename_fig));
fprintf('\nFigure saved: %s, %s\n', filename_png, filename_fig);

end
