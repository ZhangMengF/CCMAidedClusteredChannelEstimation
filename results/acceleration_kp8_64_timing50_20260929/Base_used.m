function [sys, chann, simu, estor, funcs] = Base(is_LoS)
%BASE Initialize the simulation environment and export utility functions.
% Input:
%   is_LoS - True for the CDL-D LoS configuration; false for CDL-A NLoS.
% Outputs:
%   sys    - OFDM and antenna-array parameters.
%   chann  - Channel-model, prior-uncertainty, and blockage parameters.
%   simu   - Monte Carlo simulation parameters.
%   estor  - Estimator initialization and optimization options.
%   funcs  - Function handles to the local simulation and estimation tools.

% ========================= System Parameters =========================
sys.N = 32;                  % Number of receive antennas.
sys.K = 256;                 % Total number of subcarriers.
sys.Deltaf = 960e3;          % Subcarrier spacing.
sys.f = ((-sys.K/2):(sys.K/2-1)).' * sys.Deltaf;
sys.fc = 28e9;
sys.c = 3e8;
sys.lambda = sys.c / sys.fc;
sys.d = sys.lambda / 2;
sys.Kp = 32;
sys.Np = sys.N;

fprintf('时延分辨率:%.1f ns, 角度分辨率:%.1f deg \n',1e9/(sys.K*sys.Deltaf),2*180/(sys.N*pi))

% ========================= Channel Parameters =========================
if is_LoS
    [delay_norm, power_dB, aoa_deg] = get_3gpp_CDL_D();
    chann.SNR_dB = 20;
    chann.c_ds_ref  = 5e-9;
    chann.c_asa_ref = 3 * pi/180;
else
    [delay_norm, power_dB, aoa_deg] = get_3gpp_CDL_A();
    chann.SNR_dB = 5;
    chann.c_ds_ref  = 11e-9;
    chann.c_asa_ref = 10 * pi/180;
end

chann.is_LoS = is_LoS;
chann.DS = 100e-9;
chann.tau = delay_norm * chann.DS;
chann.theta = deg2rad(aoa_deg);
chann.sigma_tau = 3e-9;
chann.sigma_theta = 3 * pi / 180;
cluster_powers = 10.^(power_dB/10);
chann.cluster_powers_NoBlockage = cluster_powers / sum(cluster_powers);
chann.M_subpath = 40;
L_chan = numel(chann.tau);

%% ===== Per-cluster spread factors =====
gamma = 0.10;
w         = chann.cluster_powers_NoBlockage .^ (-gamma);
A_ds      = chann.c_ds_ref  * sum(w) / sum(chann.cluster_powers_NoBlockage .^ (1 - gamma));
A_asa     = chann.c_asa_ref * sum(w) / sum(chann.cluster_powers_NoBlockage .^ (1 - gamma));
chann.c_ds_vec  = A_ds  .* w / sum(w);
chann.c_asa_vec = A_asa .* w / sum(w);

% For LoS: cluster 1 is the specular component (no intra-cluster spread)
if is_LoS && L_chan >= 1
    chann.c_ds_vec(1)  = 0;
    chann.c_asa_vec(1) = 0;
end

DoF_delay = 2*chann.c_ds_ref*(sys.K*sys.Deltaf);
DoF_AoA   = 4*chann.c_asa_ref/(2/sys.N);
fprintf('平均每簇参数自由度: 时延-%.1f, 角度-%.1f, 总-%.1f \n', ...
    DoF_delay, DoF_AoA, DoF_delay*DoF_AoA)

% Construct a delay-based empirical power profile for a baseline estimator.
emp_powers = exp(-chann.tau / chann.DS);
chann.cluster_powers_emp = emp_powers / sum(emp_powers);
% Default blockage probability.
chann.ClusterBlockProb = 0;

chann.SpreadTruncFactor = 2;

rt = chann.SpreadTruncFactor;
s_ds  = sqrt(1 - rt^2 * exp(rt) / (exp(rt) - 1)^2);
s_asa = sqrt(1 - (rt^2 + sqrt(2)*rt) / (exp(sqrt(2)*rt) - 1));
% Per-cluster equivalent spread (after truncation)
chann.c_ds_eq_vec  = s_ds  * chann.c_ds_vec;
chann.c_asa_eq_vec = s_asa * chann.c_asa_vec;

% ---- Coherence bandwidth
% Intra-cluster delay statistics under truncated exponential PDP
%   Use the reference (first NLoS) cluster's spread for RMS computation
%   f(x) = (1/C)*exp(-x/c_ds),  0 <= x < rt*c_ds
%   C = c_ds*(1 - exp(-rt))
mu_intra  = chann.c_ds_ref * (1 - (rt+1)*exp(-rt)) / (1 - exp(-rt));
M2_intra  = 2*chann.c_ds_ref^2 * (1 - (0.5*rt^2+rt+1)*exp(-rt)) / (1 - exp(-rt));
var_intra = M2_intra - mu_intra^2;

% Inter-cluster delay dispersion (weighted variance of cluster-center delays)
w_l = chann.cluster_powers_NoBlockage(:);  % normalized power weights
tau_avg = sum(w_l .* chann.tau(:));
var_inter = sum(w_l .* (chann.tau(:) - tau_avg).^2);

% Overall RMS delay spread (Law of Total Variance)
sigma_RMS_sq = chann.sigma_tau^2 + var_intra + var_inter;
chann.sigma_RMS = sqrt(sigma_RMS_sq);

% Coherence bandwidth (50% correlation level)
chann.Bc_50 = 1 / (5 * chann.sigma_RMS);
% Coherence bandwidth (90% correlation level)
chann.Bc_90 = 1 / (50 * chann.sigma_RMS);

fprintf('RMS时延扩展:%.2f ns, 相干带宽(50%%):%.2f MHz, 相干带宽(90%%):%.2f MHz\n', ...
    chann.sigma_RMS*1e9, chann.Bc_50/1e6, chann.Bc_90/1e6)

% ========================= Simulation Parameters =========================
simu.PilotPlacementNum = 1;
simu.MCNum = 1;

% ========================= LMMSE-Estimator Parameters =========================
% Independent switches: approximate likelihood; exact FFT-PCG final solve.
estor.parameter_backend = 'beam_delay_diag';
estor.final_solver = 'fft_pcg';
estor.pilot_only_method = 'dft';
estor.pcg_opts = struct('tol', 1e-8, 'max_iter', 5000, 'verbose', false);
estor.sigma_tau_hat = chann.sigma_tau;
estor.sigma_theta_hat = chann.sigma_theta;
% Assumed uncertainties used by CPM/Fixed-profile/Proposed. Channel generation
% and Oracle use chann.sigma_* (actual); equal defaults preserve matched runs.

estor.ao_opts = struct( ...
    'ao_max_iter',    3, ...
    'ao_tol',         1e-4, ...
    'em_max_iter',    3, ...
    'em_tol',         1e-12, ...
    'lbfgs_max_iter', 2, ...
    'lbfgs_tol',      1e-6, ...
    'lbfgs_m',        5);

% ========================= Export Functions =========================
% Original Monte Carlo entry plus shared building blocks for paper experiments.
funcs.CalcuNMSEsByMonteCarlo = @CalcuNMSEsByMonteCarlo;
funcs.GenerateChannelObservation = @Gen_hfull_hpilot_Observation;
funcs.PilotObservationIndices = @pilot_subcarrier_to_obs_indices;
funcs.BuildClusterCovarianceFactors = @build_cluster_R_factors_hp_yy;
funcs.EstimateChannelFromFactors = @Esti_h_by_yp_ClusterPower;
funcs.EstimateProposedParameters = @EstimateProposedParameters;
funcs.MomentMatching = @ClusterPowersEsti_MomentMatching;
funcs.TraceAO = @AlterOptiClusterPowerAndSpreadFactor;
funcs.TraceEM = @ClusterPowersEsti_SQUAREM;
funcs.ComputeNLL = @compute_nll;
funcs.BuildFFTOperator = @build_bttb_fft_operator;
funcs.ApplyFFTOperator = @apply_bttb_fft_operator;
funcs.DrawTrial = @DrawTrial;
funcs.PilotOnlyEstimate = @PilotOnlyEstimate;
funcs.SafeLineSearch = @SafeLineSearch;
end

function [delay_norm, power_dB, aoa_deg] = get_3gpp_CDL_A()
%GET_3GPP_CDL_A Return the 3GPP CDL-A cluster parameters.
% Inputs:
%   None.
% Outputs:
%   delay_norm - Normalized cluster delays.
%   power_dB   - Cluster powers in dB.
%   aoa_deg    - Legacy receive-angle profile: standard AOD column, retained
%                for the reciprocal uplink convention used in the experiments.

delay_norm = [
    0.0000; 0.3819; 0.4025; 0.5868; 0.4610;
    0.5375; 0.6708; 0.5750; 0.7618; 1.5375;
    1.8978; 2.2242; 2.1718; 2.4942; 2.5119;
    3.0582; 4.0810; 4.4579; 4.5695; 4.7966;
    5.0066; 5.3043; 9.6586];

power_dB = [
    -13.4; 0; -2.2; -4.0; -6.0;
    -8.2; -9.9; -10.5; -7.5; -15.9;
    -6.6; -16.7; -12.4; -15.2; -10.8;
    -11.3; -12.7; -16.2; -18.3; -18.9;
    -16.6; -19.9; -29.7
    ];

aoa_deg = [
    -178.1; -4.2; -4.2; -4.2; 90.2;
    90.2; 90.2; 121.5; -81.7; 158.4;
    -83; 134.8; -153; -172; -129.9;
    -136; 165.4; 148.4; 132.7; -118.6;
    -154.1; 126.5; -56.2
    ];
end

function [delay_norm, power_dB, aoa_deg] = get_3gpp_CDL_D()
%GET_3GPP_CDL_D Return the 3GPP CDL-D cluster parameters.
% Inputs:
%   None.
% Outputs:
%   delay_norm - Normalized cluster delays.
%   power_dB   - Cluster powers in dB.
%   aoa_deg    - Legacy receive-angle profile: standard AOD column, retained
%                for the reciprocal uplink convention used in the experiments.
delay_norm = [
    0; 0; 0.035; 0.612; 1.363;
    1.405; 1.804; 2.596; 1.775; 4.042;
    7.937; 9.424; 9.708; 12.525
    ];
power_dB = [
    -0.2; -13.5; -18.8; -21; -22.8;
    -17.9; -20.1; -21.9; -22.9; -27.8;
    -23.6; -24.8; -30.0; -27.7
    ];
aoa_deg = [
    0; 0; 89.2; 89.2; 89.2;
    13; 13; 13; 34.6; -64.5;
    -32.9; 52.6; -132.1; 77.2
    ];
end

function C = kron_fast(A, B)
%KRON_FAST Compute kron(A,B) using reshape-based implicit expansion.
% Inputs:
%   A, B - Input matrices.
% Output:
%   C    - Kronecker product of A and B.
[M,N] = size(A);
[P,Q] = size(B);
C = reshape(reshape(B,P,1,Q,1) .* reshape(A,1,M,1,N), P*M, Q*N);
end

function obs_idx = pilot_subcarrier_to_obs_indices(pilot_sc, pilot_ant, N)
%PILOT_SUBCARRIER_TO_OBS_INDICES Map pilot locations to stacked-vector indices.
% Inputs:
%   pilot_sc  - Selected subcarrier indices.
%   pilot_ant - Selected antenna indices within each pilot subcarrier.
%   N         - Total number of receive antennas.
% Output:
%   obs_idx   - Linear indices in the antenna-fast stacked channel vector.

Kp_local = numel(pilot_sc);
Np_local = numel(pilot_ant);
obs_idx = zeros(Np_local * Kp_local, 1);
pilot_ant = pilot_ant(:);
for k_idx = 1:Kp_local
    base = (pilot_sc(k_idx) - 1) * N;
    obs_idx((k_idx - 1) * Np_local + (1:Np_local)) = base + pilot_ant;
end
end

function R_thetal_full = gaussian_sin_charfun(AntNum, theta, sigma_theta, c_asa, rt)
%GAUSSIAN_SIN_CHARFUN Build a cluster spatial covariance matrix.
% The Jacobi-Anger expansion averages the array response over Gaussian
% cluster-angle uncertainty and Laplacian intra-cluster angular spread.
% Inputs:
%   AntNum      - Number of array elements.
%   theta       - Prior mean cluster AoA in radians.
%   sigma_theta - Standard deviation of the AoA prior in radians.
%   c_asa       - Intra-cluster angular spread factor in radians.
%   rt          - Optional truncation factor; zero selects infinite tails.
% Output:
%   R_thetal_full - AntNum-by-AntNum normalized spatial covariance matrix.

if nargin < 5 || isempty(rt), rt = 0; end

ant_idx = (0:AntNum-1).';
ant_phase = pi * ant_idx;
max_order = max(50, ceil(max(ant_phase) + 40));
R_thetal_full = zeros(AntNum,AntNum);

for m = 0:max_order
    % Evaluate the angular-spread characteristic function.
    denom_asa = 1 + 0.5 * m^2 * c_asa^2;
    if rt > 0 && c_asa > 0
        % Characteristic function of the truncated Laplacian profile.
        spread_factor = (1 - exp(-sqrt(2)*rt) * (cos(m*rt*c_asa) ...
            - m*c_asa/sqrt(2) * sin(m*rt*c_asa))) ...
            / ((1 - exp(-sqrt(2)*rt)) * denom_asa);
    else
        % Characteristic function of the infinite-tail Laplacian profile.
        spread_factor = 1 / denom_asa;
    end

    c_m = exp(1j * m * theta - 0.5 * (m * sigma_theta)^2) * spread_factor;

    col = besselj(m, ant_phase);
    row = (-1)^m*col;
    Jm = toeplitz(col, row);

    if m == 0
        % m=0: c_{-0} = c_0, so no doubling needed
        R_thetal_full = R_thetal_full + Jm * c_m;
    else
        % Fold m and -m together using c_{-m} = conj(c_m):
        %   J_m(x)*c_m + J_{-m}(x)*c_{-m} = J_m(x)*(c_m + (-1)^m*conj(c_m))
        R_thetal_full = R_thetal_full + Jm * (c_m + (-1)^m * conj(c_m));
    end
end

end

function nll = compute_nll(p, R_factors_pp, sigma2, likelihood_data, d, parameter_backend)
%COMPUTE_NLL Evaluate the Gaussian negative log-likelihood.
% Inputs:
%   p            - Current cluster-power vector.
%   R_factors_pp - Per-cluster pilot-domain covariance factors.
%   sigma2       - Noise variance.
%   likelihood_data - Sample covariance (dense) or DFT power (diagonal).
%   d            - Pilot observation dimension.
%   parameter_backend - 'dense' or 'beam_delay_diag'.
% Output:
%   nll          - log(det(C_yy))+trace(C_yy^{-1}*R_hat_yy).
if nargin < 6, parameter_backend = 'dense'; end
switch parameter_backend
    case 'dense'
        L_loc = numel(p);
        Cyy = sigma2 * eye(d);
        for ll = 1:L_loc
            Cyy = Cyy + p(ll) * kron_fast(R_factors_pp.tau{ll}, R_factors_pp.theta{ll});
        end
        Cyy = (Cyy + Cyy') / 2;
        [Lc, flag] = chol(Cyy, 'lower');
        if flag ~= 0
            nll = inf; return;
        end
        Cinv = Lc' \ (Lc \ eye(d));
        nll = 2*sum(log(diag(Lc))) + ...
            real(sum(Cinv .* likelihood_data.', 'all'));
    case 'beam_delay_diag'
        spectra = build_beam_delay_spectra(R_factors_pp);
        variance = combine_beam_delay_spectra(spectra, p, sigma2);
        nll = real(sum(log(variance) + likelihood_data ./ variance));
    otherwise
        error('Base:UnknownParameterBackend', ...
            'Unsupported parameter_backend: %s', parameter_backend);
end
end

function h_est = Esti_h_by_yp_ClusterPower(y_p, p_est, sigma2, ...
    R_factors_pp, R_factors_hp, final_solver, pcg_opts)
%ESTI_H_BY_YP_CLUSTERPOWER Apply the covariance-based LMMSE estimator.
% Inputs:
%   y_p          - Stacked pilot observation vector.
%   p_est        - Estimated or known cluster powers.
%   sigma2       - Noise variance.
%   R_factors_pp - Per-cluster pilot auto-covariance factors.
%   R_factors_hp - Per-cluster full-to-pilot cross-covariance factors.
%   final_solver - 'dense_chol' or 'fft_pcg'.
%   pcg_opts     - PCG tolerance and maximum iteration count.
% Output:
%   h_est        - Estimated full-band stacked channel vector.
L = length(p_est);
Np = size(R_factors_pp.theta{1}, 1);
Kp = size(R_factors_pp.tau{1}, 1);

switch final_solver
    case 'dense_chol'
        C_yy_est = sigma2 * eye(Np * Kp);
        for l = 1:L
            C_yy_est = C_yy_est + p_est(l) * ...
                kron_fast(R_factors_pp.tau{l}, R_factors_pp.theta{l});
        end
        C_yy_est = (C_yy_est + C_yy_est') / 2;
        [L_chol, chol_flag] = chol(C_yy_est, 'lower');
        if chol_flag ~= 0
            error('Base:CholeskyFailure', 'The LMMSE covariance is not positive definite.');
        end
        x = L_chol' \ (L_chol \ y_p);

    case 'fft_pcg'
        spectra = build_beam_delay_spectra(R_factors_pp);
        variance = combine_beam_delay_spectra(spectra, p_est, sigma2);
        operator_spectrum = build_bttb_fft_operator(p_est, sigma2, R_factors_pp);
        apply_A = @(v) apply_bttb_fft_operator(v, operator_spectrum, Np, Kp);
        apply_Minv = @(v) solve_beam_delay_diagonal(v, variance, Np, Kp);
        x0 = apply_Minv(y_p);
        [x, pcg_flag, relres, pcg_iter] = pcg(apply_A, y_p, ...
            pcg_opts.tol, pcg_opts.max_iter, apply_Minv, [], x0);
        if isfield(pcg_opts, 'verbose') && pcg_opts.verbose
            fprintf('PCG diagnostic: flag=%d, iter=%d, relres=%.3e\n', ...
                pcg_flag, pcg_iter, relres);
        end
        if pcg_flag ~= 0
            error('Base:PCGNotConverged', ...
                'PCG flag=%d after %s iterations (relative residual %.3e).', ...
                pcg_flag, mat2str(pcg_iter), relres);
        end

    otherwise
        error('Base:UnknownFinalSolver', ...
            'Unsupported final_solver: %s', final_solver);
end

% Apply the separable cross-covariance factors without forming C_hy.
X = reshape(x, Np, Kp);
K = size(R_factors_hp.tau{1},1);
N = size(R_factors_hp.theta{1},1);
H_est = zeros(N, K);
for l = 1:L
    H_est = H_est + p_est(l) * (R_factors_hp.theta{l} * X * R_factors_hp.tau{l}.');
end
h_est = H_est(:);
end


function [R_factors_hp,R_factors_yy] = build_cluster_R_factors_hp_yy(is_LoS,N,tau,sigma_tau,c_ds_vec,theta,sigma_theta,c_asa_vec,f_full,pilot_sc,pilot_ant,rt)
%BUILD_CLUSTER_R_FACTORS_HP_YY Build separable covariance factors.
% Inputs:
%   is_LoS       - True for the CDL-D LoS special-cluster treatment.
%   N            - Total number of receive antennas.
%   tau          - Prior mean cluster delays.
%   sigma_tau    - Standard deviation of the delay prior.
%   c_ds_vec     - Per-cluster delay spread factors.
%   theta        - Prior mean cluster AoAs in radians.
%   sigma_theta  - Standard deviation of the AoA prior in radians.
%   c_asa_vec    - Per-cluster angular spread factors in radians.
%   f_full       - Full-band subcarrier-frequency vector.
%   pilot_sc     - Pilot subcarrier indices.
%   pilot_ant    - Pilot antenna indices.
%   rt           - Optional spread truncation factor; default is zero.
% Outputs:
%   R_factors_hp - Full-to-pilot cross-covariance factors.
%   R_factors_yy - Pilot auto-covariance factors.

if nargin < 12 || isempty(rt), rt = 0; end

L = length(theta);
c_ds_vec  = c_ds_vec(:);
c_asa_vec = c_asa_vec(:);

% Build per-cluster spatial correlation matrices using cluster-specific c_asa_l
R_theta_cells = cell(L, 1);
l_start = 1;
if is_LoS && L >= 1
    R_theta_cells{1} = gaussian_sin_charfun(N, theta(1), sigma_theta, 0);
    l_start = 2;
end
for l = l_start:L
    R_theta_cells{l} = gaussian_sin_charfun(N, theta(l), sigma_theta, c_asa_vec(l), rt);
end

f_pilot = f_full(pilot_sc);
R_factors_hp = build_cluster_factors(f_full, f_pilot, tau, sigma_tau, c_ds_vec, pilot_ant, R_theta_cells, 1, is_LoS, rt);
R_factors_yy = build_cluster_factors(f_full, f_pilot, tau, sigma_tau, c_ds_vec, pilot_ant, R_theta_cells, 2, is_LoS, rt);

end

function R_factors = build_cluster_factors(f_full, f_pilot, tau_hat, sigma_tau, c_ds_vec, pilot_ant, R_theta_cells, CrossOrPilotFlag, is_LoS, rt)
%BUILD_CLUSTER_FACTORS Assemble frequency and spatial factors by cluster.
% Inputs:
%   f_full           - Full-band subcarrier-frequency vector.
%   f_pilot          - Pilot-subcarrier frequency vector.
%   tau_hat          - Prior mean cluster delays.
%   sigma_tau        - Standard deviation of the delay prior.
%   c_ds_vec         - Per-cluster delay spread factors.
%   pilot_ant        - Pilot antenna indices.
%   R_theta_cells    - Full-array spatial covariance matrices by cluster.
%   CrossOrPilotFlag - 1 for full-to-pilot factors, 2 for pilot factors.
%   is_LoS           - True for the CDL-D LoS special-cluster treatment.
%   rt               - Optional delay-spread truncation factor.
% Output:
%   R_factors        - Structure containing per-cluster tau/theta factors.

if nargin < 10 || isempty(rt), rt = 0; end

L_local  = numel(tau_hat);
c_ds_vec = c_ds_vec(:);
R_factors.tau   = cell(L_local, 1);
R_factors.theta = cell(L_local, 1);

if CrossOrPilotFlag == 1
    freq_diff = f_full(:) - f_pilot.';
else
    freq_diff = f_pilot - f_pilot.';
end

% First-cluster index for the main spread loop
if is_LoS
    if L_local >= 1
        % Cluster 1: pure phase, no spread
        R_factors.tau{1} = exp(-1j * 2 * pi * freq_diff * tau_hat(1));
    end
    if L_local >= 2
        % Cluster 2: delay spread, no Gaussian smoothing (first-arrival NLoS cluster)
        c2 = c_ds_vec(2);
        jw2 = 1j * 2 * pi * freq_diff * c2;
        denom2 = 1 + jw2;
        if rt > 0 && c2 > 0
            ds2 = (1 - exp(-rt * denom2)) ./ ((1 - exp(-rt)) * denom2);
        else
            ds2 = 1 ./ denom2;
        end
        R_factors.tau{2} = exp(-1j * 2 * pi * freq_diff * tau_hat(2)) .* ds2;
    end
    l_start = 3;
else
    if L_local >= 1
        c1 = c_ds_vec(1);
        jw1 = 1j * 2 * pi * freq_diff * c1;
        denom1 = 1 + jw1;
        if rt > 0 && c1 > 0
            ds1 = (1 - exp(-rt * denom1)) ./ ((1 - exp(-rt)) * denom1);
        else
            ds1 = 1 ./ denom1;
        end
        R_factors.tau{1} = exp(-1j * 2 * pi * freq_diff * tau_hat(1)) .* ds1;
    end
    l_start = 2;
end

gauss_smooth = exp(-0.5 * (2 * pi * sigma_tau * freq_diff).^2);
for l = l_start:L_local
    cl = c_ds_vec(l);
    jwl = 1j * 2 * pi * freq_diff * cl;
    denoml = 1 + jwl;
    if rt > 0 && cl > 0
        dsl = (1 - exp(-rt * denoml)) ./ ((1 - exp(-rt)) * denoml);
    else
        dsl = 1 ./ denoml;
    end
    R_factors.tau{l} = exp(-1j * 2 * pi * freq_diff * tau_hat(l)) ...
        .* gauss_smooth .* dsl;
end

for l = 1:L_local
    if CrossOrPilotFlag == 1
        R_factors.theta{l} = R_theta_cells{l}(:, pilot_ant);
    else
        R_factors.theta{l} = R_theta_cells{l}(pilot_ant, pilot_ant);
    end
end
end

function [h_true, y_p] = Gen_hfull_hpilot_Observation(theta,tau,c_asa_vec,c_ds_vec,M_sub,cluster_powers,sigma2,f_full,N,pilot_sc,pilot_ant,SpreadTruncFactor,is_LoS)
%GEN_HFULL_HPILOT_OBSERVATION Generate one channel and pilot observation.
% Inputs:
%   theta, tau        - True cluster AoAs and delays.
%   c_asa_vec         - Per-cluster angular spread factors.
%   c_ds_vec          - Per-cluster delay spread factors.
%   M_sub             - Number of subpaths per diffuse cluster.
%   cluster_powers    - Instantaneous cluster powers.
%   sigma2            - Pilot-noise variance.
%   f_full            - Full-band subcarrier-frequency vector.
%   N                 - Number of receive antennas.
%   pilot_sc          - Pilot subcarrier indices.
%   pilot_ant         - Pilot antenna indices.
%   SpreadTruncFactor - Delay/angle spread truncation factor.
%   is_LoS            - True when the first cluster is specular LoS.
% Outputs:
%   h_true            - Full-band stacked channel vector.
%   y_p               - Noisy stacked pilot observation vector.

L = length(cluster_powers);
c_asa_vec = c_asa_vec(:);
c_ds_vec  = c_ds_vec(:);
h_true = zeros(length(f_full)*N, 1);
l_start = 1;
if is_LoS
    % Generate the deterministic-geometry specular component.
    a_vec = exp(1j * pi * (0:N-1).' * sin(theta(1)));
    f_vec = exp(-1j * 2 * pi * f_full * tau(1));
    alpha_sub = sqrt(cluster_powers(1)) * exp(1j*2*pi*rand(1));
    h_true = h_true + alpha_sub * reshape(a_vec * f_vec.', [], 1);
    l_start = 2;
end
for l = l_start:L
    % Sample subpath locations on the truncated delay-angle support.
    c_ds_l  = c_ds_vec(l);
    c_asa_l = c_asa_vec(l);
    delta_tau   = SpreadTruncFactor*c_ds_l  * rand(M_sub, 1);
    delta_tau   = delta_tau - min(delta_tau);
    delta_theta = 2*SpreadTruncFactor*c_asa_l * (rand(M_sub, 1) - 0.5);
    if c_ds_l > 0 && c_asa_l > 0
        SubPathPowerProps = exp(-delta_tau/c_ds_l - sqrt(2)*abs(delta_theta)/c_asa_l);
    else
        SubPathPowerProps = ones(M_sub, 1);
    end
    % Normalize the subpath powers to the prescribed cluster power.
    alpha_sub = sqrt(SubPathPowerProps/sum(SubPathPowerProps)*cluster_powers(l)).*exp(1j*2*pi*rand(M_sub,1));
    for m = 1:M_sub
        a_vec = exp(1j * pi * (0:N-1).' * sin(theta(l) + delta_theta(m)));
        f_vec = exp(-1j * 2 * pi * f_full * (tau(l) + delta_tau(m)));
        h_true = h_true + alpha_sub(m) * reshape(a_vec * f_vec.', [], 1);
    end
end

% Extract the selected pilot entries and add complex Gaussian noise.
pilot_obs_idx = pilot_subcarrier_to_obs_indices(pilot_sc, pilot_ant, N);
d_pilot = numel(pilot_obs_idx);
n_p = sqrt(sigma2/2) * (randn(d_pilot, 1) + 1j * randn(d_pilot, 1));
y_p = h_true(pilot_obs_idx) + n_p;

end

function ypSamps = Gen_hpilot_Hist(theta,tau,sigma_tau,sigma_theta,c_asa_vec,c_ds_vec,M_sub,cluster_powers_NoBlockage,sigma2,f_full,N,pilot_sc,pilot_ant,SpreadTruncFactor,is_LoS,SampNum,BlockProb)
%GEN_HPILOT_HIST Generate independent historical pilot snapshots.
% Inputs:
%   theta, tau               - Nominal cluster AoAs and delays.
%   sigma_tau, sigma_theta   - Delay and AoA variation standard deviations.
%   c_asa_vec, c_ds_vec      - Per-cluster angular and delay spread factors.
%   M_sub                    - Number of subpaths per diffuse cluster.
%   cluster_powers_NoBlockage - Unblocked cluster-power vector.
%   sigma2                   - Pilot-noise variance.
%   f_full                   - Full-band subcarrier-frequency vector.
%   N                        - Number of receive antennas.
%   pilot_sc, pilot_ant      - Pilot subcarrier and antenna indices.
%   SpreadTruncFactor        - Delay/angle spread truncation factor.
%   is_LoS                   - True when the first cluster is specular LoS.
%   SampNum                  - Number of historical snapshots.
%   BlockProb                - Independent blockage probability per cluster.
% Output:
%   ypSamps                  - Pilot observations, one snapshot per column.

L = length(cluster_powers_NoBlockage);
c_asa_vec = c_asa_vec(:);
c_ds_vec  = c_ds_vec(:);
pilot_obs_idx = pilot_subcarrier_to_obs_indices(pilot_sc, pilot_ant, N);
d_pilot = numel(pilot_obs_idx);

ypSamps = zeros(d_pilot, SampNum);
f_p = f_full(pilot_sc);             % pilot frequencies only
ant_p = (pilot_ant(:) - 1);         % 0-indexed pilot antenna indices
for s = 1:SampNum
    % Draw independent blockage states and cluster-center perturbations.
    BlockedIdxs = rand(L,1) < BlockProb;
    cluster_powers = cluster_powers_NoBlockage;
    cluster_powers(BlockedIdxs) = cluster_powers(BlockedIdxs) * 0.001;

    theta_s = theta + sigma_theta * randn(L, 1);
    if is_LoS
        tau_noise = sigma_tau * randn(L, 1);
        tau_noise(1) = 0;
        if L >= 2, tau_noise(2) = 0; end
        tau_s = tau + tau_noise;
    else
        tau_s = tau + [0; sigma_tau * randn(L-1, 1)];
    end

    h_pilot_s = zeros(d_pilot, 1);
    ls = 1;
    if is_LoS
        % Add the specular LoS contribution directly in the pilot domain.
        a_s = exp(1j * pi * ant_p * sin(theta_s(1)));          % Np x 1
        f_s = exp(-1j * 2 * pi * f_p * tau_s(1));              % Kp x 1
        al_s = sqrt(cluster_powers(1)) * exp(1j*2*pi*rand(1));
        h_pilot_s = h_pilot_s + al_s * reshape(a_s * f_s.', [], 1);
        ls = 2;
    end
    for l = ls:L
        % Generate each diffuse cluster on the truncated delay-angle support.
        c_ds_l  = c_ds_vec(l);
        c_asa_l = c_asa_vec(l);
        dt  = SpreadTruncFactor*c_ds_l  * rand(M_sub, 1);
        dt  = dt - min(dt);
        dth = 2*SpreadTruncFactor*c_asa_l * (rand(M_sub, 1) - 0.5);
        if c_ds_l > 0 && c_asa_l > 0
            spp = exp(-dt/c_ds_l - sqrt(2)*abs(dth)/c_asa_l);
        else
            spp = ones(M_sub, 1);
        end
        % Normalize subpath powers and accumulate the pilot-domain channel.
        al = sqrt(spp/sum(spp)*cluster_powers(l)) .* exp(1j*2*pi*rand(M_sub,1));
        for m = 1:M_sub
            a_s = exp(1j * pi * ant_p * sin(theta_s(l) + dth(m)));
            f_s = exp(-1j * 2 * pi * f_p * (tau_s(l) + dt(m)));
            h_pilot_s = h_pilot_s + al(m) * reshape(a_s * f_s.', [], 1);
        end
    end
    n_s = sqrt(sigma2/2) * (randn(d_pilot, 1) + 1j * randn(d_pilot, 1));
    ypSamps(:, s) = h_pilot_s + n_s;
end

end

function [NMSE_Oracle_equa,NMSE_sigmaEqu0_cEqu0,NMSE_cEqu0,NMSE_emp, ...
    NMSE_esti,NMSE_Oracle,runtime_parameter,runtime_channel,NMSE_PilotOnly,trial_stats] = CalcuNMSEsByMonteCarlo( ...
    sys, chann, simu, estor, pilot_sc, pilot_ant, S_hist, MCNMSEsSwitch)
%CALCUNMSESBYMONTECARLO Evaluate selected estimators by Monte Carlo trials.
% Inputs:
%   sys              - OFDM and antenna-array parameters.
%   chann            - Channel, prior, spread, SNR, and blockage parameters.
%   simu             - Simulation settings including the trial count.
%   estor            - Estimator optimization settings.
%   pilot_sc         - Pilot subcarrier indices.
%   pilot_ant        - Pilot antenna indices.
%   S_hist           - Historical snapshot count; zero selects online mode.
%   MCNMSEsSwitch    - Up-to-seven-element logical vector selecting estimators.
%                      Element 7 enables LS plus estor.pilot_only_method.
% Outputs:
%   NMSE_Oracle_equa     - NMSE with true powers and equivalent spreads.
%   NMSE_sigmaEqu0_cEqu0 - NMSE without prior uncertainty or spread.
%   NMSE_cEqu0           - CPM-aided NMSE with zero intra-cluster spread.
%   NMSE_emp             - NMSE using empirical powers and nominal spreads.
%   NMSE_esti            - Proposed estimator NMSE.
%   NMSE_Oracle          - Oracle NMSE using truncated covariance models.
%   runtime_parameter    - Proposed parameter-estimation runtime in seconds.
%   runtime_channel      - Final channel-estimation runtime in seconds.
%   NMSE_PilotOnly       - Pilot-only LS plus selected-interpolation NMSE.
%   trial_stats          - Per-realization error energies and channel energies.

sigma2 = 10^(-chann.SNR_dB / 10);
L = length(chann.tau);
MCNMSEsSwitch = logical(MCNMSEsSwitch(:).');
if numel(MCNMSEsSwitch) < 7
    MCNMSEsSwitch(7) = false;
end

parameter_backend = estor.parameter_backend;
final_solver = estor.final_solver;

% Accumulate error and channel energy before normalization to reduce
% finite-sample variability in the Monte Carlo NMSE ratio.
MSE_Oracle_equa = 0; MSE_sigmaEqu0_cEqu0 = 0; MSE_cEqu0 = 0;
MSE_emp = 0; MSE_esti = 0; MSE_Oracle = 0; MSE_PilotOnly = 0;
h_power_total = 0;
trial_stats.error = nan(simu.MCNum,7);
trial_stats.energy = zeros(simu.MCNum,1);
runtime_parameter = 0;
runtime_channel = 0;

fprintf('开始MonteCarlo仿真...... \n')
total_timer = tic;
for mc = 1:simu.MCNum

    % Draw independent per-cluster blockage states.
    BlockedIdxs = rand(L,1)<chann.ClusterBlockProb;
    cluster_powers_WithBlockage = chann.cluster_powers_NoBlockage;
    cluster_powers_WithBlockage(BlockedIdxs) = cluster_powers_WithBlockage(BlockedIdxs)*0.001;

    % Generate the true channel and its noisy pilot observation.
    [h_true, y_p] = Gen_hfull_hpilot_Observation( ...
        chann.theta, chann.tau, chann.c_asa_vec, chann.c_ds_vec, chann.M_subpath, ...
        cluster_powers_WithBlockage, sigma2, ...
        sys.f, sys.N, pilot_sc, pilot_ant, chann.SpreadTruncFactor, ...
        chann.is_LoS);

    h_power_mc = sum(abs(h_true).^2);
    h_power_total = h_power_total + h_power_mc;
    trial_stats.energy(mc) = h_power_mc;

    if MCNMSEsSwitch(7)  % Pilot-only LS with the selected interpolation.
        h_est = PilotOnlyEstimate(y_p, sys.K, sys.N, pilot_sc, estor.pilot_only_method);
        trial_stats.error(mc,7) = sum(abs(h_true - h_est).^2);
        MSE_PilotOnly = MSE_PilotOnly + trial_stats.error(mc,7);
    end

    % Perturb the nominal cluster parameters to obtain the CCM priors.
    theta_hat = chann.theta + chann.sigma_theta * randn(L, 1);
    if chann.is_LoS
        tau_noise = chann.sigma_tau * randn(L, 1);
        tau_noise(1) = 0;
        if L >= 2, tau_noise(2) = 0; end
        tau_hat = chann.tau + tau_noise;
    else
        tau_hat = chann.tau + [0; chann.sigma_tau * randn(L-1, 1)];
    end

    if MCNMSEsSwitch(1)  % Oracle-equal spread factor (per-cluster truncated equivalent)
        [R_factors_hp, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
            chann.is_LoS, sys.N, tau_hat, chann.sigma_tau, chann.c_ds_eq_vec, ...
            theta_hat, chann.sigma_theta, chann.c_asa_eq_vec, sys.f, pilot_sc, pilot_ant);
        h_est = Esti_h_by_yp_ClusterPower(y_p, cluster_powers_WithBlockage, ...
            sigma2, R_factors_pp, R_factors_hp, final_solver, estor.pcg_opts);
        trial_stats.error(mc,1) = sum(abs(h_true - h_est).^2);
        MSE_Oracle_equa = MSE_Oracle_equa + trial_stats.error(mc,1);
    end

    if MCNMSEsSwitch(2)  % sigma=0, c=0
        [R_factors_hp, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
            chann.is_LoS, sys.N, tau_hat, 0, zeros(L,1), ...
            theta_hat, 0, zeros(L,1), sys.f, pilot_sc, pilot_ant);
        h_est = Esti_h_by_yp_ClusterPower(y_p, cluster_powers_WithBlockage, ...
            sigma2, R_factors_pp, R_factors_hp, final_solver, estor.pcg_opts);
        trial_stats.error(mc,2) = sum(abs(h_true - h_est).^2);
        MSE_sigmaEqu0_cEqu0 = MSE_sigmaEqu0_cEqu0 + trial_stats.error(mc,2);
    end

    if MCNMSEsSwitch(3)  % sigma!=0, c=0
        [R_factors_hp, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
            chann.is_LoS, sys.N, tau_hat, estor.sigma_tau_hat, zeros(L,1), ...
            theta_hat, estor.sigma_theta_hat, zeros(L,1), sys.f, pilot_sc, pilot_ant);
        h_est = Esti_h_by_yp_ClusterPower(y_p, cluster_powers_WithBlockage, ...
            sigma2, R_factors_pp, R_factors_hp, final_solver, estor.pcg_opts);
        trial_stats.error(mc,3) = sum(abs(h_true - h_est).^2);
        MSE_cEqu0 = MSE_cEqu0 + trial_stats.error(mc,3);
    end

    if MCNMSEsSwitch(4)  % Empirical delay-profile power; spread model selected below.
        baseline_c_ds = chann.c_ds_ref * ones(L,1);
        baseline_c_asa = chann.c_asa_ref * ones(L,1);
        if chann.is_LoS
            baseline_c_ds(1) = 0;
            baseline_c_asa(1) = 0;
        end
        [R_factors_hp, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
            chann.is_LoS, sys.N, tau_hat, estor.sigma_tau_hat, baseline_c_ds, ...
            theta_hat, estor.sigma_theta_hat, baseline_c_asa, sys.f, pilot_sc, pilot_ant);
        h_est = Esti_h_by_yp_ClusterPower(y_p, chann.cluster_powers_emp, ...
            sigma2, R_factors_pp, R_factors_hp, final_solver, estor.pcg_opts);
        trial_stats.error(mc,4) = sum(abs(h_true - h_est).^2);
        MSE_emp = MSE_emp + trial_stats.error(mc,4);
    end

    if MCNMSEsSwitch(5)  % AO estimated power + spread
        parameter_timer = tic;
        % Generate an independent historical data set for each trial.
        if S_hist == 0
            ypUsed = y_p;
        else
            ypUsed = Gen_hpilot_Hist( ...
                chann.theta, chann.tau, chann.sigma_tau, chann.sigma_theta, ...
                chann.c_asa_vec, chann.c_ds_vec, chann.M_subpath, ...
                chann.cluster_powers_NoBlockage, sigma2, ...
                sys.f, sys.N, pilot_sc, pilot_ant, ...
                chann.SpreadTruncFactor, chann.is_LoS, S_hist, chann.ClusterBlockProb);
        end
        [p_hat, c_ds_vec_hat, c_asa_vec_hat] = EstimateProposedParameters( ...
            ypUsed, sigma2, sys, chann, estor, tau_hat, theta_hat, ...
            pilot_sc, pilot_ant, parameter_backend);
        runtime_parameter = runtime_parameter + toc(parameter_timer);

        channel_timer = tic;
        [R_factors_hp, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
            chann.is_LoS, sys.N, tau_hat, estor.sigma_tau_hat, c_ds_vec_hat, ...
            theta_hat, estor.sigma_theta_hat, c_asa_vec_hat, sys.f, pilot_sc, pilot_ant);
        h_est = Esti_h_by_yp_ClusterPower(y_p, p_hat, sigma2, ...
            R_factors_pp, R_factors_hp, final_solver, estor.pcg_opts);
        trial_stats.error(mc,5) = sum(abs(h_true - h_est).^2);
        MSE_esti = MSE_esti + trial_stats.error(mc,5);
        runtime_channel = runtime_channel + toc(channel_timer);
    end

    if MCNMSEsSwitch(6) % Oracle (truncated CF with true per-cluster spread)
        [R_factors_hp, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
            chann.is_LoS, sys.N, tau_hat, chann.sigma_tau, chann.c_ds_vec, ...
            theta_hat, chann.sigma_theta, chann.c_asa_vec, sys.f, pilot_sc, pilot_ant, chann.SpreadTruncFactor);
        h_est = Esti_h_by_yp_ClusterPower(y_p, cluster_powers_WithBlockage, ...
            sigma2, R_factors_pp, R_factors_hp, final_solver, estor.pcg_opts);
        trial_stats.error(mc,6) = sum(abs(h_true - h_est).^2);
        MSE_Oracle = MSE_Oracle + trial_stats.error(mc,6);
    end

    fprintf('mc_idx = %d/%d, 耗时 %.2fs \n', mc, simu.MCNum, toc(total_timer));

end

% Normalize after all trials: NMSE=sum||h-h_est||^2/sum||h||^2.
NMSE_Oracle_equa     = MSE_Oracle_equa     / h_power_total;
NMSE_sigmaEqu0_cEqu0 = MSE_sigmaEqu0_cEqu0 / h_power_total;
NMSE_cEqu0           = MSE_cEqu0           / h_power_total;
NMSE_emp             = MSE_emp             / h_power_total;
NMSE_esti            = MSE_esti            / h_power_total;
NMSE_Oracle          = MSE_Oracle          / h_power_total;
NMSE_PilotOnly       = MSE_PilotOnly       / h_power_total;

end

function cluster_powers_esti = ClusterPowersEsti_MomentMatching(y_p, R_factors_pp, sigma2)
%CLUSTERPOWERSESTI_MOMENTMATCHING Initialize cluster powers by moments.
% Inputs:
%   y_p          - Pilot observations of size d-by-S.
%   R_factors_pp - Per-cluster pilot covariance factors.
%   sigma2       - Noise variance.
% Output:
%   cluster_powers_esti - Nonnegative moment-matching power estimates.

L = numel(R_factors_pp.tau);
R_trace = zeros(L, 1);
R_frob2 = zeros(L, 1);
Np = size(R_factors_pp.theta{1}, 1);
Kp = size(R_factors_pp.tau{1}, 1);
S  = size(y_p, 2);
Y_3d = reshape(y_p, Np, Kp, S);  % Np x Kp x S

for l = 1:L
    R_trace(l) = real(trace(R_factors_pp.theta{l})) * real(trace(R_factors_pp.tau{l}));
    R_frob2(l) = real(sum(abs(R_factors_pp.theta{l}(:)).^2)) * real(sum(abs(R_factors_pp.tau{l}(:)).^2));
end
p_est = zeros(L, 1);
for l = 1:L
    temp = pagemtimes(R_factors_pp.theta{l}, Y_3d);
    temp = pagemtimes(temp, R_factors_pp.tau{l}.');
    quad_term = real(sum(conj(Y_3d) .* temp, 'all')) / S;
    p_est(l) = max(quad_term - sigma2*R_trace(l), 0) / max(R_frob2(l), eps);
end
cluster_powers_esti = max(p_est, 1e-12);

end

function [p_hat, obj_history, p_history] = ClusterPowersEsti_SQUAREM(y_p, R_factors_pp, sigma2, p_init, opts)
%CLUSTERPOWERSESTI_SQUAREM Estimate cluster powers using accelerated EM.
% Inputs:
%   y_p          - Pilot observations of size d-by-S.
%   R_factors_pp - Per-cluster pilot covariance factors.
%   sigma2       - Noise variance.
%   p_init       - Initial cluster-power vector.
%   opts         - Structure with max_iter and tol.
% Outputs:
%   p_hat        - SQUAREM-accelerated cluster-power estimate.
%   obj_history  - Evaluated negative log-likelihood values.
%   p_history    - Powers at evaluated objectives (including rejected trials).

L_local = numel(p_init);
Np = size(R_factors_pp.theta{1}, 1);
Kp = size(R_factors_pp.tau{1}, 1);
d = Np * Kp;
S = size(y_p, 2);

if ~isfield(opts, 'parameter_backend'), opts.parameter_backend = 'dense'; end
parameter_backend = opts.parameter_backend;

% Precompute the representation reused by all SQUAREM fixed-point steps.
if strcmp(parameter_backend, 'dense')
    eye_vec = reshape(eye(d), d*d, 1);
    R_mat = zeros(d*d, L_local);
    for l = 1:L_local
        temp = kron_fast(R_factors_pp.tau{l}, R_factors_pp.theta{l});
        R_mat(:, l) = temp(:);
    end
    spectra = [];
    sample_power = [];
elseif strcmp(parameter_backend, 'beam_delay_diag')
    eye_vec = [];
    R_mat = [];
    spectra = build_beam_delay_spectra(R_factors_pp);
    sample_power = beam_delay_sample_power(y_p, Np, Kp);
else
    error('Base:UnknownParameterBackend', ...
        'Unsupported parameter_backend: %s', parameter_backend);
end

    function [p_next, obj_curr] = em_update_with_obj(p_curr)
        %EM_UPDATE_WITH_OBJ Perform one EM update and evaluate its objective.
        % Input:
        %   p_curr  - Current cluster-power vector.
        % Outputs:
        %   p_next  - Cluster-power vector after one EM update.
        %   obj_curr - Negative log-likelihood at p_curr.
        if strcmp(parameter_backend, 'dense')
            C_yy_vec = R_mat * p_curr + sigma2 * eye_vec;
            C_yy = reshape(C_yy_vec, d, d);
            C_yy = (C_yy + C_yy') / 2;

            [L_chol, chol_flag] = chol(C_yy, 'lower');
            if chol_flag ~= 0
                obj_curr = inf;
                p_next = p_curr;
                return;
            end

            C_inv_Y = L_chol' \ (L_chol \ y_p);   % d x S
            C_inv   = L_chol' \ (L_chol \ eye(d));

            obj_curr = 2 * sum(log(diag(L_chol))) ...
                + real(sum(conj(C_inv_Y) .* y_p, 'all')) / S;

            tr_terms = real(C_inv(:)' * R_mat).';
            Z_3d = reshape(C_inv_Y, Np, Kp, S);
            quad_terms = zeros(L_local, 1);
            for ll = 1:L_local
                tmp = pagemtimes(R_factors_pp.theta{ll}, Z_3d);
                tmp = pagemtimes(tmp, R_factors_pp.tau{ll}.');
                quad_terms(ll) = real(sum(conj(Z_3d) .* tmp, 'all')) / S;
            end
        else
            variance = combine_beam_delay_spectra(spectra, p_curr, sigma2);
            inv_variance = 1 ./ variance;
            obj_curr = real(sum(log(variance) + sample_power .* inv_variance));
            tr_terms = real(spectra.' * inv_variance);
            quad_terms = real(spectra.' * ...
                (sample_power .* inv_variance.^2));
        end

        p_next = p_curr + (p_curr.^2 / d) .* (quad_terms - tr_terms);
        p_next = max(p_next, 1e-12);
    end

% SQUAREM extrapolation loop with reused EM objective evaluations.
p = p_init(:);

obj_history = zeros(opts.max_iter * 2 + 1, 1);
p_history = zeros(L_local, opts.max_iter * 2 + 1);

[p1, obj_p] = em_update_with_obj(p);

cost_idx = 1;
obj_history(cost_idx) = obj_p;
p_history(:, cost_idx) = p;

for iter = 1:opts.max_iter
    [p2, obj_p1] = em_update_with_obj(p1);

    cost_idx = cost_idx + 1;
    obj_history(cost_idx) = obj_p1;
    p_history(:, cost_idx) = p1;

    r = p1 - p;
    v = (p2 - p1) - r;

    if norm(v) < 1e-12
        % Skip extrapolation when the second-order direction vanishes.
        p = p1;
        p1 = p2;
    else
        % Compute the SQUAREM extrapolation step.
        alpha = -norm(r) / norm(v);

        % Form a nonnegative extrapolated power vector.
        p_ext = p - 2 * alpha * r + (alpha^2) * v;
        p_ext = max(p_ext, 1e-12);

        % Evaluate the extrapolated point and cache its next EM iterate.
        [p_ext_next, obj_ext] = em_update_with_obj(p_ext);

        cost_idx = cost_idx + 1;
        obj_history(cost_idx) = obj_ext;
        p_history(:, cost_idx) = p_ext;

        % Reject extrapolation if it increases the negative log-likelihood.
        if obj_ext > obj_p1
            p = p1;
            p1 = p2;
        else
            % Accept the safeguarded extrapolation.
            p = p_ext;
            p1 = p_ext_next;
        end
    end

    % Use the relative EM fixed-point residual as the stopping criterion.
    if norm(p1 - p) / (norm(p) + eps) < opts.tol
        p = p1;
        break;
    end
end

if length(obj_history) > cost_idx
    obj_history = obj_history(1:cost_idx);
    p_history = p_history(:, 1:cost_idx);
end

p_hat = p;
end

function [c_ds_vec_hat, c_asa_vec_hat, obj_history, c_history] = SpreadFactorEsti(y_p, p_est, sigma2, c_ds_vec_init, c_asa_vec_init, tau_hat, sigma_tau, sigma_theta, theta_hat, f_pilot, pilot_ant, N, is_LoS, opts)
%SPREADFACTORESTI Jointly estimate per-cluster delay and angular spreads.
% Uses log-domain L-BFGS with a strong-Wolfe line search to optimize all
% {c_ds,l,c_asa,l} while the cluster powers remain fixed.
%
% Inputs:
%   y_p              - Pilot observations of size d-by-S.
%   p_est            - Fixed cluster-power estimates.
%   sigma2           - Noise variance.
%   c_ds_vec_init    - Initial delay spread factors.
%   c_asa_vec_init   - Initial angular spread factors.
%   tau_hat          - Prior mean cluster delays.
%   sigma_tau        - Standard deviation of the delay prior.
%   sigma_theta      - Standard deviation of the AoA prior.
%   theta_hat        - Prior mean cluster AoAs in radians.
%   f_pilot          - Pilot-subcarrier frequency vector.
%   pilot_ant        - Pilot antenna indices.
%   N                - Total number of receive antennas.
%   is_LoS           - True for the CDL-D LoS special-cluster treatment.
%   opts             - Structure with max_iter, tol, and lbfgs_m.
%
% Outputs:
%   c_ds_vec_hat     - Estimated delay spread factors.
%   c_asa_vec_hat    - Estimated angular spread factors.
%   obj_history      - Negative log-likelihood history.
%   c_history        - Interleaved spread-factor history by iteration.

% Apply default optimizer settings when fields are omitted.
if ~isfield(opts, 'max_iter'), opts.max_iter = 50;  end
if ~isfield(opts, 'tol'),      opts.tol      = 1e-6; end
if ~isfield(opts, 'lbfgs_m'),  opts.lbfgs_m  = 5;    end
if ~isfield(opts, 'parameter_backend'), opts.parameter_backend = 'dense'; end
parameter_backend = opts.parameter_backend;

L  = numel(p_est);
Kp = numel(f_pilot);
Np = numel(pilot_ant);
d  = Kp * Np;
m_mem = opts.lbfgs_m;

c_ds_vec_init  = c_ds_vec_init(:);
c_asa_vec_init = c_asa_vec_init(:);

% Cache only the observation statistic required by the selected backend.
if strcmp(parameter_backend, 'dense')
    likelihood_data = (y_p * y_p') / size(y_p, 2);
elseif strcmp(parameter_backend, 'beam_delay_diag')
    likelihood_data = beam_delay_sample_power(y_p, Np, Kp);
else
    error('Base:UnknownParameterBackend', ...
        'Unsupported parameter_backend: %s', parameter_backend);
end

% Pilot frequency difference matrix (Kp x Kp)
freq_diff = f_pilot(:) - f_pilot(:).';

% Precompute the iteration-invariant Bessel basis.
ant_phase = pi * (0:N-1)';
max_order = max(50, ceil(max(ant_phase) + 40));
m_vals = (0:max_order)';
Jm_pilot = zeros(Np, Np, max_order+1);  % Np x Np x (M+1)
for mm = 0:max_order
    col = besselj(mm, ant_phase);
    Jfull = toeplitz(col, (-1)^mm * col);
    Jm_pilot(:,:,mm+1) = Jfull(pilot_ant, pilot_ant);
end
% Precompute phase weights determined by the AoA prior.
basis_w = zeros(L, max_order+1);  % complex
for li = 1:L
    for mm = 0:max_order
        cph = exp(1j*mm*theta_hat(li) - 0.5*(mm*sigma_theta)^2);
        if mm == 0, basis_w(li,mm+1) = cph;
        else,       basis_w(li,mm+1) = cph + (-1)^mm * conj(cph); end
    end
end
% LoS cluster 1 has c_asa=0; its R_a is fixed
if is_LoS && L>=1
    Ra_LoS = sum(Jm_pilot .* reshape(basis_w(1,:),1,1,[]), 3);
    dRa_LoS = zeros(Np);
end

% Precompute frequency-domain phase and prior-smoothing terms.
phase_f = cell(L,1);
gauss_smooth = exp(-0.5*(2*pi*sigma_tau*freq_diff).^2);
for li = 1:L
    ph = exp(-1j*2*pi*freq_diff*tau_hat(li));
    if is_LoS && li == 1
        phase_f{li} = ph;           % LoS: pure phase, no c_ds dependence
    elseif (is_LoS && li == 2) || (~is_LoS && li == 1)
        phase_f{li} = ph;           % first NLoS: no Gaussian smoothing
    else
        phase_f{li} = ph .* gauss_smooth;
    end
end

% Initialize interleaved log-spreads:
% c_tilde(2l-1)=log(c_ds_l), c_tilde(2l)=log(c_asa_l).
bounded = isfield(opts, 'spread_bounds');
if bounded
    log_bounds = log(opts.spread_bounds); % Rows: delay/angle; columns: lower/upper.
    lower_log = repmat(log_bounds(:,1), L, 1);
    upper_log = repmat(log_bounds(:,2), L, 1);
end
c_tilde = zeros(2*L, 1);
for li = 1:L
    c_tilde(2*li-1) = log(max(c_ds_vec_init(li),  1e-15));
    c_tilde(2*li)   = log(max(c_asa_vec_init(li), 1e-15));
end
% LoS cluster 1 c_ds=0 and c_asa=0 are fixed; keep them out of optimization
% by clamping their gradient to zero in obj_grad

% L-BFGS curvature storage (2L-dimensional)
S_store = zeros(2*L, 0);
Y_store = zeros(2*L, 0);
rho_store = [];

obj_history = zeros(opts.max_iter + 1, 1);
c_history   = zeros(2*L, opts.max_iter + 1);
[obj, grad] = obj_grad(c_tilde);
assert(isfinite(obj) && all(isfinite(grad)), 'Base:InvalidInitialObjective', ...
    'Initial spread objective and gradient must be finite.');
obj_history(1) = obj;
c_history(:, 1) = exp(c_tilde);

for iter = 1:opts.max_iter
    if norm(grad) < opts.tol || all(grad == 0)
        obj_history = obj_history(1:iter);
        c_history = c_history(:, 1:iter);
        break;
    end
    % --- L-BFGS two-loop recursion ---
    dir = lbfgs_dir(grad, S_store, Y_store, rho_store);

    % Strong Wolfe when available; otherwise verified Armijo descent only.
    [alpha, obj_new, grad_new, dir, ls_status, reset] = ...
        SafeLineSearch(@obj_grad, c_tilde, dir, obj, grad);
    if reset || strcmp(ls_status, 'armijo')
        S_store = zeros(2*L,0); Y_store = zeros(2*L,0); rho_store = [];
    end
    if alpha == 0
        obj_history = obj_history(1:iter);
        c_history = c_history(:, 1:iter);
        warning('Base:LineSearchFailed', ...
            'No finite sufficient-decrease step found; spread parameters unchanged.');
        break;
    end

    c_tilde_new = c_tilde + alpha * dir;
    s_k = c_tilde_new - c_tilde;
    y_k = grad_new - grad;
    sy  = s_k' * y_k;

    if isfinite(sy) && sy > 1e-15
        if size(S_store, 2) >= m_mem
            S_store   = S_store(:, 2:end);
            Y_store   = Y_store(:, 2:end);
            rho_store = rho_store(2:end);
        end
        S_store   = [S_store, s_k]; %#ok<AGROW>
        Y_store   = [Y_store, y_k]; %#ok<AGROW>
        rho_store = [rho_store; 1/sy]; %#ok<AGROW>
    end

    c_tilde = c_tilde_new;
    obj  = obj_new;
    grad = grad_new;
    obj_history(iter + 1) = obj;
    c_history(:, iter + 1) = exp(c_tilde);

    if norm(grad) < opts.tol
        obj_history = obj_history(1:iter+1);
        c_history   = c_history(:, 1:iter+1);
        break;
    end
    if iter == opts.max_iter
        obj_history = obj_history(1:iter+1);
        c_history   = c_history(:, 1:iter+1);
    end
end

% Extract per-cluster estimates
c_ds_vec_hat  = zeros(L, 1);
c_asa_vec_hat = zeros(L, 1);
for li = 1:L
    c_ds_vec_hat(li)  = exp(c_tilde(2*li-1));
    c_asa_vec_hat(li) = exp(c_tilde(2*li));
end

    function [fval, g] = obj_grad(ct)
        %OBJ_GRAD Evaluate the spread objective and analytical gradient.
        % Input:
        %   ct   - Interleaved log-domain delay and angular spreads.
        % Outputs:
        %   fval - Gaussian negative log-likelihood.
        %   g    - Gradient with respect to ct.

        % Optional experiment-only box: reject infeasible steps, never clip results.
        if bounded && any(ct < lower_log | ct > upper_log)
            fval = inf; g = zeros(size(ct)); return;
        end
        % Extract per-cluster spread factors
        cds_vec  = exp(ct(1:2:end));   % L x 1
        casa_vec = exp(ct(2:2:end));   % L x 1
        if any(~isfinite(ct)) || any(~isfinite(cds_vec)) || any(~isfinite(casa_vec))
            fval = inf; g = zeros(2*L,1); return;
        end

        jfreq = -1j*2*pi*freq_diff;

        % --- Freq-domain R_f, dR_f (per-cluster, each uses its own c_ds_l) ---
        Rf  = cell(L,1);
        dRf = cell(L,1);
        if is_LoS && L>=1
            Rf{1}  = phase_f{1};
            dRf{1} = zeros(Kp);   % LoS: c_ds fixed at 0
        end
        ls_f = 1 + is_LoS;
        for ll = ls_f:L
            inv_denom_ll = 1 ./ (1 + 1j*2*pi*freq_diff*cds_vec(ll));
            Rf{ll}  = phase_f{ll} .* inv_denom_ll;
            dRf{ll} = phase_f{ll} .* jfreq .* inv_denom_ll.^2;
        end

        % --- Spatial R_a, dR_a (per-cluster, each uses its own c_asa_l) ---
        Ra  = cell(L,1);
        dRa = cell(L,1);
        ls_a = 1;
        if is_LoS && L>=1
            Ra{1}  = Ra_LoS;
            dRa{1} = dRa_LoS;   % LoS: c_asa fixed at 0
            ls_a = 2;
        end
        for ll = ls_a:L
            casa_ll   = casa_vec(ll);
            denoms_ll = 1 + 0.5 * m_vals.^2 * casa_ll^2;   % (M+1) x 1
            if any(~isfinite(denoms_ll))
                fval = inf; g = zeros(2*L,1); return;
            end
            scales_ll  = 1 ./ denoms_ll;
            dscales_ll = -m_vals.^2 * casa_ll ./ denoms_ll.^2;
            w_l = basis_w(ll,:).';
            Ra{ll}  = sum(Jm_pilot .* reshape(w_l .* scales_ll,  1,1,[]), 3);
            dRa{ll} = sum(Jm_pilot .* reshape(w_l .* dscales_ll, 1,1,[]), 3);
        end

        % --- Evaluate the selected covariance representation ---
        g = zeros(2*L, 1);
        if strcmp(parameter_backend, 'dense')
            Cyy = sigma2 * eye(d);
            for ll = 1:L
                Cyy = Cyy + p_est(ll) * kron_fast(Rf{ll}, Ra{ll});
            end
            Cyy = (Cyy + Cyy') / 2;

            if any(~isfinite(Cyy), 'all')
                fval = inf; g = zeros(2*L,1); return;
            end
            [Lc, flag] = chol(Cyy, 'lower');
            if flag ~= 0
                fval = inf; g = zeros(2*L,1); return;
            end
            Cinv = Lc' \ (Lc \ eye(d));

            fval = 2*sum(log(diag(Lc))) + ...
                real(sum(Cinv .* likelihood_data.', 'all'));

            % Q = C^{-1}(C - Rhat)C^{-1}, shared by all gradients.
            Q = Cinv - Cinv * likelihood_data * Cinv;
            for ll = 1:L
                dCyy_ll_cds  = p_est(ll) * kron_fast(dRf{ll}, Ra{ll});
                dCyy_ll_casa = p_est(ll) * kron_fast(Rf{ll},  dRa{ll});
                g(2*ll-1) = cds_vec(ll)  * ...
                    real(sum(Q .* dCyy_ll_cds.',  'all'));
                g(2*ll)   = casa_vec(ll) * ...
                    real(sum(Q .* dCyy_ll_casa.', 'all'));
            end
        else
            spectra = zeros(d, L);
            spectra_cds = zeros(d, L);
            spectra_casa = zeros(d, L);
            for ll = 1:L
                sf = dft_covariance_diagonal(Rf{ll}, true);
                sa = dft_covariance_diagonal(Ra{ll}, true);
                dsf = dft_covariance_diagonal(dRf{ll}, false);
                dsa = dft_covariance_diagonal(dRa{ll}, false);
                spectra(:, ll) = kron(sf, sa);
                spectra_cds(:, ll) = kron(dsf, sa);
                spectra_casa(:, ll) = kron(sf, dsa);
            end
            variance = combine_beam_delay_spectra(spectra, p_est, sigma2);
            spectral_residual = 1 ./ variance - likelihood_data ./ variance.^2;
            fval = real(sum(log(variance) + likelihood_data ./ variance));
            for ll = 1:L
                g(2*ll-1) = cds_vec(ll) * p_est(ll) * ...
                    real(sum(spectral_residual .* spectra_cds(:, ll)));
                g(2*ll) = casa_vec(ll) * p_est(ll) * ...
                    real(sum(spectral_residual .* spectra_casa(:, ll)));
            end
        end
        % Zero out gradient for fixed clusters (LoS specular: l=1)
        if is_LoS && L >= 1
            g(1) = 0;  % c_ds_1 fixed
            g(2) = 0;  % c_asa_1 fixed
        end
    end

    function r = lbfgs_dir(g, Ss, Ys, rhos)
        %LBFGS_DIR Compute the L-BFGS descent direction.
        % Inputs:
        %   g     - Current objective gradient.
        %   Ss    - Stored parameter-displacement vectors.
        %   Ys    - Stored gradient-displacement vectors.
        %   rhos  - Reciprocal curvature products.
        % Output:
        %   r     - Approximate inverse-Hessian descent direction.

        q = g;
        k = size(Ss, 2);
        a_arr = zeros(k, 1);
        for ii = k:-1:1
            a_arr(ii) = rhos(ii) * (Ss(:,ii)' * q);
            q = q - a_arr(ii) * Ys(:,ii);
        end
        if k > 0
            gamma = (Ss(:,k)'*Ys(:,k)) / (Ys(:,k)'*Ys(:,k));
            r = gamma * q;
        else
            r = q;
        end
        for ii = 1:k
            b = rhos(ii) * (Ys(:,ii)' * r);
            r = r + (a_arr(ii) - b) * Ss(:,ii);
        end
        r = -r;
    end


end

function [p_hat, c_ds_vec_hat, c_asa_vec_hat, history] = AlterOptiClusterPowerAndSpreadFactor( ...
    y_p, sigma2, p_init, c_ds_vec_init, c_asa_vec_init, ...
    tau_hat, sigma_tau, sigma_theta, theta_hat, ...
    f_full, pilot_sc, pilot_ant, N, is_LoS, opts)
%ALTEROPTICLUSTERPOWERANDSPREADFACTOR Jointly estimate powers and spreads.
% Alternates a SQUAREM-accelerated EM power update with a joint L-BFGS
% update of all per-cluster delay and angular spread factors.
%
% Inputs:
%   y_p              - Pilot observations of size d-by-S.
%   sigma2           - Noise variance.
%   p_init           - Initial cluster-power vector.
%   c_ds_vec_init    - Initial delay spread factors.
%   c_asa_vec_init   - Initial angular spread factors.
%   tau_hat          - Prior mean cluster delays.
%   sigma_tau        - Standard deviation of the delay prior.
%   sigma_theta      - Standard deviation of the AoA prior.
%   theta_hat        - Prior mean cluster AoAs in radians.
%   f_full           - Full-band subcarrier-frequency vector.
%   pilot_sc         - Pilot subcarrier indices.
%   pilot_ant        - Pilot antenna indices.
%   N                - Total number of receive antennas.
%   is_LoS           - True for the CDL-D LoS special-cluster treatment.
%   opts             - Options struct:
%       .ao_max_iter, .ao_tol, .em_max_iter, .em_tol,
%       .lbfgs_max_iter, .lbfgs_tol, .lbfgs_m, .parameter_backend.
%       .fit_powers, .fit_spreads: optional ablation flags, both true by default.
%
% Outputs:
%   p_hat            - Final cluster-power estimates.
%   c_ds_vec_hat     - Final delay spread estimates.
%   c_asa_vec_hat    - Final angular spread estimates.
%   history          - Struct:
%       .p           - Power estimates after each optimization stage.
%       .c_ds_vec    - Delay spread estimates after each stage.
%       .c_asa_vec   - Angular spread estimates after each stage.
%       .nll         - Negative log-likelihood after each stage.
%       .stage       - Stage labels ('EM' or 'LBFGS').

L  = numel(p_init);
Kp = numel(pilot_sc);
Np = numel(pilot_ant);
d  = Kp * Np;
S  = size(y_p, 2);

if ~isfield(opts, 'parameter_backend'), opts.parameter_backend = 'dense'; end
parameter_backend = opts.parameter_backend;

% Cache the sufficient observation statistic for NLL evaluations.
if strcmp(parameter_backend, 'dense')
    likelihood_data = (y_p * y_p') / S;
elseif strcmp(parameter_backend, 'beam_delay_diag')
    likelihood_data = beam_delay_sample_power(y_p, Np, Kp);
else
    error('Base:UnknownParameterBackend', ...
        'Unsupported parameter_backend: %s', parameter_backend);
end

% Initialize the alternating-optimization state.
if ~isfield(opts, 'fit_powers'), opts.fit_powers = true; end
if ~isfield(opts, 'fit_spreads'), opts.fit_spreads = true; end
p_cur         = p_init(:);
c_ds_vec_cur  = c_ds_vec_init(:);
c_asa_vec_cur = c_asa_vec_init(:);

% Store one history entry after each EM and L-BFGS stage.
max_entries   = 2 * opts.ao_max_iter;
hist_p        = zeros(L, max_entries);
hist_cds_vec  = cell(1, max_entries);
hist_casa_vec = cell(1, max_entries);
hist_nll      = zeros(1, max_entries);
hist_stage    = cell(1, max_entries);
idx = 0;

prev_nll = inf;

start = tic;
for ao_iter = 1:opts.ao_max_iter
    % Update cluster powers with the current spread factors fixed.
    [~, R_factors_pp] = build_cluster_R_factors_hp_yy( ...
        is_LoS, N, tau_hat, sigma_tau, c_ds_vec_cur, ...
        theta_hat, sigma_theta, c_asa_vec_cur, f_full, pilot_sc, pilot_ant);

    em_opts = struct('max_iter', opts.em_max_iter, 'tol', opts.em_tol, ...
        'parameter_backend', parameter_backend);
    if opts.fit_powers
        [p_cur, ~, ~] = ClusterPowersEsti_SQUAREM(y_p, R_factors_pp, sigma2, p_cur, em_opts);
    end

    nll_em = compute_nll(p_cur, R_factors_pp, sigma2, ...
        likelihood_data, d, parameter_backend);

    idx = idx + 1;
    hist_p(:, idx)        = p_cur;
    hist_cds_vec{idx}     = c_ds_vec_cur;
    hist_casa_vec{idx}    = c_asa_vec_cur;
    hist_nll(idx)         = nll_em;
    hist_stage{idx}       = 'EM';
    if ~opts.fit_powers, hist_stage{idx} = 'PowerFixed'; end

    fprintf('AO iter %d - EM: NLL=%.4f, %.2fs \n', ao_iter, nll_em, toc(start));

    % Jointly update all spread factors with the powers fixed.
    lbfgs_opts = struct('max_iter', opts.lbfgs_max_iter, ...
        'tol', opts.lbfgs_tol, ...
        'lbfgs_m', opts.lbfgs_m, ...
        'parameter_backend', parameter_backend);
    if isfield(opts, 'spread_bounds'), lbfgs_opts.spread_bounds = opts.spread_bounds; end
    if opts.fit_spreads
        [c_ds_vec_cur, c_asa_vec_cur, ~, ~] = SpreadFactorEsti( ...
        y_p, p_cur, sigma2, c_ds_vec_cur, c_asa_vec_cur, ...
        tau_hat, sigma_tau, sigma_theta, theta_hat, ...
        f_full(pilot_sc), pilot_ant, N, is_LoS, lbfgs_opts);
    end

    % Rebuild R_factors_pp with updated per-cluster spreads for NLL evaluation
    [~, R_factors_pp_new] = build_cluster_R_factors_hp_yy( ...
        is_LoS, N, tau_hat, sigma_tau, c_ds_vec_cur, ...
        theta_hat, sigma_theta, c_asa_vec_cur, f_full, pilot_sc, pilot_ant);
    nll_lbfgs = compute_nll(p_cur, R_factors_pp_new, sigma2, ...
        likelihood_data, d, parameter_backend);

    idx = idx + 1;
    hist_p(:, idx)        = p_cur;
    hist_cds_vec{idx}     = c_ds_vec_cur;
    hist_casa_vec{idx}    = c_asa_vec_cur;
    hist_nll(idx)         = nll_lbfgs;
    hist_stage{idx}       = 'LBFGS';
    if ~opts.fit_spreads, hist_stage{idx} = 'SpreadFixed'; end

    fprintf('AO iter %d - LBFGS: NLL=%.4f  mean_cds=%.3e  mean_casa=%.3f deg, %.2fs \n', ...
        ao_iter, nll_lbfgs, mean(c_ds_vec_cur), rad2deg(mean(c_asa_vec_cur)), toc(start));

    % Stop when the relative outer-loop objective change is sufficiently small.
    if abs(prev_nll - nll_lbfgs) / (abs(prev_nll) + eps) < opts.ao_tol
        fprintf('AO converged at iter %d (relative NLL change < %.1e)\n', ...
            ao_iter, opts.ao_tol);
        break;
    end
    prev_nll = nll_lbfgs;
end

% Trim history
history.p         = hist_p(:, 1:idx);
history.c_ds_vec  = hist_cds_vec(1:idx);
history.c_asa_vec = hist_casa_vec(1:idx);
history.nll       = hist_nll(1:idx);
history.stage     = hist_stage(1:idx);

p_hat         = p_cur;
c_ds_vec_hat  = c_ds_vec_cur;
c_asa_vec_hat = c_asa_vec_cur;
end

function spectrum = dft_covariance_diagonal(R, enforce_nonnegative)
%DFT_COVARIANCE_DIAGONAL Return diag(F*R*F'), with unitary DFT matrix F.
if nargin < 2, enforce_nonnegative = true; end
R = (R + R') / 2;
R_dft = ifft(fft(R, [], 1), [], 2);
spectrum = real(diag(R_dft));
if enforce_nonnegative
    spectrum = max(spectrum, 0);
end
end

function spectra = build_beam_delay_spectra(R_factors_pp)
%BUILD_BEAM_DELAY_SPECTRA Build per-cluster diagonal beam-delay spectra.
L = numel(R_factors_pp.tau);
Np = size(R_factors_pp.theta{1}, 1);
Kp = size(R_factors_pp.tau{1}, 1);
spectra = zeros(Np * Kp, L);
for l = 1:L
    spectrum_f = dft_covariance_diagonal(R_factors_pp.tau{l}, true);
    spectrum_a = dft_covariance_diagonal(R_factors_pp.theta{l}, true);
    spectra(:, l) = kron(spectrum_f, spectrum_a);
end
end

function variance = combine_beam_delay_spectra(spectra, p, sigma2)
%COMBINE_BEAM_DELAY_SPECTRA Form the positive aggregate spectral variance.
variance = real(sigma2 + spectra * p(:));
variance = max(variance, max(eps, real(sigma2) * 1e-12));
end

function sample_power = beam_delay_sample_power(y_p, Np, Kp)
%BEAM_DELAY_SAMPLE_POWER Return the snapshot-averaged unitary-DFT power.
S = size(y_p, 2);
Y = reshape(y_p, Np, Kp, S);
Y_dft = fft2(Y) / sqrt(Np * Kp);
sample_power = mean(abs(Y_dft).^2, 3);
sample_power = sample_power(:);
end

function x = solve_beam_delay_diagonal(rhs, variance, Np, Kp)
%SOLVE_BEAM_DELAY_DIAGONAL Apply F' * diag(1./variance) * F.
num_rhs = size(rhs, 2);
RHS = reshape(rhs, Np, Kp, num_rhs);
variance_grid = reshape(variance, Np, Kp);
X = ifft2(fft2(RHS) ./ variance_grid);
x = reshape(X, Np * Kp, num_rhs);
end

function embedding = toeplitz_circulant_embedding(R)
%TOEPLITZ_CIRCULANT_EMBEDDING Embed a Toeplitz matrix in a 2n circulant.
R = (R + R') / 2;
embedding = [R(:, 1); 0; R(1, end:-1:2).'];
end

function operator_spectrum = build_bttb_fft_operator(p, sigma2, R_factors_pp)
%BUILD_BTTB_FFT_OPERATOR Build the exact 2-D circulant-embedding spectrum.
L = numel(p);
Np = size(R_factors_pp.theta{1}, 1);
Kp = size(R_factors_pp.tau{1}, 1);
operator_spectrum = zeros(2 * Np, 2 * Kp);
for l = 1:L
    embed_a = toeplitz_circulant_embedding(R_factors_pp.theta{l});
    embed_f = toeplitz_circulant_embedding(R_factors_pp.tau{l});
    operator_spectrum = operator_spectrum + p(l) * ...
        (fft(embed_a) * fft(embed_f).');
end
operator_spectrum = operator_spectrum + sigma2;
end

function y = apply_bttb_fft_operator(x, operator_spectrum, Np, Kp)
%APPLY_BTTB_FFT_OPERATOR Apply the original BTTB matrix by 2-D FFT.
num_rhs = size(x, 2);
X = reshape(x, Np, Kp, num_rhs);
X_pad = zeros(2 * Np, 2 * Kp, num_rhs, 'like', x);
X_pad(1:Np, 1:Kp, :) = X;
Y_pad = ifft2(fft2(X_pad) .* operator_spectrum);
Y = Y_pad(1:Np, 1:Kp, :);
y = reshape(Y, Np * Kp, num_rhs);
end

function [p_hat, c_ds_vec_hat, c_asa_vec_hat] = EstimateProposedParameters( ...
    yp_used, sigma2, sys, chann, estor, tau_hat, theta_hat, ...
    pilot_sc, pilot_ant, parameter_backend)
%ESTIMATEPROPOSEDPARAMETERS Run the shared initialization and selected AO backend.
L = length(chann.tau);
[~, R_factors_pp_init] = build_cluster_R_factors_hp_yy( ...
    chann.is_LoS, sys.N, tau_hat, estor.sigma_tau_hat, ...
    chann.c_ds_ref * ones(L, 1), theta_hat, estor.sigma_theta_hat, ...
    chann.c_asa_ref * ones(L, 1), sys.f, pilot_sc, pilot_ant);
p_init = ClusterPowersEsti_MomentMatching(yp_used, R_factors_pp_init, sigma2);

ao_opts = estor.ao_opts;
ao_opts.parameter_backend = parameter_backend;
[p_hat, c_ds_vec_hat, c_asa_vec_hat, ~] = ...
    AlterOptiClusterPowerAndSpreadFactor( ...
    yp_used, sigma2, p_init, chann.c_ds_ref * ones(L, 1), ...
    chann.c_asa_ref * ones(L, 1), tau_hat, estor.sigma_tau_hat, ...
    estor.sigma_theta_hat, theta_hat, sys.f, pilot_sc, pilot_ant, ...
    sys.N, chann.is_LoS, ao_opts);
end

function [h, y, tau_hat, theta_hat] = DrawTrial(sys, chann, funcs)
%DRAWTRIAL Original channel/noise/prior draw order, shared by all estimators.
sc = (1:sys.K/sys.Kp:sys.K).';
ant = (1:sys.Np).';
L = numel(chann.tau);
p = chann.cluster_powers_NoBlockage;
blocked = rand(L,1) < chann.ClusterBlockProb;
p(blocked) = 0.001*p(blocked);
[h,y] = funcs.GenerateChannelObservation(chann.theta,chann.tau, ...
    chann.c_asa_vec,chann.c_ds_vec,chann.M_subpath,p, ...
    10^(-chann.SNR_dB/10),sys.f,sys.N,sc,ant, ...
    chann.SpreadTruncFactor,chann.is_LoS);
theta_hat = chann.theta + chann.sigma_theta*randn(L,1);
if chann.is_LoS
    e = chann.sigma_tau*randn(L,1); e(1:2) = 0;
else
    e = [0;chann.sigma_tau*randn(L-1,1)];
end
tau_hat = chann.tau + e;
end

function h = PilotOnlyEstimate(y,K,N,sc,method)
%PILOTONLYESTIMATE Pilot LS plus frequency interpolation, without CCM/DS/CP.
% Uniform comb pilots start at subcarrier 1; every antenna is observed.
% Real and imaginary components are interpolated separately for PCHIP.
sc = sc(:); Kp = numel(sc); D = K/Kp;
Y = reshape(y,N,Kp);
switch lower(method)
    case {'linear','pchip','spline'}
        x = [sc-K;sc;sc+K]; V = [Y,Y,Y].';
        H = (interp1(x,real(V),(1:K).',method) ...
            + 1j*interp1(x,imag(V),(1:K).',method)).';
    case 'quadratic'
        k = 0:K-1; m = floor(k/D); a = k/D-m;
        H = Y(:,mod(m-1,Kp)+1).*(a.*(a-1)/2) ...
            + Y(:,m+1).*(1-a.^2) ...
            + Y(:,mod(m+1,Kp)+1).*(a.*(a+1)/2);
    case 'fir'
        % Windowed-sinc low-pass: 8D+1 taps, cutoff pi/D; zero-phase
        % circular convolution. Normalize each polyphase DC gain.
        n = -4*D:4*D; t = n/D;
        s = ones(size(t)); nz = t~=0; s(nz) = sin(pi*t(nz))./(pi*t(nz));
        taps = s.*(0.54+0.46*cos(pi*n/(4*D)));
        b = accumarray(mod(n(:),K)+1,taps(:),[K,1]).';
        z = zeros(N,K); z(:,sc) = Y;
        mask = zeros(1,K); mask(sc) = 1;
        gain = real(ifft(fft(mask).*fft(b)));
        H = ifft(fft(z,[],2).*fft(b),[],2)./gain;
    case 'dft'
        % Causal IDFT -> append zeros -> DFT. No privileged delay truncation.
        H = fft(ifft(Y,[],2),K,2);
    otherwise
        error('Unknown pilot-only method: %s',method);
end
h = H(:);
end

function [al,f_al,g_al,p,status,reset] = SafeLineSearch(fun,x,p,f0,g0,opts)
%SAFELINESEARCH Strong Wolfe, verified Armijo fallback, or an unchanged point.
% The returned direction can differ from p: non-descent/invalid directions
% are reset to -g0. All returned values refer to exactly x + al*p.
if nargin<6, opts=struct(); end
if ~isfield(opts,'max_ls'), opts.max_ls=20; end
if ~isfield(opts,'max_zoom'), opts.max_zoom=10; end
if ~isfield(opts,'max_backtrack'), opts.max_backtrack=40; end
assert(isfinite(f0) && all(isfinite(x)) && all(isfinite(g0)), ...
    'SafeLineSearch:InvalidStart','Initial point, value and gradient must be finite.');
c1=1e-4; c2=.9;
reset=any(~isfinite(p)) || ~isfinite(g0'*p) || g0'*p>=0;
if reset, p=-g0; end
al=0; f_al=f0; g_al=g0; status='no_step';
if all(g0==0), status='stationary'; return; end
dg0=g0'*p;
if ~isfinite(dg0) || dg0>=0, return; end
best_a=0; best_f=f0; best_g=g0;
a_prev=0; f_prev=f0; dg_prev=dg0; a_cur=1;
for k=1:opts.max_ls
    [f,g,valid]=probe(a_cur); dg=g'*p;
    if ~valid || f>f0+c1*a_cur*dg0 || (k>1 && f>=f_prev)
        if zoom(a_prev,a_cur,f_prev,f,dg_prev), return; end
        break;
    end
    if abs(dg)<=c2*abs(dg0) && f<f0
        accept(a_cur,f,g,'strong_wolfe'); return;
    end
    if dg>=0
        if zoom(a_cur,a_prev,f,f_prev,dg), return; end
        break;
    end
    a_prev=a_cur; f_prev=f; dg_prev=dg;
    if a_cur>=100, break; end
    a_cur=min(2*a_cur,100);
end
% Never return an unevaluated step or an unchecked last zoom iterate.
if best_a>0
    accept(best_a,best_f,best_g,'armijo'); return;
end
a=1;
for k=1:opts.max_backtrack
    probe(a);
    if best_a>0
        accept(best_a,best_f,best_g,'armijo'); return;
    end
    a=a/2;
end
% Failure deliberately returns zero with the original value and gradient.

    function [f,g,valid]=probe(a)
        trial=x+a*p;
        f=inf; g=zeros(size(g0));
        if all(isfinite(trial)), [f,g]=fun(trial); end
        valid=isfinite(f) && all(isfinite(g));
        if ~valid, f=inf; g=zeros(size(g0)); end
        if valid && f<f0 && f<=f0+c1*a*dg0 && f<best_f
            best_a=a; best_f=f; best_g=g;
        end
    end

    function accept(a,f,g,label)
        al=a; f_al=f; g_al=g; status=label;
    end

    function found=zoom(lo,hi,flo,fhi,dglo)
        found=false;
        for j=1:opts.max_zoom
            da=hi-lo; denom=fhi-flo-dglo*da;
            if isfinite(denom) && abs(da)>1e-20 && denom>1e-30*abs(da^2)
                a=lo-.5*dglo*da^2/denom;
                a=max(min(lo,hi)+.1*abs(da),min(a,max(lo,hi)-.1*abs(da)));
            else
                a=(lo+hi)/2;
            end
            [f,g,valid]=probe(a); dg=g'*p;
            if ~valid || f>f0+c1*a*dg0 || f>=flo
                hi=a; fhi=f;
            else
                if abs(dg)<=c2*abs(dg0) && f<f0
                    accept(a,f,g,'strong_wolfe'); found=true; return;
                end
                if dg*(hi-lo)>=0, hi=lo; fhi=flo; end
                lo=a; flo=f; dglo=dg;
            end
        end
    end
end
