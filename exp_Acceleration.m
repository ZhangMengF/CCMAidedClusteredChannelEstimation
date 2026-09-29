function exp_Acceleration(mc_num,folder,resume,grid,seed_base)
%EXP_ACCELERATION Paired parameter-only timing; exact final-solver check.
% Run alone, with eight MATLAB computational threads. Timing is hardware-
% dependent; seeds, outputs, per-trial times, and thread count are saved.
if nargin<1, mc_num=[20,20,20,5,5]; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results'); end
if nargin<3, resume=false; end
if nargin<4, grid=[8,16,32,64,128]; end
if nargin<5
    seed_base=24090000;
    [valid,seed_points]=ismember(grid,[8,16,32,64,128]);
    assert(all(valid),'Unsupported pilot count.');
    seed_rule='24090000 + 1000*original_grid_point + trial';
else
    seed_points=grid;
    seed_rule=sprintf('%d + 1000*Kp + trial',seed_base);
end
if isscalar(mc_num), mc_num=repmat(mc_num,1,numel(grid)); end
assert(numel(mc_num)==numel(grid) && all(mc_num>0) && all(mc_num==floor(mc_num)), ...
    'Provide a positive integer trial count per pilot count.');
assert(all(ismember(grid,[8,16,32,64,128])),'Unsupported pilot count.');
ng=numel(grid);
max_mc=max(mc_num);
old_threads=maxNumCompThreads(8);
restore=onCleanup(@() maxNumCompThreads(old_threads));
[s,c,~,e,f]=Base(false);
backends={'dense','beam_delay_diag'}; ant=(1:s.N).'; sigma2=10^(-c.SNR_dB/10);
result=struct('sys',s,'channel',c,'estor',e,'grid',grid,'mc_num',mc_num, ...
    'seed_rule',seed_rule,'seed_base',seed_base,'seed_points',seed_points,'backends',{backends}, ...
    'error',zeros(max_mc,2,ng),'energy',zeros(max_mc,ng), ...
    'parameter_seconds',zeros(max_mc,2,ng), ...
    'solver_seconds',zeros(max_mc,2,2,ng), ...
    'solver_relative_error',zeros(max_mc,2,ng), ...
    'dense_final_error',zeros(max_mc,2,ng),'completed_trials',zeros(1,ng));
result.values=nan(ng,7);
result.columns={'Kp','nmse_dense_dB','nmse_accelerated_dB','difference_dB', ...
    'parameter_speedup','max_solver_relative_error','max_solver_nmse_difference_dB'};
if resume
    saved=load(fullfile(folder,'acceleration.mat'),'result'); result=saved.result;
    assert(isequal(result.grid,grid),'Checkpoint pilot grid mismatch.');
    if isfield(result,'seed_base')
        assert(result.seed_base==seed_base && isequal(result.seed_points,seed_points), ...
            'Checkpoint seed mismatch.');
    else
        assert(nargin<5,'Legacy checkpoint requires the original seed rule.');
    end
    result.mc_num=mc_num;
end
for j=1:ng
    s.Kp=grid(j); sc=(1:s.K/s.Kp:s.K).';
    % Warm both paths on the first paired realization, outside the timers.
    if result.completed_trials(j)<mc_num(j)
        fprintf('Warming both backends: Kp=%d (excluded from timing).\n',s.Kp);
        rng(seed_base+1000*seed_points(j)+1,'twister'); [~,y,tau,theta]=f.DrawTrial(s,c,f);
        for b=1:2
            evalc('[p,ds,as]=f.EstimateProposedParameters(y,sigma2,s,c,e,tau,theta,sc,ant,backends{b});');
            [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,ds,theta,c.sigma_theta,as,s.f,sc,ant);
            f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'fft_pcg',e.pcg_opts);
            f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'dense_chol',e.pcg_opts);
        end
    end
    for mc=result.completed_trials(j)+1:mc_num(j)
        rng(seed_base+1000*seed_points(j)+mc,'twister'); [h,y,tau,theta]=f.DrawTrial(s,c,f);
        result.energy(mc,j)=sum(abs(h).^2);
        order=1:2; if mod(mc,2)==0, order=2:-1:1; end
        for b=order
            timer=tic;
            evalc('[p,ds,as]=f.EstimateProposedParameters(y,sigma2,s,c,e,tau,theta,sc,ant,backends{b});');
            result.parameter_seconds(mc,b,j)=toc(timer);
            [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,ds,theta,c.sigma_theta,as,s.f,sc,ant);
            timer=tic; hp=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'fft_pcg',e.pcg_opts);
            result.solver_seconds(mc,1,b,j)=toc(timer);
            timer=tic; hd=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'dense_chol',e.pcg_opts);
            result.solver_seconds(mc,2,b,j)=toc(timer);
            rel=norm(hp-hd)/norm(hd);
            assert(rel<1e-6,'FFT-PCG/dense estimate disagreement');
            result.solver_relative_error(mc,b,j)=rel;
            result.error(mc,b,j)=sum(abs(hp-h).^2);
            result.dense_final_error(mc,b,j)=sum(abs(hd-h).^2);
        end
        fprintf('Acceleration Kp=%d trial=%d/%d, parameter times=%s s\n', ...
            s.Kp,mc,mc_num(j),mat2str(result.parameter_seconds(mc,:,j),4));
        result.completed_trials(j)=mc;
        SaveResult(folder,'acceleration',result);
    end
    use=1:mc_num(j);
    db=10*log10(sum(result.error(use,:,j),1)/sum(result.energy(use,j)));
    dense_db=10*log10(sum(result.dense_final_error(use,:,j),1)/sum(result.energy(use,j)));
    dt=mean(result.parameter_seconds(use,:,j),1);
    result.values(j,:)=[s.Kp,db,db(2)-db(1),dt(1)/dt(2), ...
        max(result.solver_relative_error(use,:,j),[],'all'),max(abs(db-dense_db))];
    assert(max(abs(db-dense_db))<.001);
    SaveResult(folder,'acceleration',result);
end
end
