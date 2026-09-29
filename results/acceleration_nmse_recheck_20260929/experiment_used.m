function exp_AccelerationNMSERecheck(mc_num,folder,workers,grid)
%EXP_ACCELERATIONNMSERECHECK Parallel paired NMSE replication, not a timer.
% Same estimators/settings as exp_Acceleration; independent fixed seed bank.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','acceleration_nmse_recheck_20260929'); end
if nargin<3, workers=8; end
if nargin<4, grid=[16,64,8,32,128]; end
if ~isfolder(folder), mkdir(folder); end
pool=gcp('nocreate');
if isempty(pool), pool=parpool('Processes',workers); end
assert(pool.NumWorkers==workers,'Requested worker count differs from existing pool.');
[s,c,~,e,f]=Base(false);
result=struct('sys',s,'channel',c,'estor',e,'grid',grid,'mc_num',mc_num, ...
    'seed_base',26092900,'seed_rule','26092900 + 1000*Kp + trial', ...
    'workers',workers,'worker_threads',4,'timing_valid',false, ...
    'backends',{{'dense','beam_delay_diag'}}, ...
    'error',nan(mc_num,2,numel(grid)),'energy',nan(mc_num,numel(grid)), ...
    'completed_trials',zeros(1,numel(grid)));
result.columns={'Kp','nmse_dense_dB','nmse_accelerated_dB','difference_dB', ...
    'bootstrap_low_dB','bootstrap_high_dB','accelerated_better_fraction'};
result.values=nan(numel(grid),7);
% Preserve the exact implementation used by this run, including safety fixes.
code_folder=fileparts(mfilename('fullpath'));
if ~isfile(fullfile(folder,'Base_used.m'))
    copyfile(fullfile(code_folder,'Base.m'),fullfile(folder,'Base_used.m'));
    copyfile([mfilename('fullpath'),'.m'],fullfile(folder,'experiment_used.m'));
else
    assert(strcmp(fileread(fullfile(folder,'Base_used.m')),fileread(fullfile(code_folder,'Base.m'))), ...
        'Base changed: use a fresh output folder.');
end
for j=1:numel(grid)
    kp=grid(j);
    for first=1:workers:mc_num
        ids=first:min(first+workers-1,mc_num);
        errors=zeros(numel(ids),2); energies=zeros(numel(ids),1);
        parfor k=1:numel(ids)
            old_threads=maxNumCompThreads(4);
            cleanup=onCleanup(@() maxNumCompThreads(old_threads));
            mc=ids(k); seed=result.seed_base+1000*kp+mc;
            file=fullfile(folder,sprintf('trial_kp%03d_mc%04d.mat',kp,mc));
            if isfile(file)
                saved=load(file,'trial'); trial=saved.trial;
                assert(trial.seed==seed,'Checkpoint seed mismatch.');
            else
                trial=paired_trial(s,c,e,f,kp,seed);
                save_trial(file,trial);
            end
            errors(k,:)=trial.error; energies(k)=trial.energy;
        end
        result.error(ids,:,j)=errors; result.energy(ids,j)=energies;
        result.completed_trials(j)=ids(end);
        use=1:ids(end); err=result.error(use,:,j);
        db=10*log10(sum(err,1)/sum(result.energy(use,j)));
        % Paired percentile bootstrap of the pooled-energy NMSE difference.
        rng(91000+kp,'twister'); bootstrap=zeros(2000,1);
        for b=1:numel(bootstrap)
            pick=randi(numel(use),numel(use),1);
            total=sum(err(pick,:),1);
            bootstrap(b)=10*log10(total(2)/total(1));
        end
        bootstrap=sort(bootstrap);
        result.values(j,:)=[kp,db,db(2)-db(1),bootstrap(50),bootstrap(1950), ...
            mean(err(:,2)<err(:,1))];
        SaveResult(folder,'acceleration_nmse_recheck',result);
        fprintf('Kp=%d paired trials=%d/%d: difference=%.4f dB, 95%% bootstrap [%.4f, %.4f]\n', ...
            kp,ids(end),mc_num,result.values(j,4:6));
    end
end
end

function trial=paired_trial(s,c,e,f,kp,seed)
s.Kp=kp; sc=(1:s.K/kp:s.K).'; ant=(1:s.N).';
sigma2=10^(-c.SNR_dB/10); backends={'dense','beam_delay_diag'};
rng(seed,'twister'); [h,y,tau,theta]=f.DrawTrial(s,c,f);
trial=struct('seed',seed,'error',zeros(1,2),'energy',sum(abs(h).^2));
for b=1:2
    evalc('[p,ds,as]=f.EstimateProposedParameters(y,sigma2,s,c,e,tau,theta,sc,ant,backends{b});');
    [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,ds,theta,c.sigma_theta,as,s.f,sc,ant);
    hp=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'fft_pcg',e.pcg_opts);
    trial.error(b)=sum(abs(hp-h).^2);
end
end

function save_trial(file,trial)
save(file,'trial','-v7');
end
