function exp_ParameterAccuracy(mc_num,folder,finite_subpaths)
%EXP_PARAMETERACCURACY Truncated finite-subpath fitting, independent validation.
% Third argument false reproduces the earlier Gaussian covariance control.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','parameter_accuracy_truncated'); end
if nargin<3, finite_subpaths=true; end
[s,c,~,e,f]=Base(false); g=GaussianCovarianceTools; grid=[32,64]; snapshots=[1,8];
L=numel(c.tau); ant=(1:s.N).'; sigma2=10^(-c.SNR_dB/10); ptrue=c.cluster_powers_NoBlockage;
% Fixed prior means; Gaussian prior uncertainty is already marginalized in R.
tau=c.tau; theta=c.theta;
rt=0; dtarget=c.c_ds_vec; atarget=c.c_asa_vec;
if finite_subpaths
    rt=c.SpreadTruncFactor; dtarget=c.c_ds_eq_vec; atarget=c.c_asa_eq_vec;
end
[~,Rtrue]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_vec,theta,c.sigma_theta,c.c_asa_vec,s.f,(1:s.K).',ant,rt);
if ~finite_subpaths, Bfull=g.Roots(Rtrue); end
result=struct('sys',s,'channel',c,'estor',e,'mc_num',mc_num,'grid',grid,'snapshots',snapshots, ...
    'seeds',280925+(1:mc_num),'model','Zero-mean Gaussian; untruncated full separable covariance; fixed prior means', ...
    'methods',{{'Initialized','Proposed','Oracle'}}, ...
    'metrics',{{'Relative power error','Delay log-spread RMSE','Angle log-spread RMSE','Relative covariance error'}}, ...
    'error',zeros(mc_num,2,2,3),'energy',zeros(mc_num,2), ...
    'metric',zeros(mc_num,2,2,2,4),'power',zeros(L,mc_num,2,2,2), ...
    'delay_spread',zeros(L,mc_num,2,2,2),'angle_spread',zeros(L,mc_num,2,2,2), ...
    'nmse_db',zeros(2,2,3),'ci',zeros(2,2,2,3),'values',[], ...
    'columns',{{'Kp','snapshots','Initial_NMSE','Proposed_NMSE','Oracle_NMSE', ...
    'Initial_power_error','Fitted_power_error','Initial_delay_log_RMSE','Fitted_delay_log_RMSE', ...
    'Initial_angle_log_RMSE','Fitted_angle_log_RMSE','Initial_covariance_error','Fitted_covariance_error'}});
result.source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']), ...
    'gaussian',fileread(which('GaussianCovarianceTools')));
result.finite_subpaths=finite_subpaths;
result.target_delay_rms=dtarget; result.target_angle_rms=atarget;
if finite_subpaths
    result.model='Manuscript finite weighted subpaths; truncated profiles; independent center uncertainty per snapshot';
    result.metrics{2}='Delay RMS-spread log error';
    result.metrics{3}='Angle RMS-spread log error';
end
for j=1:2
    s.Kp=grid(j); sc=(1:s.K/s.Kp:s.K).'; idx=f.PilotObservationIndices(sc,ant,s.N);
    [Rhtrue,Rpp]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_vec,theta,c.sigma_theta,c.c_asa_vec,s.f,sc,ant,rt);
    if ~finite_subpaths, Bpilot=g.Roots(Rpp); end
    d0=c.c_ds_ref*ones(L,1); a0=c.c_asa_ref*ones(L,1);
    [Rh0,Rp0]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,d0,theta,e.sigma_theta_hat,a0,s.f,sc,ant);
    [~,Rfull0]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,d0,theta,e.sigma_theta_hat,a0,s.f,(1:s.K).',ant);
    for mc=1:mc_num
        rng(result.seeds(mc),'twister');
        if finite_subpaths
            [h,ytest]=finite_snapshot(s,c,f,tau,theta,sc,ant,sigma2);
            ytrain=zeros(numel(idx),8);
            for snap=1:8
                [~,ytrain(:,snap)]=finite_snapshot(s,c,f,tau,theta,sc,ant,sigma2);
            end
        else
            h=g.Sample(Bfull,ptrue,1);
            ytest=h(idx)+sqrt(sigma2/2)*(randn(numel(idx),1)+1j*randn(numel(idx),1));
            ytrain=g.Sample(Bpilot,ptrue,8)+sqrt(sigma2/2)*(randn(numel(idx),8)+1j*randn(numel(idx),8));
        end
        result.energy(mc,j)=sum(abs(h).^2);
        horacle=f.EstimateChannelFromFactors(ytest,ptrue,sigma2,Rpp,Rhtrue,e.final_solver,e.pcg_opts);
        for si=1:2
            y=ytrain(:,1:snapshots(si));
            p0=f.MomentMatching(y,Rp0,sigma2);
            opts=e.ao_opts; opts.parameter_backend=e.parameter_backend;
            hist=struct;
            evalc('[ph,dh,ah,hist]=f.TraceAO(y,sigma2,p0,d0,a0,tau,e.sigma_tau_hat,e.sigma_theta_hat,theta,s.f,sc,ant,s.N,false,opts);');
            result.ao_history{mc,j,si}=hist;
            for stage=1:2
                if stage==1
                    p=p0; d=d0; a=a0; Rh=Rh0; R=Rp0; Rfull=Rfull0;
                else
                    p=ph; d=dh; a=ah;
                    [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,d,theta,e.sigma_theta_hat,a,s.f,sc,ant);
                    [~,Rfull]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,d,theta,e.sigma_theta_hat,a,s.f,(1:s.K).',ant);
                end
                he=f.EstimateChannelFromFactors(ytest,p,sigma2,R,Rh,e.final_solver,e.pcg_opts);
                result.error(mc,j,si,stage)=sum(abs(he-h).^2);
                result.metric(mc,j,si,stage,:)=[norm(p-ptrue)/norm(ptrue), ...
                    sqrt(mean(log(d./dtarget).^2)),sqrt(mean(log(a./atarget).^2)), ...
                    g.RelativeError(p,Rfull,ptrue,Rtrue)];
                result.power(:,mc,j,si,stage)=p;
                result.delay_spread(:,mc,j,si,stage)=d;
                result.angle_spread(:,mc,j,si,stage)=a;
            end
            result.error(mc,j,si,3)=sum(abs(horacle-h).^2);
        end
        if mod(mc,10)==0 || mc==mc_num, fprintf('Accuracy Kp=%d: %d/%d\n',s.Kp,mc,mc_num); end
    end
    for si=1:2
        err=reshape(result.error(:,j,si,:),mc_num,3);
        [db,ci,boot]=ResponseStatistics(err,result.energy(:,j),true);
        result.nmse_db(j,si,:)=db; result.ci(j,si,:,:)=ci; result.bootstrap{j,si}=boot;
        metrics=reshape(result.metric(:,j,si,:,:),mc_num,8);
        [avg,mci,mb]=ResponseStatistics(metrics);
        result.metric_mean(j,si,:,:)=reshape(avg,2,4);
        result.metric_ci(j,si,:,:,:)=reshape(mci,2,2,4);
        result.metric_bootstrap{j,si}=mb;
        result.values(end+1,:)=[s.Kp,snapshots(si),db,avg];
    end
    SaveResult(folder,'parameter_accuracy',result);
end
end

function [h,y]=finite_snapshot(s,c,f,tau,theta,sc,ant,sigma2)
% Same finite-subpath generator as Base, with the historical-mode center law.
ts=tau+[0;c.sigma_tau*randn(numel(tau)-1,1)];
as=theta+c.sigma_theta*randn(numel(theta),1);
[h,y]=f.GenerateChannelObservation(as,ts,c.c_asa_vec,c.c_ds_vec, ...
    c.M_subpath,c.cluster_powers_NoBlockage,sigma2,s.f,s.N,sc,ant, ...
    c.SpreadTruncFactor,false);
end
