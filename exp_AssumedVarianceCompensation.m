function exp_AssumedVarianceCompensation(mc_num,folder,workers,snr_db,kp)
%EXP_ASSUMEDVARIANCECOMPENSATION Fixed data and powers; scan assumed prior SD.
if nargin<1, mc_num=20; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','assumed_variance_compensation'); end
if nargin<3, workers=4; end
if nargin<4, snr_db=15; end
if nargin<5, kp=64; end
[s,c,~,e,f]=Base(false); s.Kp=kp; c.SNR_dB=snr_db;
grid=[.25,.5,1,2,4]; L=numel(c.tau); J=numel(grid);
e.ao_opts.spread_bounds=[c.c_ds_ref;c.c_asa_ref]*[.05,4];
result=struct('sys',s,'channel',c,'estor',e,'grid',grid,'seed',260927, ...
    'mc_num',mc_num,'workers',workers,'methods',{{'Spread refit (fixed power)','Oracle','Frozen at a = 1'}}, ...
    'error',zeros(mc_num,2,J,3),'energy',zeros(mc_num,1), ...
    'power',zeros(L,mc_num),'delay_spread',zeros(L,mc_num,2,J), ...
    'angle_spread',zeros(L,mc_num,2,J),'anchor_delay',zeros(L,mc_num), ...
    'anchor_angle',zeros(L,mc_num),'nmse_db',zeros(2,J,3),'ci',zeros(2,J,2,3), ...
    'spread_mean',zeros(2,J,2),'spread_ci',zeros(2,J,2,2), ...
    'values',[],'columns',{{'axis','a_actual_over_assumed_SD','Spread_refit','Oracle','Frozen', ...
    'mean_delay_ns','mean_angle_deg','gain_db','gain_ci_low','gain_ci_high'}});
result.source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']));
rng(result.seed,'twister'); data=cell(mc_num,1); trials=cell(mc_num,1);
for mc=1:mc_num
    [data{mc}.h,data{mc}.y,data{mc}.tau,data{mc}.theta]=f.DrawTrial(s,c,f);
end
if workers>0, pool=parpool('Processes',workers); end
parfor (mc=1:mc_num,workers)
    trials{mc}=fit_trial(data{mc},s,c,e,f,grid);
    fprintf('Assumed-variance trial %d finished\n',mc);
end
if workers>0, delete(pool); end
for mc=1:mc_num
    t=trials{mc}; result.energy(mc)=t.energy; result.error(mc,:,:,:)=t.error;
    result.power(:,mc)=t.power; result.delay_spread(:,mc,:,:)=t.delay;
    result.angle_spread(:,mc,:,:)=t.angle; result.anchor_delay(:,mc)=t.anchor_delay;
    result.anchor_angle(:,mc)=t.anchor_angle; result.histories{mc}=t.histories;
    result.calibration{mc}=t.calibration; result.console{mc}=t.console;
end
assert(isequal(result.error(:,:,3,1),result.error(:,:,3,3)));
assert(isequal(result.error(:,1,3,:),result.error(:,2,3,:)));
for axis=1:2
    for j=1:J
        assert(isequal(result.error(:,axis,j,2),result.error(:,1,1,2)));
        err=reshape(result.error(:,axis,j,:),mc_num,3);
        [db,ci,boot]=ResponseStatistics(err,result.energy,true);
        result.nmse_db(axis,j,:)=db; result.ci(axis,j,:,:)=ci;
        gain=db(3)-db(1); bg=sort(boot(:,3)-boot(:,1)); gain_ci=bg([50,1950]);
        result.gain_db(axis,j)=gain; result.gain_ci(axis,j,:)=gain_ci;
        result.bootstrap{axis,j}=boot;
        ds=result.delay_spread(:,:,axis,j)*1e9; as=result.angle_spread(:,:,axis,j)*180/pi;
        [mu,ci]=ResponseStatistics([mean(ds,1).',mean(as,1).'],ones(mc_num,1),false);
        result.spread_mean(axis,j,:)=mu; result.spread_ci(axis,j,:,:)=ci;
        result.values(end+1,:)=[axis,grid(j),db,mu,gain,gain_ci.'];
    end
end
SaveResult(folder,'assumed_variance_compensation',result);
plot_AssumedVarianceCompensation(folder);
disp(result.values);
end

function t=fit_trial(data,s,c,e,f,grid)
maxNumCompThreads(2);
h=data.h; y=data.y; tau=data.tau; theta=data.theta;
L=numel(c.tau); J=numel(grid); sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).';
sigma2=10^(-c.SNR_dB/10); dinit=c.c_ds_ref*ones(L,1); ainit=c.c_asa_ref*ones(L,1);
opts=e.ao_opts; opts.parameter_backend=e.parameter_backend;
[~,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,dinit,theta,e.sigma_theta_hat,ainit,s.f,sc,ant);
pinit=f.MomentMatching(y,R,sigma2);
t=struct('energy',sum(abs(h).^2),'error',zeros(2,J,3), ...
    'delay',zeros(L,2,J),'angle',zeros(L,2,J),'histories',{cell(2,J)},'console',{cell(2,J)});
% Matched-uncertainty calibration sets the shared power. It is never refitted below.
evalc('[p,~,~,t.calibration]=f.TraceAO(y,sigma2,pinit,dinit,ainit,tau,e.sigma_tau_hat,e.sigma_theta_hat,theta,s.f,sc,ant,s.N,false,opts);');
t.power=p; opts.fit_powers=false;
% Same reference-spread initialization and iteration budget at every a.
for axis=1:2
    for j=1:J
        st=c.sigma_tau; sa=c.sigma_theta;
        if axis==1, st=st/grid(j); else, sa=sa/grid(j); end
        hist=struct; pf=p; d=dinit; a=ainit;
        t.console{axis,j}=evalc('[pf,d,a,hist]=f.TraceAO(y,sigma2,p,dinit,ainit,tau,st,sa,theta,s.f,sc,ant,s.N,false,opts);');
        assert(isequal(pf,p));
        assert(all(d>=opts.spread_bounds(1,1)*(1-1e-12)&d<=opts.spread_bounds(1,2)*(1+1e-12)));
        assert(all(a>=opts.spread_bounds(2,1)*(1-1e-12)&a<=opts.spread_bounds(2,2)*(1+1e-12)));
        t.delay(:,axis,j)=d; t.angle(:,axis,j)=a; t.histories{axis,j}=hist;
    end
end
t.anchor_delay=t.delay(:,1,3); t.anchor_angle=t.angle(:,1,3);
[Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_vec, ...
    theta,c.sigma_theta,c.c_asa_vec,s.f,sc,ant,c.SpreadTruncFactor);
he=f.EstimateChannelFromFactors(y,c.cluster_powers_NoBlockage,sigma2,R,Rh,e.final_solver,e.pcg_opts);
oracle_error=sum(abs(he-h).^2);
for axis=1:2
    for j=1:J
        st=c.sigma_tau; sa=c.sigma_theta;
        if axis==1, st=st/grid(j); else, sa=sa/grid(j); end
        for m=[1,3]
            d=t.delay(:,axis,j); a=t.angle(:,axis,j);
            if m==3, d=t.anchor_delay; a=t.anchor_angle; end
            [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,st,d,theta,sa,a,s.f,sc,ant);
            he=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,e.final_solver,e.pcg_opts);
            t.error(axis,j,m)=sum(abs(he-h).^2);
        end
        t.error(axis,j,2)=oracle_error;
    end
end
end
