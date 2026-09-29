function exp_CCMVarianceCompensation(mc_num,folder,workers)
%EXP_CCMVARIANCECOMPENSATION Refit versus per-trial profiles frozen at a=1.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','ccm_variance_compensation'); end
if nargin<3, workers=4; end % Use zero for serial execution.
[s,c,~,e,f]=Base(false); grid=[.5,1,2,4]; L=numel(c.tau);
assert(c.ClusterBlockProb==0); % Oracle powers below correspond to no blockage.
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).'; sigma2=10^(-c.SNR_dB/10);
result=struct('sys',s,'channel',c,'estor',e,'mc_num',mc_num,'grid',grid, ...
    'methods',{{'Proposed','Oracle','Frozen at a = 1'}},'seed',260925, ...
    'error',zeros(mc_num,2,4,3),'energy',zeros(mc_num,1), ...
    'power',zeros(L,mc_num,2,4),'delay_spread',zeros(L,mc_num,2,4), ...
    'angle_spread',zeros(L,mc_num,2,4),'anchor_power',zeros(L,mc_num), ...
    'anchor_delay_spread',zeros(L,mc_num),'anchor_angle_spread',zeros(L,mc_num), ...
    'anchor_tau',zeros(L,mc_num),'anchor_theta',zeros(L,mc_num), ...
    'nmse_db',zeros(2,4,3),'ci',zeros(2,4,2,3),'values',[], ...
    'columns',{{'axis','scale','Proposed','Oracle','Frozen', ...
    'mean_delay_ns','mean_angle_deg','weighted_delay_ns','weighted_angle_deg'}});
result.source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']));
rng(result.seed,'twister');
data=cell(mc_num,1); trials=cell(mc_num,1);
for mc=1:mc_num
    [data{mc}.h,data{mc}.y,data{mc}.tau,data{mc}.theta]=f.DrawTrial(s,c,f);
end
if workers>0, pool=parpool('Processes',workers); end
parfor (mc=1:mc_num,workers)
    trials{mc}=fit_trial(data{mc},s,c,e,f,sc,ant,sigma2,grid);
    fprintf('Compensation trial %d finished\n',mc);
end
if workers>0, delete(pool); end
for mc=1:mc_num
    t=trials{mc};
    result.energy(mc)=t.energy; result.error(mc,:,:,:)=t.error;
    result.anchor_tau(:,mc)=t.anchor_tau; result.anchor_theta(:,mc)=t.anchor_theta;
    result.anchor_power(:,mc)=t.anchor_power;
    result.anchor_delay_spread(:,mc)=t.anchor_delay_spread;
    result.anchor_angle_spread(:,mc)=t.anchor_angle_spread;
    result.power(:,mc,:,:)=t.power; result.delay_spread(:,mc,:,:)=t.delay_spread;
    result.angle_spread(:,mc,:,:)=t.angle_spread;
end
result.workers=workers;
assert(isequal(result.error(:,:,2,1),result.error(:,:,2,3)));
assert(isequal(result.error(:,1,2,:),result.error(:,2,2,:)));
result.spread_mean=zeros(2,4,2); result.spread_ci=zeros(2,4,2,2);
result.spread_weighted_mean=zeros(2,4,2);
for axis=1:2
    for j=1:4
        err=reshape(result.error(:,axis,j,:),mc_num,3);
        [db,ci,boot]=ResponseStatistics(err,result.energy,true);
        result.nmse_db(axis,j,:)=db; result.ci(axis,j,:,:)=ci;
        result.bootstrap{axis,j}=boot;
        result.refit_gain_db(axis,j)=db(3)-db(1);
        gains=sort(boot(:,3)-boot(:,1)); result.refit_gain_ci(axis,j,:)=gains([50,1950]);
        ds=result.delay_spread(:,:,axis,j)*1e9; as=result.angle_spread(:,:,axis,j)*180/pi;
        v=[mean(ds,1).',mean(as,1).'];
        [mu,ci]=ResponseStatistics(v,ones(mc_num,1),false);
        result.spread_mean(axis,j,:)=mu; result.spread_ci(axis,j,:,:)=ci;
        p=result.power(:,:,axis,j); weights=p./sum(p,1);
        weighted=[mean(sum(weights.*ds,1)),mean(sum(weights.*as,1))];
        result.spread_weighted_mean(axis,j,:)=weighted;
        result.values(end+1,:)=[axis,grid(j),db,mu,weighted];
    end
end
SaveResult(folder,'variance_compensation',result);
plot_CCMVarianceCompensation(folder);
disp(result.values);
end

function t=fit_trial(data,s,c,e,f,sc,ant,sigma2,grid)
maxNumCompThreads(2);
h=data.h; y=data.y; tau0=data.tau; theta0=data.theta;
L=numel(c.tau);
t=struct('error',zeros(2,4,3),'power',zeros(L,2,4), ...
    'delay_spread',zeros(L,2,4),'angle_spread',zeros(L,2,4));
    t.energy=sum(abs(h).^2);
    t.anchor_tau=tau0; t.anchor_theta=theta0;
    evalc('[p0,d0,a0]=f.EstimateProposedParameters(y,sigma2,s,c,e,tau0,theta0,sc,ant,e.parameter_backend);');
    t.anchor_power=p0; t.anchor_delay_spread=d0; t.anchor_angle_spread=a0;
    for axis=1:2
        for j=1:4
            tau=tau0; theta=theta0; st=c.sigma_tau; sa=c.sigma_theta;
            if axis==1
                tau=c.tau+grid(j)*(tau0-c.tau); st=grid(j)*st;
            else
                theta=c.theta+grid(j)*(theta0-c.theta); sa=grid(j)*sa;
            end
            if grid(j)==1
                p=p0; d=d0; a=a0;
            else
                evalc('[p,d,a]=f.EstimateProposedParameters(y,sigma2,s,c,e,tau,theta,sc,ant,e.parameter_backend);');
            end
            t.power(:,axis,j)=p;
            t.delay_spread(:,axis,j)=d; t.angle_spread(:,axis,j)=a;
            for m=1:3
                pm=p; dm=d; am=a; smt=e.sigma_tau_hat; sma=e.sigma_theta_hat; trunc=0;
                if m==2
                    pm=c.cluster_powers_NoBlockage; dm=c.c_ds_vec; am=c.c_asa_vec;
                    smt=st; sma=sa; trunc=c.SpreadTruncFactor;
                elseif m==3
                    pm=p0; dm=d0; am=a0;
                end
                [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,smt,dm,theta,sma,am,s.f,sc,ant,trunc);
                he=f.EstimateChannelFromFactors(y,pm,sigma2,R,Rh,e.final_solver,e.pcg_opts);
                t.error(axis,j,m)=sum(abs(he-h).^2);
            end
        end
    end
end
