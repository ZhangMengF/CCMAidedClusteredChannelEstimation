function exp_ParameterAblation(mc_num,folder)
%EXP_PARAMETERABLATION Freeze power or spread blocks without retuning updates.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','response_supplement_20260925'); end
[s,c,~,e,f]=Base(false); grid=[32,64]; L=numel(c.tau);
sigma2=10^(-c.SNR_dB/10); ant=(1:s.N).';
result=struct('sys',s,'channel',c,'estor',e,'mc_num',mc_num,'grid',grid, ...
    'methods',{{'Fixed-profile','Power only','Spread only','Proposed','Oracle'}}, ...
    'seeds',270925+(1:mc_num),'error',zeros(mc_num,2,5), ...
    'energy',zeros(mc_num,2),'power',zeros(L,mc_num,2,5), ...
    'delay_spread',zeros(L,mc_num,2,5),'angle_spread',zeros(L,mc_num,2,5), ...
    'power_init',zeros(L,mc_num,2),'nmse_db',zeros(2,5),'ci',zeros(2,2,5), ...
    'values',[],'columns',{{'Kp','Fixed_profile','Power_only','Spread_only','Proposed','Oracle'}});
result.source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']));
for j=1:2
    s.Kp=grid(j); sc=(1:s.K/s.Kp:s.K).';
    for mc=1:mc_num
        rng(result.seeds(mc),'twister'); [h,y,tau,theta]=f.DrawTrial(s,c,f);
        result.energy(mc,j)=sum(abs(h).^2);
        ds=c.c_ds_ref*ones(L,1); as=c.c_asa_ref*ones(L,1);
        [~,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,ds,theta,e.sigma_theta_hat,as,s.f,sc,ant);
        p0=f.MomentMatching(y,R,sigma2); result.power_init(:,mc,j)=p0;
        for m=1:5
            p=c.cluster_powers_emp; d=ds; a=as; rt=0;
            if m>=2 && m<=4
                opts=e.ao_opts; opts.parameter_backend=e.parameter_backend;
                opts.fit_powers=(m~=3); opts.fit_spreads=(m~=2);
                if opts.fit_powers, p=p0; end
                hist=struct;
                evalc('[p,d,a,hist]=f.TraceAO(y,sigma2,p,d,a,tau,e.sigma_tau_hat,e.sigma_theta_hat,theta,s.f,sc,ant,s.N,false,opts);');
                if m==2, assert(isequal(d,ds)&&isequal(a,as)); end
                if m==3, assert(isequal(p,c.cluster_powers_emp)); end
                result.ao_history{mc,j,m}=hist;
            elseif m==5
                p=c.cluster_powers_NoBlockage; d=c.c_ds_vec; a=c.c_asa_vec; rt=c.SpreadTruncFactor;
            end
            [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,d,theta,e.sigma_theta_hat,a,s.f,sc,ant,rt);
            he=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,e.final_solver,e.pcg_opts);
            result.error(mc,j,m)=sum(abs(he-h).^2);
            result.power(:,mc,j,m)=p; result.delay_spread(:,mc,j,m)=d; result.angle_spread(:,mc,j,m)=a;
        end
        if mod(mc,10)==0 || mc==mc_num, fprintf('Ablation Kp=%d: %d/%d\n',s.Kp,mc,mc_num); end
    end
    err=reshape(result.error(:,j,:),mc_num,5);
    [db,ci,boot]=ResponseStatistics(err,result.energy(:,j),true);
    result.nmse_db(j,:)=db; result.ci(j,:,:)=ci; result.bootstrap{j}=boot;
    result.gain_db(j,:)=db(1)-db;
    b=sort(boot(:,1)-boot,1); result.gain_ci(j,:,:)=b([50,1950],:);
    result.values(end+1,:)=[s.Kp,db];
    SaveResult(folder,'parameter_ablation',result);
end
end
