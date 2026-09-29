function inspect_CCMVarianceCompensation(folder)
%INSPECT_CCMVARIANCECOMPENSATION Replay selected extreme-spread trials, unchanged AO.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','ccm_variance_compensation'); end
load(fullfile(folder,'variance_compensation.mat'),'result'); r=result;
[~,~,~,~,f]=Base(false); s=r.sys; c=r.channel; e=r.estor; L=numel(c.tau);
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).'; sigma2=10^(-c.SNR_dB/10);
cases=[97,1,1;42,1,4;60,2,4;64,2,4]; traces=cell(4,1);
rng(r.seed,'twister');
for mc=1:max(cases(:,1))
    [~,y,tau0,theta0]=f.DrawTrial(s,c,f);
    for k=find(cases(:,1)==mc).'
        axis=cases(k,2); a=cases(k,3); tau=tau0; theta=theta0;
        if axis==1, tau=c.tau+a*(tau0-c.tau); else, theta=c.theta+a*(theta0-c.theta); end
        ds=c.c_ds_ref*ones(L,1); as=c.c_asa_ref*ones(L,1);
        [~,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,e.sigma_tau_hat,ds,theta,e.sigma_theta_hat,as,s.f,sc,ant);
        p0=f.MomentMatching(y,R,sigma2); opts=e.ao_opts; opts.parameter_backend=e.parameter_backend;
        fprintf('\nReplay trial %d, axis %d, a=%g\n',mc,axis,a);
        [p,d,ang,hist]=f.TraceAO(y,sigma2,p0,ds,as,tau,e.sigma_tau_hat,e.sigma_theta_hat,theta,s.f,sc,ant,s.N,false,opts);
        j=find(r.grid==a);
        assert(norm(d-r.delay_spread(:,mc,axis,j))/norm(d)<1e-6);
        assert(norm(ang-r.angle_spread(:,mc,axis,j))/norm(ang)<1e-6);
        traces{k}=struct('power',p,'delay',d,'angle',ang,'history',hist);
    end
end
save(fullfile(folder,'extreme_spread_diagnostics.mat'),'cases','traces');
end
