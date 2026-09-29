function exp_InitializationConvergence(mc_num,folder)
%EXP_INITIALIZATIONCONVERGENCE Individually moment-matched spread starts.
if nargin<1, mc_num=20; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results'); end
[s,c,~,e,f]=Base(false); s.Kp=64;
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).'; L=numel(c.tau);
sigma2=10^(-c.SNR_dB/10); scales=[.5,1,1.5];
opts=e.ao_opts; opts.ao_max_iter=30; opts.ao_tol=0;
opts.parameter_backend=e.parameter_backend;
result=struct('sys',s,'channel',c,'estor',e,'opts',opts,'scales',scales, ...
    'seeds',170900+(1:mc_num),'error',zeros(mc_num,31,3), ...
    'nll',zeros(mc_num,31,3),'energy',zeros(mc_num,1), ...
    'power_init',zeros(L,3,mc_num));
for mc=1:mc_num
    rng(result.seeds(mc),'twister'); [h,y,tau,theta]=f.DrawTrial(s,c,f);
    result.energy(mc)=sum(abs(h).^2);
    Z=fft2(reshape(y,s.N,s.Kp))/sqrt(s.N*s.Kp); power=abs(Z(:)).^2;
    for si=1:3
        ds=scales(si)*c.c_ds_ref*ones(L,1); as=scales(si)*c.c_asa_ref*ones(L,1);
        [~,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,ds,theta,c.sigma_theta,as,s.f,sc,ant);
        p=f.MomentMatching(y,R,sigma2); result.power_init(:,si,mc)=p;
        nll0=f.ComputeNLL(p,R,sigma2,power,s.N*s.Kp,e.parameter_backend);
        evalc('[~,~,~,hist]=f.TraceAO(y,sigma2,p,ds,as,tau,c.sigma_tau,c.sigma_theta,theta,s.f,sc,ant,s.N,false,opts);');
        assert(numel(hist.nll)==60);
        result.nll(mc,:,si)=[nll0,hist.nll(2:2:end)];
        for k=0:30
            if k==0, pk=p; dk=ds; ak=as;
            else, pk=hist.p(:,2*k); dk=hist.c_ds_vec{2*k}; ak=hist.c_asa_vec{2*k}; end
            [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,dk,theta,c.sigma_theta,ak,s.f,sc,ant);
            he=f.EstimateChannelFromFactors(y,pk,sigma2,R,Rh,e.final_solver,e.pcg_opts);
            result.error(mc,k+1,si)=sum(abs(he-h).^2);
        end
    end
    fprintf('Initialization trial %d/%d\n',mc,mc_num);
end
result.nmse_db=10*log10(squeeze(sum(result.error,1))/sum(result.energy));
result.mean_nll=squeeze(mean(result.nll,1))/(s.N*s.Kp);
result.values=[(0:30).',result.mean_nll,result.nmse_db];
result.columns={'iteration','nll_0p5','nll_1','nll_1p5','nmse_0p5','nmse_1','nmse_1p5'};
SaveResult(folder,'initialization',result);
end
