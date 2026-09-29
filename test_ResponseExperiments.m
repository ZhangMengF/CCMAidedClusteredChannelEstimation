function test_ResponseExperiments
%TEST_RESPONSEEXPERIMENTS Compatibility, frozen blocks and Gaussian covariance.
[s,c,sim,e,f]=Base(false); s.N=4; s.Np=4; s.K=16; s.Kp=4;
s.f=(-s.K/2:s.K/2-1)'*s.Deltaf; sim.MCNum=1;
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).'; L=numel(c.tau); sigma2=10^(-c.SNR_dB/10);
rng(419,'twister'); [h,y,tau,theta]=f.DrawTrial(s,c,f);
d=c.c_ds_ref*ones(L,1); a=c.c_asa_ref*ones(L,1);
[~,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,d,theta,c.sigma_theta,a,s.f,sc,ant);
p=f.MomentMatching(y,R,sigma2); opts=e.ao_opts; opts.parameter_backend=e.parameter_backend;
evalc('[p0,d0,a0,h0]=f.TraceAO(y,sigma2,p,d,a,tau,c.sigma_tau,c.sigma_theta,theta,s.f,sc,ant,s.N,false,opts);');
opts.fit_powers=true; opts.fit_spreads=true;
evalc('[p1,d1,a1,h1]=f.TraceAO(y,sigma2,p,d,a,tau,c.sigma_tau,c.sigma_theta,theta,s.f,sc,ant,s.N,false,opts);');
assert(isequal(p0,p1)&&isequal(d0,d1)&&isequal(a0,a1)&&isequal(h0,h1));
for fit=[false true;true false;false false].'
    opts.fit_powers=fit(1); opts.fit_spreads=fit(2);
    evalc('[ph,dh,ah]=f.TraceAO(y,sigma2,p,d,a,tau,c.sigma_tau,c.sigma_theta,theta,s.f,sc,ant,s.N,false,opts);');
    if ~opts.fit_powers, assert(isequal(p,ph)); end
    if ~opts.fit_spreads, assert(isequal(d,dh)&&isequal(a,ah)); end
end
% Actual-error changes do not modify pilot-only data; assumed-error changes
% do affect the map-aware estimators without changing the generated channel.
flags=[false false true true true true true];
rng(420,'twister'); evalc('[~,~,~,~,~,~,~,~,~,t0]=f.CalcuNMSEsByMonteCarlo(s,c,sim,e,sc,ant,0,flags);');
e2=e; e2.sigma_tau_hat=2*e.sigma_tau_hat; e2.sigma_theta_hat=2*e.sigma_theta_hat;
rng(420,'twister'); evalc('[~,~,~,~,~,~,~,~,~,t1]=f.CalcuNMSEsByMonteCarlo(s,c,sim,e2,sc,ant,0,flags);');
assert(isequal(t0.energy,t1.energy)&&isequal(t0.error(:,6:7),t1.error(:,6:7)));
assert(any(abs(t0.error(:,3:5)-t1.error(:,3:5))>1e-8));
g=GaussianCovarianceTools;
[~,Rf]=f.BuildClusterCovarianceFactors(false,s.N,c.tau,c.sigma_tau,c.c_ds_vec,c.theta,c.sigma_theta,c.c_asa_vec,s.f,(1:s.K).',ant);
[~,Tf]=f.BuildClusterCovarianceFactors(false,s.N,c.tau,c.sigma_tau,d,c.theta,c.sigma_theta,a,s.f,(1:s.K).',ant);
q=c.cluster_powers_NoBlockage; p=c.cluster_powers_emp;
C=zeros(s.N*s.K); D=C;
for l=1:L
    C=C+q(l)*kron(Rf.tau{l},Rf.theta{l});
    D=D+p(l)*kron(Tf.tau{l},Tf.theta{l});
end
err=g.RelativeError(p,Tf,q,Rf);
assert(abs(err-norm(D-C,'fro')/norm(C,'fro'))<1e-10);
B=g.Roots(Rf);
for l=1:L
    assert(norm(B.tau{l}*B.tau{l}'-Rf.tau{l},'fro')/norm(Rf.tau{l},'fro')<1e-10);
end
rng(421,'twister'); H=g.Sample(B,q,5000);
empirical_error=norm(H*H'/size(H,2)-C,'fro')/norm(C,'fro');
assert(empirical_error<.08);
before=rng; ResponseStatistics(abs(H(1:2,:).').^2); assert(isequal(before,rng));
fprintf('PASS: default AO, frozen blocks, assumed/actual priors, Gaussian sampling (%.3g), covariance norm and isolated bootstrap RNG.\n',empirical_error);
end
