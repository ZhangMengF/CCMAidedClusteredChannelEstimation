function test_Numerics
%TEST_NUMERICS Fast deterministic checks; no external data/toolboxes.
[s,c,sim,e,f]=Base(false); s.Kp=8; sim.MCNum=1;
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).';
rng(981,'twister'); [h,y,tau,theta]=f.DrawTrial(s,c,f);
rng(981,'twister'); [h2,y2,tau2,theta2]=f.DrawTrial(s,c,f);
assert(isequal(h,h2)&&isequal(y,y2)&&isequal(tau,tau2)&&isequal(theta,theta2));
sigma2=10^(-c.SNR_dB/10); p=c.cluster_powers_NoBlockage;
[Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_vec,theta,c.sigma_theta,c.c_asa_vec,s.f,sc,ant);
A=sigma2*eye(numel(y));
for l=1:numel(p), A=A+p(l)*kron(R.tau{l},R.theta{l}); end
S=f.BuildFFTOperator(p,sigma2,R);
fft_error=norm(f.ApplyFFTOperator(y,S,s.N,s.Kp)-A*y)/norm(A*y);
assert(fft_error<1e-12);
hp=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'fft_pcg',e.pcg_opts);
hd=f.EstimateChannelFromFactors(y,p,sigma2,R,Rh,'dense_chol',e.pcg_opts);
solver_error=norm(hp-hd)/norm(hd); assert(solver_error<1e-6);
methods={'linear','quadratic','fir','pchip','spline','dft'};
for m=1:numel(methods)
    he=reshape(f.PilotOnlyEstimate(y,s.K,s.N,sc,methods{m}),s.N,s.K);
    z=he(:,sc); assert(norm(z(:)-y)/norm(y)<1e-12);
    z=f.PilotOnlyEstimate(ones(size(y)),s.K,s.N,sc,methods{m});
    assert(max(abs(z-1))<1e-12);
    z=f.PilotOnlyEstimate(h,s.K,s.N,(1:s.K).',methods{m});
    assert(norm(z-h)/norm(h)<1e-12);
    % The scan's keyword must reach the Monte Carlo estimator unchanged.
    e.pilot_only_method=upper(methods{m}); rng(981,'twister');
    evalc('[~,~,~,~,~,~,~,~,pilot_nmse]=f.CalcuNMSEsByMonteCarlo(s,c,sim,e,sc,ant,0,[false,false,false,false,false,false,true]);');
    direct=f.PilotOnlyEstimate(y,s.K,s.N,sc,methods{m});
    assert(abs(pilot_nmse-sum(abs(direct-h).^2)/sum(abs(h).^2))<1e-12);
end
for history=[0,1]
    rng(701,'twister');
    evalc('[~,~,~,~,nmse]=f.CalcuNMSEsByMonteCarlo(s,c,sim,e,sc,ant,history,[false,false,false,false,true,false]);');
    assert(isfinite(nmse)&&nmse>0);
end
% Copula preserves one-sided exponential and symmetric Laplacian marginals.
rng(917,'twister'); M=100000; X=randn(M,1); Y=randn(M,1);
u=.5*erfc(-X/sqrt(2)); sign_a=2*(rand(M,1)>.5)-1;
for r=[0,.4,.8,1]
    v=.5*erfc(-(r*X+sqrt(1-r^2)*Y)/sqrt(2));
    dt=-log(1-u*(1-exp(-2)));
    da=-sign_a/sqrt(2).*log(1-v*(1-exp(-2*sqrt(2))));
    assert(all(dt>=0 & dt<=2) && all(abs(da)<=2));
    assert(abs(mean(da))<.01);
    assert(max(abs(sort(v)-(1:M)'/M))<.01);
    if r==1
        assert(max(abs(v-u))<1e-12);
    end
end
fprintf('PASS: FFT matvec relative error %.3g; FFT-PCG estimate error %.3g; interpolation, seeds, history and copula checks passed.\n',fft_error,solver_error);
end
