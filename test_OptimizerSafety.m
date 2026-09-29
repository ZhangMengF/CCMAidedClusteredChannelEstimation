function test_OptimizerSafety
%TEST_OPTIMIZERSAFETY Failure-path tests and real dense/diagonal AO checks.
% Read-only tests: do not regenerate or overwrite published experiment data.
[~,~,~,~,utilities]=Base(false);
q=@(x) deal(.5*sum(x.^2),x);
for direction=[-1,1,NaN,Inf]
    [a,f,g,p,status,reset]=utilities.SafeLineSearch(q,1,direction,.5,1);
    check_step(q,1,.5,1,a,f,g,p);
    assert(strcmp(status,'strong_wolfe') && a>0);
    assert(reset==(direction~=-1));
end
[a,f,g,~,status]=utilities.SafeLineSearch(q,0,0,0,0);
assert(a==0 && f==0 && g==0 && strcmp(status,'stationary'));

% Exhausted bracketing: the returned value must match the evaluated step.
linear=@(x) deal(x,1);
[a,f,g,p,status]=utilities.SafeLineSearch(linear,0,-1,0,1,struct('max_ls',1));
check_step(linear,0,0,1,a,f,g,p);
assert(a==1 && strcmp(status,'armijo'));

% Ten safeguarded zoom steps cannot resolve this steep quadratic.
steep=@(x) deal(.5e12*x.^2,1e12*x);
[a,f,g,p,status]=utilities.SafeLineSearch(steep,1,-1e12,.5e12,1e12);
check_step(steep,1,.5e12,1e12,a,f,g,p);
assert(a>0 && strcmp(status,'armijo'));
[a,f,g,~,status]=utilities.SafeLineSearch(steep,1,-1e12,.5e12,1e12, ...
    struct('max_backtrack',0));
assert(a==0 && f==.5e12 && g==1e12 && strcmp(status,'no_step'));

% Nonfinite objective OR gradient must be rejected, never accepted.
for kind=1:3
    fun=@(x) restricted_objective(x,kind);
    [f0,g0]=fun(1);
    [a,f,g,p]=utilities.SafeLineSearch(fun,1,-g0,f0,g0);
    check_step(fun,1,f0,g0,a,f,g,p); assert(a>0);
end
rejected=false;
try
    utilities.SafeLineSearch(q,1,-1,NaN,1);
catch ex
    rejected=strcmp(ex.identifier,'SafeLineSearch:InvalidStart');
end
assert(rejected);

[s,c,~,e,funcs]=Base(false);
s.N=4; s.Np=4; s.K=16; s.Kp=4;
s.f=(-s.K/2:s.K/2-1)'*s.Deltaf;
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).'; sigma2=10^(-c.SNR_dB/10);
L=numel(c.tau); checked=0;
for snapshots=[1,4]
    rng(5600+snapshots,'twister');
    [~,y,tau,theta]=funcs.DrawTrial(s,c,funcs);
    for j=2:snapshots, [~,y(:,j)]=funcs.DrawTrial(s,c,funcs); end
    for backend={'dense','beam_delay_diag'}
        opts=e.ao_opts; opts.parameter_backend=backend{1};
        opts.ao_max_iter=5; opts.ao_tol=0;
        if strcmp(backend{1},'dense')
            statistic=y*y'/snapshots;
        else
            z=fft2(reshape(y,s.N,s.Kp,snapshots))/sqrt(s.N*s.Kp);
            statistic=reshape(mean(abs(z).^2,3),[],1);
        end
        for scale=[.5,1,1.5,2]
            ds=scale*c.c_ds_ref*ones(L,1); as=scale*c.c_asa_ref*ones(L,1);
            [~,R]=funcs.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau, ...
                ds,theta,c.sigma_theta,as,s.f,sc,ant);
            power=funcs.MomentMatching(y,R,sigma2);
            em_opts=struct('max_iter',3,'tol',1e-12,'parameter_backend',backend{1});
            [~,cost,states]=funcs.TraceEM(y,R,sigma2,power,em_opts);
            for j=1:numel(cost)
                actual=funcs.ComputeNLL(states(:,j),R,sigma2,statistic,s.N*s.Kp,backend{1});
                assert(abs(actual-cost(j))<1e-9*max(1,abs(actual)), ...
                    'EM objective and power history are misaligned.');
            end
            initial=funcs.ComputeNLL(power,R,sigma2,statistic,s.N*s.Kp,backend{1});
            evalc('[~,~,~,history]=funcs.TraceAO(y,sigma2,power,ds,as,tau,c.sigma_tau,c.sigma_theta,theta,s.f,sc,ant,s.N,false,opts);');
            values=[initial,history.nll];
            assert(all(isfinite(values)) && all(diff(values)<=1e-9*max(1,max(abs(values)))));
            for j=1:numel(history.nll)
                [~,Rj]=funcs.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau, ...
                    history.c_ds_vec{j},theta,c.sigma_theta,history.c_asa_vec{j},s.f,sc,ant);
                actual=funcs.ComputeNLL(history.p(:,j),Rj,sigma2,statistic,s.N*s.Kp,backend{1});
                assert(abs(actual-history.nll(j))<1e-9*max(1,abs(actual)));
            end
            checked=checked+1;
        end
    end
end
fprintf('PASS: line-search failures, nonfinite trials, EM record alignment and %d dense/diagonal AO traces.\n',checked);
end

function check_step(fun,x,f0,g0,a,f,g,p)
[actual,gradient]=fun(x+a*p);
assert(isfinite(f) && all(isfinite(g)));
assert(abs(actual-f)<=1e-12*max(1,abs(f)) && norm(gradient-g)<=1e-12*max(1,norm(g)));
assert(f<=f0 && (a==0 || f<=f0+1e-4*a*(g0'*p)));
end

function [f,g]=restricted_objective(x,kind)
f=(x+1)^2; g=2*(x+1);
if x<0
    if kind==1, f=Inf; elseif kind==2, f=NaN; else, g=NaN; end
end
end
