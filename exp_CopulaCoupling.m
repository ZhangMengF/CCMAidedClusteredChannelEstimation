function exp_CopulaCoupling(mc_num,folder)
%EXP_COPULACOUPLING Gaussian dependence of delay / |angle|, fixed marginals.
% Equal-power random rays; Oracle knows the population joint PDAP, not rays.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results'); end
[s,c,~,e,f]=Base(false); grid=0:.2:.8;
result=struct('sys',s,'channel',c,'estor',e,'grid',grid,'seed',20260923, ...
    'mc_num',mc_num,'methods',{{'Oracle','Fixed-profile','Proposed'}}, ...
    'error',zeros(mc_num,3,5),'energy',zeros(mc_num,5), ...
    'nmse_db',nan(5,3),'quadrature_relative_error',zeros(1,5));
result.values=nan(5,6);
result.columns={'r','oracle_dB','fixed_dB','proposed_dB','proposed_gap_dB','fixed_gap_dB'};
for j=1:5
    [err,energy,qerror]=one_point(grid(j),mc_num,s,c,e,f);
    db=10*log10(sum(err,1)/sum(energy));
    result.error(:,:,j)=err; result.energy(:,j)=energy;
    result.nmse_db(j,:)=db; result.quadrature_relative_error(j)=qerror;
    result.values(j,:)=[grid(j),db,db(3)-db(1),db(2)-db(1)];
    SaveResult(folder,'copula',result);
    fprintf('Copula r=%.1f NMSE=%s dB\n',grid(j),mat2str(db,5));
end
end

function [err,energy,qerror]=one_point(r,mc_num,s,c,e,f)
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).';
idx=f.PilotObservationIndices(sc,ant,s.N); sigma2=10^(-c.SNR_dB/10);
L=numel(c.tau); M=c.M_subpath;
fixed=exp(-c.tau/c.DS); fixed=fixed/sum(fixed);
df=(-(s.K-1):s.K-1)*s.Deltaf; conditional_F=cell(2,L);
if r>0
    for level=1:2
        nq=128*level; [u,~]=uniform_quadrature(nq);
        za=sqrt(2)*erfinv(2*u-1); [z,w]=normal_quadrature(nq);
        ud=0.5*erfc(-(r*za(:)+sqrt(1-r^2)*z.')/sqrt(2));
        for l=1:L
            delays=delay_quantile(ud,c.c_ds_vec(l)); F=zeros(nq,numel(df));
            for k=1:nq, F=F+w(k)*exp(-1j*2*pi*delays(:,k)*df); end
            conditional_F{level,l}=F;
        end
    end
end
err=zeros(mc_num,3); energy=zeros(mc_num,1); qerror=0;
rng(20260923,'twister'); timer=tic;
for mc=1:mc_num
    h=zeros(s.N*s.K,1);
    for l=1:L
        x=randn(M,1); y=randn(M,1);
        u=0.5*erfc(-x/sqrt(2)); v=0.5*erfc(-(r*x+sqrt(1-r^2)*y)/sqrt(2));
        dt=delay_quantile(u,c.c_ds_vec(l));
        da=sign(rand(M,1)-.5).*magnitude_quantile(v,c.c_asa_vec(l));
        a=exp(1j*pi*(0:s.N-1).'*sin(c.theta(l)+da.'));
        b=exp(-1j*2*pi*s.f*(c.tau(l)+dt.'));
        gain=sqrt(c.cluster_powers_NoBlockage(l)/M)*exp(1j*2*pi*rand(M,1));
        h=h+reshape((a.*gain.')*b.',[],1);
    end
    noise=(randn(numel(idx),1)+1j*randn(numel(idx),1))/sqrt(2);
    yp=h(idx)+sqrt(sigma2)*noise;
    theta=c.theta+c.sigma_theta*randn(L,1);
    te=c.sigma_tau*randn(L,1); te(1)=0; tau=c.tau+te;
    kernel=joint_kernel(r,128,conditional_F(1,:),s,c,f,tau,theta);
    ho=kernel_estimate(kernel,yp,s.N,s.K,idx,sigma2);
    if mc==1
        if r>0
            k2=joint_kernel(r,256,conditional_F(2,:),s,c,f,tau,theta);
            h2=kernel_estimate(k2,yp,s.N,s.K,idx,sigma2);
            qerror=norm(ho-h2)/norm(h2);
            assert(qerror<1e-3,'Oracle quadrature requires refinement');
        else
            [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_vec,theta,c.sigma_theta,c.c_asa_vec,s.f,sc,ant,2);
            h2=f.EstimateChannelFromFactors(yp,c.cluster_powers_NoBlockage,sigma2,R,Rh,e.final_solver,e.pcg_opts);
            qerror=norm(ho-h2)/norm(h2);
            assert(qerror<1e-6,'Zero-coupling Oracle mismatch');
        end
    end
    [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_ref*ones(L,1),theta,c.sigma_theta,c.c_asa_ref*ones(L,1),s.f,sc,ant);
    hf=f.EstimateChannelFromFactors(yp,fixed,sigma2,R,Rh,e.final_solver,e.pcg_opts);
    evalc('[p,ds,as]=f.EstimateProposedParameters(yp,sigma2,s,c,e,tau,theta,sc,ant,e.parameter_backend);');
    [Rh,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,ds,theta,c.sigma_theta,as,s.f,sc,ant);
    hp=f.EstimateChannelFromFactors(yp,p,sigma2,R,Rh,e.final_solver,e.pcg_opts);
    err(mc,:)=[sum(abs(h-ho).^2),sum(abs(h-hf).^2),sum(abs(h-hp).^2)];
    energy(mc)=sum(abs(h).^2);
    if mod(mc,10)==0 || mc==mc_num
        fprintf('Copula r=%.1f trial=%d/%d elapsed=%.1fs\n',r,mc,mc_num,toc(timer));
    end
end
end

function kernel=joint_kernel(r,nq,F,s,c,f,tau,theta)
dn=(-(s.N-1):s.N-1).'; df=(-(s.K-1):s.K-1)*s.Deltaf;
kernel=zeros(2*s.N-1,2*s.K-1); L=numel(tau);
if r==0
    [~,R]=f.BuildClusterCovarianceFactors(false,s.N,tau,c.sigma_tau,c.c_ds_vec,theta,c.sigma_theta,c.c_asa_vec,s.f,(1:s.K).',(1:s.N).',2);
    for l=1:L
        a=R.theta{l}(:,1); b=R.tau{l}(:,1);
        kernel=kernel+c.cluster_powers_NoBlockage(l)*[conj(a(end:-1:2));a]*[conj(b(end:-1:2));b].';
    end
else
    [u,w]=uniform_quadrature(nq); [z,wz]=normal_quadrature(12);
    for l=1:L
        a=magnitude_quantile(u,c.c_asa_vec(l)); A=zeros(numel(dn),nq);
        for g=1:numel(z)
            A=A+0.5*wz(g)*(exp(1j*pi*dn*sin(theta(l)+a+c.sigma_theta*z(g))) ...
                +exp(1j*pi*dn*sin(theta(l)-a+c.sigma_theta*z(g))));
        end
        st=c.sigma_tau*(l~=1);
        C=(A.*w)*F{l}.*exp(-1j*2*pi*tau(l)*df-0.5*(2*pi*st*df).^2);
        kernel=kernel+c.cluster_powers_NoBlockage(l)*C;
    end
end
end

function d=delay_quantile(u,c)
d=-c*log(1-u*(1-exp(-2)));
end
function a=magnitude_quantile(u,c)
a=-c/sqrt(2)*log(1-u*(1-exp(-2*sqrt(2))));
end
function [x,w]=normal_quadrature(n)
[V,D]=eig(diag(sqrt(1:n-1),1)+diag(sqrt(1:n-1),-1));
[x,ix]=sort(diag(D));w=V(1,ix).^2;
end
function [u,w]=uniform_quadrature(n)
% Gauss-Legendre split at the Laplace median to integrate its cusp.
k=1:n/2-1; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
[x,ix]=sort(diag(D)); v=V(1,ix).^2;
u=[(x.'+1)/4,(x.'+3)/4]; w=[v/2,v/2];
end
function h=kernel_estimate(kernel,y,N,K,idx,sigma2)
na=2^nextpow2(2*N-1); nf=2^nextpow2(2*K-1);
E=zeros(na,nf);
E(mod(-(N-1):N-1,na)+1,mod(-(K-1):K-1,nf)+1)=kernel;
S=fft2(E);
apply=@(v) pilot_apply(v,S,N,K,idx,sigma2);
[x,flag,relres]=pcg(apply,y,1e-8,5000);
assert(flag==0,'Joint Oracle PCG failed: %g',relres);
z=zeros(N*K,1);z(idx)=x;
h=full_apply(z,S,N,K);
end
function y=pilot_apply(x,S,N,K,idx,sigma2)
z=zeros(N*K,1);z(idx)=x;
z=full_apply(z,S,N,K);y=z(idx)+sigma2*x;
end
function y=full_apply(x,S,N,K)
z=ifft2(fft2(reshape(x,N,K),size(S,1),size(S,2)).*S);
y=reshape(z(1:N,1:K),[],1);
end
