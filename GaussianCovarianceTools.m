function f = GaussianCovarianceTools
%GAUSSIANCOVARIANCETOOLS Separable complex-Gaussian sampling and exact norms.
f.Roots=@roots; f.Sample=@sample; f.RelativeError=@relative_error;
end

function B=roots(R)
for l=1:numel(R.tau)
    B.tau{l}=psd_root(R.tau{l});
    B.theta{l}=psd_root(R.theta{l});
end
end

function B=psd_root(A)
A=(A+A')/2; [V,D]=eig(A,'vector'); D=real(D);
% Remove roundoff only, never regularize a materially indefinite covariance.
assert(min(D)>=-1e-10*max(1,max(D)),'Gaussian covariance is not PSD.');
B=V.*sqrt(max(D,0)).';
end

function h=sample(B,p,S)
N=size(B.theta{1},1); K=size(B.tau{1},1);
h=zeros(N*K,S);
for l=1:numel(p)
    for k=1:S
        Z=(randn(N,K)+1j*randn(N,K))/sqrt(2);
        H=sqrt(p(l))*B.theta{l}*Z*B.tau{l}.';
        h(:,k)=h(:,k)+H(:);
    end
end
end

function err=relative_error(p,R,q,T)
den=inner(q,T,q,T);
num=inner(p,R,p,R)+den-2*inner(p,R,q,T);
err=sqrt(max(real(num),0)/real(den));
end

function v=inner(p,R,q,T)
% Full-grid factors are Hermitian Toeplitz. Weight each lag by multiplicity.
F=toeplitz_gram(R.tau,T.tau); A=toeplitz_gram(R.theta,T.theta);
v=real(p(:).'*(F.*A)*q(:));
end

function G=toeplitz_gram(R,T)
n=size(R{1},1); X=zeros(n,numel(R)); Y=zeros(n,numel(T));
for k=1:numel(R), X(:,k)=R{k}(:,1); end
for k=1:numel(T), Y(:,k)=T{k}(:,1); end
G=2*real(X'*((n:-1:1)'.*Y))-n*real(X(1,:)'*Y(1,:));
end
