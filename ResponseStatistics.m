function [point,ci,boot] = ResponseStatistics(numerator,denominator,in_db)
%RESPONSESTATISTICS Paired percentile bootstrap, independent of simulation RNG.
% Columns are methods/metrics, rows are independent Monte Carlo trials.
if nargin<2, denominator=ones(size(numerator,1),1); end
if nargin<3, in_db=false; end
n=size(numerator,1); B=2000;
stream=RandStream('mt19937ar','Seed',260925);
ix=randi(stream,n,n,B);
den=sum(reshape(denominator(ix),n,B),1).';
boot=zeros(B,size(numerator,2));
for j=1:size(numerator,2)
    v=numerator(:,j);
    boot(:,j)=sum(reshape(v(ix),n,B),1).'./den;
end
point=sum(numerator,1)/sum(denominator);
if in_db, point=10*log10(point); boot=10*log10(boot); end
ordered=sort(boot,1); ci=ordered([50,1950],:);
end
