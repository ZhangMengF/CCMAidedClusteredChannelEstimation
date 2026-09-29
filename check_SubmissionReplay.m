function checks=check_SubmissionReplay(folder)
%CHECK_SUBMISSIONREPLAY Compare fresh audit runs with saved submitted evidence.
root=fullfile(fileparts(mfilename('fullpath')),'results');
if nargin<1, folder=fullfile(root,'release_audit'); end
original={'initialization.mat','pilot_only.mat', ...
    'response_supplement_20260925/response_los.mat', ...
    'response_supplement_20260925/parameter_ablation.mat', ...
    'parameter_accuracy_truncated/parameter_accuracy.mat','copula.mat'};
fresh={'initialization.mat','pilot_only.mat','response/response_los.mat', ...
    'response/parameter_ablation.mat','accuracy/parameter_accuracy.mat','copula.mat'};
checks=struct;
for j=1:numel(original)
    if ~isfile(fullfile(folder,fresh{j})), continue; end
    x=load(fullfile(root,original{j}),'result'); a=x.result;
    x=load(fullfile(folder,fresh{j}),'result'); b=x.result;
    [~,name]=fileparts(original{j});
    if isfield(a,'trials')
        av=cellfun(@(t) t.error,a.trials,'UniformOutput',false);
        bv=cellfun(@(t) t.error,b.trials,'UniformOutput',false);
    else
        av=a.error; bv=b.error;
    end
    [delta,scale]=compare_prefix(av,bv);
    checks.(name)=struct('max_absolute_error_difference',delta, ...
        'relative_to_max_error',delta/max(scale,eps),'old_source',original{j},'replay',fresh{j});
    if isfield(a,'nmse_db') && isequal(size(av),size(bv)) && ~iscell(av)
        checks.(name).max_nmse_difference_dB=max(abs(a.nmse_db-b.nmse_db),[],'all');
    end
    fprintf('%s: raw squared-error difference %.6g (relative %.6g)\n',name,delta,delta/max(scale,eps));
end
save(fullfile(folder,'replay_checks.mat'),'checks');
end

function [delta,scale]=compare_prefix(a,b)
if iscell(a)
    delta=0; scale=0;
    for k=1:numel(a)
        [d,s]=compare_prefix(a{k},b{k}); delta=max(delta,d); scale=max(scale,s);
    end
else
    dims=repmat({':'},1,ndims(a)); dims{1}=1:size(b,1); a=a(dims{:});
    mask=isfinite(a)&isfinite(b);
    assert(isequal(isfinite(a),isfinite(b)),'Nonfinite pattern changed.');
    delta=max(abs(a(mask)-b(mask)),[],'all'); scale=max(abs(a(mask)),[],'all');
end
end
