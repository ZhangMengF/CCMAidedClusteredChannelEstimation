function exp_CCMPriorMismatch(mc_num,folder)
%EXP_CCMPRIORMISMATCH Actual versus assumed Gaussian CCM error scale.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','response_supplement_20260925'); end
[s,c,sim,e,f]=Base(false); sim.MCNum=mc_num;
sc=(1:s.K/s.Kp:s.K).'; ant=(1:s.N).';
grid=[.5,1,2,4];
result=struct('sys',s,'channel',c,'estor',e,'mc_num',mc_num, ...
    'grid',grid,'axes',{{'Delay uncertainty','Angle uncertainty'}}, ...
    'methods',{{'CPM-aided','Fixed-profile','Proposed','Oracle','Pilot-only DFT'}}, ...
    'seed',260925,'trials',{cell(2,4)},'nmse_db',zeros(2,4,5), ...
    'ci',zeros(2,4,2,5),'values',[], ...
    'columns',{{'axis','scale','CPM','Fixed_profile','Proposed','Oracle','Pilot_only'}});
result.source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']));
for axis=1:2
    for j=1:numel(grid)
        ch=c;
        if axis==1, ch.sigma_tau=c.sigma_tau*grid(j);
        else, ch.sigma_theta=c.sigma_theta*grid(j); end
        % Same draws at every scale; estimators retain e.sigma_*_hat.
        rng(result.seed,'twister');
        fprintf('Prior mismatch: axis %d, actual/assumed SD %.1f\n',axis,grid(j));
        [~,~,~,~,~,~,~,~,~,trial]=f.CalcuNMSEsByMonteCarlo( ...
            s,ch,sim,e,sc,ant,0,[false false true true true true true]);
        [db,ci,boot]=ResponseStatistics(trial.error(:,3:7),trial.energy,true);
        result.trials{axis,j}=trial; result.bootstrap{axis,j}=boot;
        result.nmse_db(axis,j,:)=db; result.ci(axis,j,:,:)=ci;
        result.values(end+1,:)=[axis,grid(j),db];
        SaveResult(folder,'prior_mismatch',result);
    end
end
% Pilot-only and true channel are independent of the map uncertainty.
anchor=result.trials{1,1};
for axis=1:2
    for j=1:numel(grid)
        t=result.trials{axis,j};
        assert(isequal(t.energy,anchor.energy) && isequal(t.error(:,7),anchor.error(:,7)));
    end
end
assert(isequaln(result.trials{1,2}.error,result.trials{2,2}.error));
end
