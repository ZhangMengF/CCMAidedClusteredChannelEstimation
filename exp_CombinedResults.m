function exp_CombinedResults(folder,workers)
%EXP_COMBINEDRESULTS Reproduce all Fig.1 data with accelerated AO and FFT-PCG.
% Same grids, trial counts and seeds as the three original scan entries.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','combined_accelerated'); end
if nargin<2, workers=8; end
if ~isfolder(folder), mkdir(folder); end
[s,c,sim,e,f]=Base(false);
e.parameter_backend='beam_delay_diag'; e.final_solver='fft_pcg';
e.pilot_only_method='dft';
source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']));
sourcefile=fullfile(folder,'source.mat');
if isfile(sourcefile)
    old=load(sourcefile,'source'); assert(isequal(old.source,source),'Use a fresh folder after source changes.');
else
    save(sourcefile,'source');
end
snrs=[-5 0 5 10 15]; kps=[8 16 32 64 128]; histories=[0 1 3 7 15]; blocks=[.5 .3 0];
% Columns: panel, scan index, history index, SNR, Kp, blockage, S, MC, seed.
jobs=zeros(25,9);
for j=1:5
    jobs(j,:)=[1 j 0 snrs(j) 32 0 0 100 100+j];
    jobs(5+j,:)=[2 j 0 5 kps(j) 0 0 50 100+j];
end
for j=1:3
    for i=1:5
        jobs(10+5*(j-1)+i,:)=[3 j i 5 64 blocks(j) histories(i) 100 1000*j+100*i];
    end
end
pool=gcp('nocreate');
if workers>0 && isempty(pool), pool=parpool('Processes',workers); end
points=cell(25,1);
parfor (j=1:25,workers)
    maxNumCompThreads(2);
    file=fullfile(folder,sprintf('point_%02d.mat',j));
    if isfile(file)
        saved=load(file,'point'); point=saved.point;
        assert(isequal(point.job,jobs(j,:)));
    else
        point=run_point(s,c,sim,e,f,jobs(j,:));
        save_point(file,point);
    end
    points{j}=point;
    fprintf('Fig.1 configuration %d/25 complete\n',j);
end
for panel=1:2
    ids=(panel-1)*5+(1:5); values=zeros(5,6); trials=cell(5,1);
    for j=1:5
        p=points{ids(j)}; values(j,:)=[jobs(ids(j),3+panel),p.nmse_db(3:7)]; trials{j}=p.trials;
    end
    names={'SNR_dB','Kp'};
    result=struct('sys',s,'channel',c,'estor',e,'jobs',jobs(ids,:), ...
        'source',source,'trials',{trials},'values',values, ...
        'columns',{{names{panel},'CPM','Fixed_profile','Proposed','Oracle','Pilot_only'}});
    SaveResult(folder,['scan_',names{panel}],result);
end
values=zeros(5,4); values(:,1)=histories(:); trials=cell(5,3);
for j=1:3
    for i=1:5
        p=points{10+5*(j-1)+i}; values(i,j+1)=p.nmse_db(5); trials{i,j}=p.trials;
    end
end
s.Kp=64;
result=struct('sys',s,'channel',c,'estor',e,'source',source,'trials',{trials}, ...
    'S_values',histories,'scan_values',blocks,'jobs',jobs(11:25,:), ...
    'values',values,'columns',{{'S','block_0p5','block_0p3','block_0'}});
SaveResult(folder,'history_BlockProb',result);
plot_CombinedResults(folder);
if workers>0, delete(pool); end
end

function save_point(file,point)
save(file,'point');
end

function point=run_point(s,c,sim,e,f,job)
c.SNR_dB=job(4); s.Kp=job(5); c.ClusterBlockProb=job(6); sim.MCNum=job(8);
switches=[false false true true true true true];
if job(1)==3, switches=[false false false false true false false]; end
rng(job(9),'twister');
evalc('[~,~,~,~,~,~,~,~,~,trials]=f.CalcuNMSEsByMonteCarlo(s,c,sim,e,1:s.K/s.Kp:s.K,(1:s.Np).'',job(7),switches);');
point=struct('job',job,'trials',trials,'nmse_db',10*log10(sum(trials.error,1)/sum(trials.energy)));
end
