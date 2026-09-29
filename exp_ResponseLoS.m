function exp_ResponseLoS(mc_num,folder)
%EXP_RESPONSELOS Existing CDL-D channel and Kp scan; no new channel algorithm.
if nargin<1, mc_num=100; end
if nargin<2, folder=fullfile(fileparts(mfilename('fullpath')),'results','response_supplement_20260925'); end
[s,c,sim,e,f]=Base(true); c.SNR_dB=5; sim.MCNum=mc_num;
result=MonteCarloNMSE_ScanParas(s,c,sim,e,f,0, ...
    [false false true true true true true],'Kp',[32,64],fullfile(folder,'los_scan'));
% Keep the existing plot, only simplify its layout for the response.
fig=gcf; ax=gca; ax.Units='normalized'; ax.Position=[.14 .17 .80 .73];
title(sprintf('CDL-D-based LoS | 5 dB | %d paired trials',mc_num));
xticks([32 64]); xlim([28 68]); ylim([-17 -2]);
ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.13,0];
ax.YLabel.Units='normalized'; ax.YLabel.Position=[-.12,.5,0];
lg=findall(fig,'Type','legend'); lg.Units='normalized'; lg.Position=[.55 .49 .36 .32];
SaveResponseFigure(fig,fullfile(folder,'los_scan'),'scan_Kp');
result.methods={'CPM-aided','Fixed-profile','Proposed','Oracle','Pilot-only DFT'};
result.ci=zeros(2,2,5); result.bootstrap=cell(2,1);
for j=1:2
    [~,ci,b]=ResponseStatistics(result.trials{j}.error(:,3:7),result.trials{j}.energy,true);
    result.ci(j,:,:)=ci; result.bootstrap{j}=b;
end
result.source=struct('Base',fileread(which('Base')),'entry',fileread([mfilename('fullpath'),'.m']));
SaveResult(folder,'response_los',result);
end
