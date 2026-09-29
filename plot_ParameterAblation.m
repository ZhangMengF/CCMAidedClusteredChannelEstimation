function plot_ParameterAblation(folder)
%PLOT_PARAMETERABLATION Replot five parameter-block controls from saved data.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','response_supplement_20260925'); end
data=load(fullfile(folder,'parameter_ablation.mat'),'result'); r=data.result;
fig=figure('Visible','off','Color','w','Position',[100 100 870 550]);
ax=axes(fig,'Position',[.14 .17 .81 .77]); hold(ax,'on');
colors=lines(5); markers={'d','^','v','s','o'}; curves=gobjects(1,5);
for m=1:5
    y=r.nmse_db(:,m); lo=r.ci(:,1,m); hi=r.ci(:,2,m);
    curves(m)=errorbar(r.grid,y,y-lo,hi-y,['-',markers{m}],'Color',colors(m,:), ...
        'LineWidth',1.7,'MarkerSize',8,'CapSize',7);
end
xlabel('Pilot subcarriers Kp'); ylabel('NMSE (dB)'); xticks(r.grid);
xlim([28 68]); grid on; box on;
ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.12,0];
ax.YLabel.Units='normalized'; ax.YLabel.Position=[-.12,.5,0];
lg=legend(curves,r.methods,'Visible','on','FontSize',15);
lg.Units='normalized'; lg.Position=[.57 .57 .32 .31];
SaveResponseFigure(fig,folder,'parameter_ablation');
end
