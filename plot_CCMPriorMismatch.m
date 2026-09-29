function plot_CCMPriorMismatch(folder)
%PLOT_CCMPRIORMISMATCH Replot saved actual/assumed uncertainty scans only.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','response_supplement_20260925'); end
data=load(fullfile(folder,'prior_mismatch.mat'),'result'); r=data.result;
fig=figure('Visible','off','Color','w','Position',[100 100 1200 500]);
colors=[0 .447 .741;.494 .184 .556;.85 .325 .098;.929 .694 .125;.35 .35 .35];
markers={'^','d','s','o','x'}; positions=[.065 .20 .32 .73;.455 .20 .32 .73];
for axis=1:2
    ax=axes(fig,'Position',positions(axis,:)); hold(ax,'on');
    curves=gobjects(1,5);
    for m=1:5
        y=reshape(r.nmse_db(axis,:,m),1,[]);
        lo=reshape(r.ci(axis,:,1,m),1,[]); hi=reshape(r.ci(axis,:,2,m),1,[]);
        curves(m)=errorbar(r.grid,y,y-lo,hi-y,['-',markers{m}],'Color',colors(m,:), ...
            'LineWidth',1.6,'MarkerSize',7,'CapSize',5);
    end
    xlabel('Actual / assumed error SD'); ylabel('NMSE (dB)');
    title(r.axes{axis}); xticks(r.grid); grid on; box on; xlim([.3 4.2]);
    ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.16,0];
    ax.YLabel.Units='normalized'; ax.YLabel.Position=[-.145,.5,0];
    if axis==1
        lg=legend(curves,r.methods,'NumColumns',1,'Visible','on','FontSize',15);
        lg.Units='normalized'; lg.Position=[.80 .35 .19 .38];
    end
end
SaveResponseFigure(fig,folder,'prior_mismatch');
end
