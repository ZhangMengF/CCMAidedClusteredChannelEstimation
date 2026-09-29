function plot_ParameterAccuracy(folder)
%PLOT_PARAMETERACCURACY Initial/fitted accuracy; independently evaluated NMSE.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','parameter_accuracy_truncated'); end
data=load(fullfile(folder,'parameter_accuracy.mat'),'result'); r=data.result;
fig=figure('Visible','off','Color','w','Position',[100 100 1200 900]);
colors=lines(2); styles={'--o','-s'};
positions=[.10 .60 .36 .33;.61 .60 .36 .33;.10 .16 .36 .33;.61 .16 .36 .33];
labels=cell(1,4);
for metric=1:4
    ax=axes(fig,'Position',positions(metric,:)); hold(ax,'on'); curves=gobjects(1,4);
    for j=1:2
        for stage=1:2
            y=reshape(r.metric_mean(j,:,stage,metric),1,[]);
            lo=reshape(r.metric_ci(j,:,1,stage,metric),1,[]);
            hi=reshape(r.metric_ci(j,:,2,stage,metric),1,[]);
            idx=2*(j-1)+stage; labels{idx}=sprintf('Kp=%d, %s',r.grid(j),r.methods{stage});
            curves(idx)=errorbar([1,2],y,y-lo,hi-y,styles{stage},'Color',colors(j,:), ...
                'LineWidth',2.3,'MarkerSize',6,'CapSize',5);
        end
    end
    xlabel('Independent fitting snapshots'); ylabel(r.metrics{metric});
    xticks([1 2]); xticklabels({'1','8'}); xlim([.85 2.15]); grid on; box on;
    ax.FontSize=13; ax.XLabel.FontSize=22.5; ax.YLabel.FontSize=22.5;
    ax.Title.FontSize=22.5;
    ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.10,0];
    ax.YLabel.Units='normalized'; ax.YLabel.Position=[-.13,.5,0];
    if metric==1
        % Place a compact, editable legend in the gap between the curves.
        for k=1:4
            legend_y=1.32-.17*(k-1);
            plot(ax,[1.03 1.12],[legend_y legend_y], ...
                'Color',curves(k).Color,'LineStyle',curves(k).LineStyle, ...
                'LineWidth',2.3,'HandleVisibility','off');
            plot(ax,1.075,legend_y,'Color',curves(k).Color, ...
                'Marker',curves(k).Marker,'MarkerSize',6,'LineWidth',2.3, ...
                'HandleVisibility','off');
            text(ax,1.16,legend_y,labels{k},'FontSize',22.5, ...
                'VerticalAlignment','middle','Interpreter','none');
        end
    end
end
% Preserve this panel's typography rather than the shared 15-pt defaults.
drawnow; fig.Visible='on'; savefig(fig,fullfile(folder,'parameter_accuracy.fig'));
fig.Visible='off'; exportgraphics(fig,fullfile(folder,'parameter_accuracy.png'),'Resolution',200);
close(fig);
fig=figure('Visible','off','Color','w','Position',[100 100 1200 530]);
colors=lines(3); positions=[.10 .22 .31 .68;.49 .22 .31 .68];
for j=1:2
    ax=axes(fig,'Position',positions(j,:)); hold(ax,'on'); curves=gobjects(1,3);
    for m=1:3
        y=reshape(r.nmse_db(j,:,m),1,[]);
        lo=reshape(r.ci(j,:,1,m),1,[]); hi=reshape(r.ci(j,:,2,m),1,[]);
        curves(m)=errorbar([1 2],y,y-lo,hi-y,'-o','Color',colors(m,:), ...
            'LineWidth',1.6,'MarkerSize',7);
    end
    xlabel('Independent fitting snapshots'); ylabel('Held-out channel NMSE (dB)');
    title(sprintf('Kp = %d',r.grid(j))); xticks([1 2]); xticklabels({'1','8'});
    xlim([.85 2.15]); grid on; box on;
    ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.16,0];
    ax.YLabel.Units='normalized'; ax.YLabel.Position=[-.19,.5,0];
    if j==1
        lg=legend(curves,r.methods,'NumColumns',1,'Visible','on','FontSize',15);
        lg.Units='normalized'; lg.Position=[.815 .42 .175 .26];
    end
end
SaveResponseFigure(fig,folder,'parameter_accuracy_nmse');
end
