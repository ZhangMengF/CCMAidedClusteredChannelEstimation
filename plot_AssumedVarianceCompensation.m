function plot_AssumedVarianceCompensation(folder)
%PLOT_ASSUMEDVARIANCECOMPENSATION Fixed geometry, varying assumed uncertainty.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','assumed_variance_compensation'); end
load(fullfile(folder,'assumed_variance_compensation.mat'),'result'); r=result;
fig=figure('Visible','off','Color','w','Position',[100 100 1450 800]);
positions=[.07 .16 .34 .75;.58 .16 .34 .75];
colors=[.85 .325 .098;.929 .694 .125;.494 .184 .556]; marks={'s','o','d'};
titles={'Assumed delay uncertainty','Assumed angle uncertainty'};
labels={'Mean fitted delay spread (ns)','Mean fitted angular spread (deg)'};
for axis=1:2
    ax=axes(fig,'Position',positions(axis,:)); yyaxis(ax,'left'); hold(ax,'on'); curves=gobjects(1,4);
    for m=1:3
        curves(m)=plot(r.grid,reshape(r.nmse_db(axis,:,m),1,[]),['-',marks{m}], ...
            'Color',colors(m,:),'LineWidth',2.5,'MarkerSize',7);
    end
    ylabel('NMSE (dB)');
    lo=floor(min(r.nmse_db(:))); hi=ceil(max(r.nmse_db(:)));
    legend_top=min(r.nmse_db(:))-.25;
    ylim([min(lo-.8,legend_top-.9),hi+.15]); yticks(lo:1:hi);
    yyaxis(ax,'right'); mu=reshape(r.spread_mean(axis,:,axis),1,[]);
    ci=reshape(r.spread_ci(axis,:,:,axis),numel(r.grid),2);
    curves(4)=errorbar(r.grid,mu,mu-ci(:,1).',ci(:,2).'-mu,'--x', ...
        'Color',[.1 .45 .65],'LineWidth',2.5,'MarkerSize',8,'CapSize',5);
    if r.channel.SNR_dB==15 && r.sys.Kp==64
        if axis==1, ylim([7.5 10]); else, ylim([6 10]); end
    else
        spread_range=max(ci(:,2))-min(ci(:,1));
        ylim([min(ci(:,1))-.6*spread_range,max(ci(:,2))+.05*spread_range]);
    end
    ylabel(labels{axis}); ax.FontSize=13;
    ax.YAxis(1).Label.FontSize=22.5; ax.YAxis(2).Label.FontSize=22.5;
    ax.YAxis(1).Color='k'; ax.YAxis(2).Color=[.1 .45 .65];
    title(titles{axis}); xlabel('Actual / assumed error SD, a');
    ax.Title.FontSize=21.45; ax.XLabel.FontSize=22.5;
    ax.XScale='log'; xticks(r.grid); xticklabels(string(r.grid)); xlim([.23,4.3]); grid on; box on;
    ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.075,0];
    ax.YAxis(1).Label.Units='normalized'; ax.YAxis(1).Label.Position=[-.09,.5,0];
    ax.YAxis(2).Label.Units='normalized'; ax.YAxis(2).Label.Position=[1.09,.5,0];
    if axis==2
        % Editable legend below Oracle, clear of both y-axis curve groups.
        yyaxis(ax,'left');
        legend_labels=[r.methods,{'Mean refitted spread (right)'}];
        for k=1:4
            legend_y=legend_top-.22*(k-1);
            plot(ax,[.30 .36],[legend_y legend_y],'Color',curves(k).Color, ...
                'LineStyle',curves(k).LineStyle,'Marker','none', ...
                'LineWidth',2.5,'HandleVisibility','off');
            plot(ax,sqrt(.30*.36),legend_y,'Color',curves(k).Color, ...
                'Marker',curves(k).Marker,'MarkerSize',7,'LineWidth',2.5, ...
                'HandleVisibility','off');
            text(ax,.39,legend_y,legend_labels{k},'FontSize',22.5, ...
                'VerticalAlignment','middle','Interpreter','none');
        end
    end
end
drawnow; fig.Visible='on'; savefig(fig,fullfile(folder,'assumed_variance_compensation.fig'));
fig.Visible='off'; exportgraphics(fig,fullfile(folder,'assumed_variance_compensation.png'),'Resolution',200);
close(fig);
end
