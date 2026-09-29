function SaveResponseFigure(fig,folder,name)
%SAVERESPONSEFIGURE Common typography; editable visible FIG and PNG only.
axes_list=findall(fig,'Type','axes');
for a=axes_list.'
    a.FontSize=13; a.XLabel.FontSize=15; a.YLabel.FontSize=15;
end
legends=findall(fig,'Type','legend');
for lg=legends.'
    position=lg.Position;
    set(lg,'FontSize',15,'Box','off','Visible','on');
    drawnow;
    lg.Position=position; % Preserve explicit placement after typography updates.
end
drawnow; fig.Visible='on'; savefig(fig,fullfile(folder,[name,'.fig']));
fig.Visible='off'; exportgraphics(fig,fullfile(folder,[name,'.png']),'Resolution',200);
close(fig);
end
