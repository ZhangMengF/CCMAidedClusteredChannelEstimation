function plot_CombinedResults(folder)
%PLOT_COMBINEDRESULTS Replace data only in the approved, frozen Fig.1 layout.
root=fileparts(mfilename('fullpath'));
if nargin<1, folder=fullfile(root,'results','combined_accelerated'); end
fig=openfig(fullfile(root,'figure_assets','Fig1_layout.fig'),'invisible');
fig.Units='pixels'; fig.Position=[100 100 1400 450];
axes_list=findall(fig,'Type','axes');
% Approved exception: expose the new -5.02 dB point without changing ticks.
ax_c=axes_list(arrayfun(@(a) startsWith(string(a.Title.String),'(c)'),axes_list));
ax_c.YLim=[-10 -5];
layout=cell(numel(axes_list),1);
for j=1:numel(axes_list)
    layout{j}=get(axes_list(j),{'Position','XLim','YLim','XTick','YTick','FontSize'});
end
sources={'scan_SNR_dB','scan_Kp'};
display_names={'CPM-aided','Fixed-profile','Proposed','Oracle','Pilot-only'};
columns={'CPM','Fixed_profile','Proposed','Oracle','Pilot_only'};
for panel=1:2
    data=load(fullfile(folder,[sources{panel},'.mat']),'result'); r=data.result;
    assert(strcmp(r.estor.parameter_backend,'beam_delay_diag'));
    assert(strcmp(r.estor.final_solver,'fft_pcg'));
    ax=axes_list(arrayfun(@(a) startsWith(string(a.Title.String),sprintf('(%c)',char('a'+panel-1))),axes_list));
    for k=1:5
        h=findobj(ax,'Type','line','DisplayName',display_names{k});
        assert(isequal(h.XData(:),r.values(:,1)),'Do not change the approved sample grid.');
        h.YData=r.values(:,strcmp(r.columns,columns{k})).';
        assert(all(h.YData>=ax.YLim(1)&h.YData<=ax.YLim(2)),'New curve exceeds the fixed plotting range.');
    end
end
data=load(fullfile(folder,'history_BlockProb.mat'),'result'); r=data.result;
assert(strcmp(r.estor.parameter_backend,'beam_delay_diag'));
assert(strcmp(r.estor.final_solver,'fft_pcg'));
ax=axes_list(arrayfun(@(a) startsWith(string(a.Title.String),'(c)'),axes_list));
curves=findall(ax,'Type','line');
offline=curves(arrayfun(@(h) numel(h.XData)>2,curves));
% The frozen layout orders the old offline endpoints as blockage 0, .3, .5.
[~,order]=sort(arrayfun(@(h) h.YData(end),offline)); offline=offline(order);
probabilities=[0 .3 .5];
for j=1:3
    h=offline(j); column=find(abs(r.scan_values-probabilities(j))<1e-12)+1;
    online=curves(arrayfun(@(q) numel(q.XData)==2 && isequal(q.Color,h.Color),curves));
    assert(isequal(h.XData(:),r.values(2:end,1)));
    h.YData=r.values(2:end,column).';
    online.YData=repmat(r.values(1,column),size(online.YData));
    assert(all([h.YData online.YData]>=ax.YLim(1)&[h.YData online.YData]<=ax.YLim(2)));
end
drawnow;
for j=1:numel(axes_list)
    assert(isequal(layout{j},get(axes_list(j),{'Position','XLim','YLim','XTick','YTick','FontSize'})), ...
        'The approved figure layout must remain unchanged.');
end
% Export before making visible: a headless 1024-pixel screen otherwise resizes it.
exportgraphics(fig,fullfile(folder,'Result.png'),'Resolution',600);
fig.Visible='on'; fig.Position=[100 100 1400 450]; drawnow;
savefig(fig,fullfile(folder,'Combined_Plots.fig'));
close(fig);
end
