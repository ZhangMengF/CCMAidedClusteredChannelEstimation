function plot_CCMVarianceCompensation(folder)
%PLOT_CCMVARIANCECOMPENSATION Left: CE error; right: mean fitted spread.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results','ccm_variance_compensation'); end
load(fullfile(folder,'variance_compensation.mat'),'result'); r=result;
fig=figure('Visible','off','Color','w','Position',[100 100 1450 650]);
positions=[.07 .29 .34 .62;.58 .29 .34 .62];
colors=[.85 .325 .098;.929 .694 .125;.494 .184 .556]; marks={'s','o','d'};
titles={'Delay-prior uncertainty','Angle-prior uncertainty'};
labels={'Mean fitted delay spread (ns)','Mean fitted angular spread (deg)'};
for axis=1:2
    ax=axes(fig,'Position',positions(axis,:)); yyaxis(ax,'left'); hold(ax,'on');
    curves=gobjects(1,4);
    for m=1:3
        curves(m)=plot(r.grid,reshape(r.nmse_db(axis,:,m),1,[]),['-',marks{m}], ...
            'Color',colors(m,:),'LineWidth',1.8,'MarkerSize',7);
    end
    ylabel('NMSE (dB)');
    lo=floor(min(r.nmse_db(:))*2)/2; hi=ceil(max(r.nmse_db(:))*2)/2;
    ylim([lo-.15,hi+.15]); yticks(lo:.5:hi);
    yyaxis(ax,'right');
    mu=reshape(r.spread_mean(axis,:,axis),1,[]);
    ci=reshape(r.spread_ci(axis,:,:,axis),4,2);
    curves(4)=errorbar(r.grid,mu,mu-ci(:,1).',ci(:,2).'-mu,'--x', ...
        'Color',[.1 .45 .65],'LineWidth',1.7,'MarkerSize',8,'CapSize',5);
    ylabel(labels{axis}); ax.YAxis(2).Label.FontSize=15;
    ax.YAxis(1).Color='k'; ax.YAxis(2).Color=[.1 .45 .65];
    title(titles{axis}); xlabel('Actual / assumed error SD, a');
    xticks(r.grid); xlim([.35,4.15]); grid on; box on;
    ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.13,0];
    ax.YAxis(1).Label.Units='normalized'; ax.YAxis(1).Label.Position=[-.105,.5,0];
    ax.YAxis(2).Label.Units='normalized'; ax.YAxis(2).Label.Position=[1.10,.5,0];
    if axis==1
        lg=legend(curves,[r.methods,{'Proposed mean spread (right)'}], ...
            'NumColumns',1,'Box','off');
        lg.Units='normalized'; lg.Position=[.31 .015 .38 .18];
    end
end
SaveResponseFigure(fig,folder,'variance_compensation');
end
