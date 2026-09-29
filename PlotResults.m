function PlotResults(folder, output_folder, names)
%PLOTRESULTS Replot saved experiments without regenerating channels.
if nargin<1, folder=fullfile(fileparts(mfilename('fullpath')),'results'); end
if nargin<2, output_folder=folder; end
if nargin<3, names={'initialization','pilot_only','copula','acceleration'}; end
if ~isfolder(output_folder), mkdir(output_folder); end
for j=1:numel(names)
    file=fullfile(folder,[names{j},'.mat']);
    if ~isfile(file), continue; end
    data=load(file,'result'); r=data.result;
    fig=figure('Visible','off','Color','w','Position',[100,100,1050,430]);
    switch names{j}
        case 'initialization'
            fig.Units='centimeters'; fig.Position=[2,2,8.8,6.5];
            last=min(15,size(r.mean_nll,1)-1); x=0:last; rows=1:last+1;
            nll=r.mean_nll(rows,:); nmse=r.nmse_db(rows,:);
            colors=[0,.447,.741; .85,.325,.098; .49,.18,.56];
            markers={'o','s','^'}; curves=gobjects(1,3);
            yyaxis left; hold on;
            for m=1:3
                curves(m)=plot(x,nll(:,m),'-','Color',colors(m,:),'LineWidth',1.1, ...
                    'Marker',markers{m},'MarkerIndices',1:5:numel(x), ...
                    'MarkerFaceColor','w','MarkerSize',3.5, ...
                    'DisplayName',sprintf('Spread init. %gx',r.scales(m)));
            end
            ylabel('Negative log-likelihood (NLL) / (N Kp)');
            % Both axes use the full plotting height; no artificial bands.
            ylim([.4,.9]);
            yticks(.4:.1:.9); ytickformat('%.1f');
            yyaxis right; hold on;
            for m=1:3
                plot(x,nmse(:,m),'--','Color',colors(m,:),'LineWidth',1.1, ...
                    'Marker',markers{m}, ...
                    'MarkerIndices',1:5:numel(x),'MarkerFaceColor','w','MarkerSize',3.5, ...
                    'HandleVisibility','off');
            end
            ylabel('Channel NMSE (dB)'); xlabel('Outer AO iteration');
            ylim([-9.5,-7.3]);
            yticks(-9.5:.5:-7.5); ytickformat('%.1f');
            styles(1)=plot(nan,nan,'k-','LineWidth',1.1,'DisplayName','NLL / (N Kp) (left axis)');
            styles(2)=plot(nan,nan,'k--','LineWidth',1.1, ...
                'DisplayName','NMSE (right axis)');
            lg=legend([curves,styles],'Location','northeast','Visible','on','Box','off','FontSize',8);
            ax=gca; ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
            ax.FontSize=8; ax.XGrid='on'; ax.YGrid='off'; ax.GridAlpha=.12;
            ax.XLabel.FontSize=9;
            ax.YAxis(1).Label.FontSize=8; ax.YAxis(2).Label.FontSize=9;
            xlim([0,last]); xticks(0:5:last); xtickangle(0); box on;
            % Fixed margins and text positions avoid headless auto-layout overlap.
            ax.Units='normalized'; ax.PositionConstraint='innerposition';
            ax.Position=[.16,.17,.66,.79];
            ax.XLabel.Units='normalized'; ax.XLabel.Position=[.5,-.085,0];
            ax.YAxis(1).Label.Units='normalized';
            ax.YAxis(1).Label.Position=[-.115,.5,0];
            ax.YAxis(2).Label.Units='normalized';
            ax.YAxis(2).Label.Position=[1.15,.5,0];
            lg.ItemTokenSize=[14,8];
            lg.Units='normalized'; lg.Position=[.20,.64,.46,.30];
        case 'pilot_only'
            display_names={'linear','quadratic','FIR','PCHIP','cubic-spline','DFT'};
            positions=[.08 .20 .39 .70; .58 .20 .39 .70];
            for scan=1:2
                ax=axes(fig,'Position',positions(scan,:));
                curves=plot(r.grids{scan},squeeze(r.nmse_db(scan,:,:)),'-o','LineWidth',1.4);
                if scan==1, xlabel('SNR (dB)'); title('Kp = 32');
                else, xlabel('Pilot subcarriers Kp'); title('SNR = 5 dB'); end
                ylabel('NMSE (dB)'); grid on;
                ax.FontSize=20; ax.XLabel.FontSize=16.5; ax.YLabel.FontSize=16.5;
                ax.Title.FontSize=16.5;
                bounds=[floor(min(r.nmse_db(scan,:,:),[],'all')/2)*2, ...
                    ceil(max(r.nmse_db(scan,:,:),[],'all')/2)*2];
                ylim(bounds); yticks(bounds(1):2:bounds(2));
                xlim([min(r.grids{scan}),max(r.grids{scan})]);
                lg=legend(curves,display_names,'FontSize',18,'Location','northeast', ...
                    'Visible','on','Box','off');
                drawnow; lg.Units='normalized';
                lg.Position=[.27+.50*(scan-1),.54,.19,.34];
            end
        case 'copula'
            curves=plot(r.grid,r.nmse_db,'-o','LineWidth',1.7); grid on;
            xlabel('Latent Gaussian correlation r'); ylabel('NMSE (dB)');
            lg=legend(curves,r.methods,'Location','northeast','Visible','on');
            lg.Units='normalized'; lg.Position=[.72,.79,.23,.16];
            title('Delay / absolute-angle coupling');
        case 'acceleration'
            tiledlayout(1,2,'TileSpacing','compact'); nexttile;
            curves=plot(r.grid,r.values(:,2:3),'-o','LineWidth',1.7); grid on;
            xlabel('Kp'); ylabel('NMSE (dB)');
            legend(curves,{'Dense parameters','Accelerated parameters'},'Location','northeast','Visible','on');
            nexttile; plot(r.grid,r.values(:,5),'-o','LineWidth',1.7); grid on;
            xlabel('Kp'); ylabel('Parameter-estimation speedup');
    end
    drawnow;
    fig.Visible='on'; % Saved FIG should open visibly in an interactive session.
    savefig(fig,fullfile(output_folder,[names{j},'.fig']));
    fig.Visible='off';
    resolution=200;
    if strcmp(names{j},'initialization'), resolution=300; end
    exportgraphics(fig,fullfile(output_folder,[names{j},'.png']),'Resolution',resolution);
    close(fig);
end
end
