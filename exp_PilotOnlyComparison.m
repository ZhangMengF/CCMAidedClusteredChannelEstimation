function exp_PilotOnlyComparison(mc_counts,folder)
%EXP_PILOTONLYCOMPARISON Paired, prior-free alternatives for paper baseline.
if nargin<1, mc_counts = [100,50]; end
if nargin<2, folder = fullfile(fileparts(mfilename('fullpath')),'results'); end
[sys,c,~,e,f] = Base(false);
methods = {'linear','quadratic','fir','pchip','spline','dft'};
grids = {[-5,0,5,10,15],[8,16,32,64,128]};
names = {'SNR_dB','Kp'};
result = struct('methods',{methods},'grids',{grids},'mc_counts',mc_counts, ...
    'sys',sys,'channel',c,'estor',e,'seed_rule','100 + scan point');
result.error = cell(2,5); result.energy = cell(2,5);
result.nmse_db = zeros(2,5,numel(methods));
result.values = []; result.columns = [{'scan','value'},methods];
for scan = 1:2
    for j = 1:5
        s = sys; ch = c;
        if scan==1, ch.SNR_dB=grids{scan}(j); else, s.Kp=grids{scan}(j); end
        sc = (1:s.K/s.Kp:s.K).'; rng(100+j,'twister');
        err = zeros(mc_counts(scan),numel(methods)); energy = zeros(mc_counts(scan),1);
        for mc = 1:mc_counts(scan)
            [h,y] = f.DrawTrial(s,ch,f);
            energy(mc) = sum(abs(h).^2);
            for m = 1:numel(methods)
                he = f.PilotOnlyEstimate(y,s.K,s.N,sc,methods{m});
                err(mc,m) = sum(abs(he-h).^2);
            end
        end
        db = 10*log10(sum(err,1)/sum(energy));
        result.error{scan,j}=err; result.energy{scan,j}=energy;
        result.nmse_db(scan,j,:)=db;
        result.values(end+1,:)=[scan,grids{scan}(j),db];
        fprintf('Pilot-only %s=%g: %s dB\n',names{scan},grids{scan}(j),mat2str(db,4));
    end
end
result.score_linear = mean(10.^(result.values(:,3:end)/10),1);
[~,best] = min(result.score_linear); result.selected_method = methods{best};
SaveResult(folder,'pilot_only',result);
fprintf('Best overall (equal weight per scan point): %s\n',result.selected_method);
end
