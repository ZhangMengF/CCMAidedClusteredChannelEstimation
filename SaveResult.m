function SaveResult(folder,name,result)
%SAVERESULT Compact raw statistics, environment, and table export.
if ~isfolder(folder), mkdir(folder); end
result.environment = struct('matlab',version,'computer',computer, ...
    'threads',maxNumCompThreads,'hostname',getenv('HOSTNAME'), ...
    'generated_at',char(datetime('now')));
save(fullfile(folder,[name,'.mat']),'result','-v7');
writetable(array2table(result.values,'VariableNames',result.columns), ...
    fullfile(folder,[name,'.csv']));
end
