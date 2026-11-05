function [L, plan] = bb_load(cfg)
% Returns log matrix L = [run phase x1 x2 x3 y ypred] and the run sheet.
if exist(cfg.logFile, 'file')
    L = dlmread(cfg.logFile, ',', 1, 0);
    if isempty(L), L = zeros(0,7); end
else
    L = zeros(0,7);
end
if exist(cfg.planFile, 'file')
    plan = dlmread(cfg.planFile, ',', 1, 0);
else
    error('No run sheet: run make_sampling_plan first');
end
end
