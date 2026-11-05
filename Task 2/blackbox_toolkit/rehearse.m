% REHEARSE  Dry-run the whole session against fake_blackbox.
% Do this at least once on the machine you will take into the room.
clear; clc;
cfg = bb_config();
if exist(cfg.logFile,'file'), delete(cfg.logFile); end
if exist('bb_pending.mat','file'), delete('bb_pending.mat'); end
make_sampling_plan(cfg, false);
plan = dlmread(cfg.planFile, ',', 1, 0);
tic;
for k = 1:size(plan,1)
    x = bb_next(); bb_record(fake_blackbox(x));
end
bb_status(false);                                    % decision point after fixed block
for k = 1:cfg.nInfill
    x = bb_next('auto'); bb_record(fake_blackbox(x));
end
for k = 1:cfg.nConfirm
    x = bb_next('confirm'); bb_record(fake_blackbox(x));
end
bb_status(true);
fprintf('Compute time for whole session: %.1f s\n', toc);

% Compare with the true optimum of the fake function
obj = @(u) -fake_blackbox(denormalise(min(max(u,0),1), cfg, false));
best = inf;
for s = 1:20
    [u, f] = fminsearch(obj, rand(1,3), optimset('Display','off'));
    if f < best, best = f; ut = min(max(u,0),1); end
end
L = dlmread(cfg.logFile, ',', 1, 0);
[ym, i] = max(L(:,6));
fprintf('\nTRUE optimum  : x = [%.4g %.4g %.4g], MCSI = %.5g\n', denormalise(ut,cfg,false), -best);
fprintf('Best observed : x = [%.4g %.4g %.4g], MCSI = %.5g  (%.3f%% below true)\n', ...
        L(i,3:5), ym, 100*(-best - ym)/(-best));
