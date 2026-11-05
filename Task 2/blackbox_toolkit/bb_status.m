function bb_status(makePlots)
% BB_STATUS  Decision point. Run after the fixed block (call 32) and
% whenever you want to re-assess during infill.
if nargin < 1, makePlots = false; end
cfg = bb_config();
L = bb_load(cfg);
fprintf('\n===== %d / %d calls used =====\n', size(L,1), cfg.budget);
[ymax, i] = max(L(:,6));
fprintf('Best observed: MCSI = %.6g at x = [%g %g %g] (run %d)\n', ymax, L(i,3:5), L(i,1));

tr = L(L(:,2)==1,:); te = L(L(:,2)==2,:);
if ~isempty(te)
    fprintf('\n--- Train on LHS (n=%d), test on hold-out (n=%d) ---', size(tr,1), size(te,1));
    compare_models(normalise(tr(:,3:5),cfg), tr(:,6), cfg, normalise(te(:,3:5),cfg), te(:,6), true);
end

fprintf('\n--- All data (repeats averaged) ---');
[m, best, R] = bb_fit(L, cfg, true);
[u, yp] = maximise_surrogate(m);
xs = denormalise(u, cfg, true);
fprintf('\nSurrogate optimum: x = [%g %g %g], predicted MCSI = %.6g\n', xs, yp);
onB = (u < 0.01) | (u > 0.99);
if any(onB), fprintf('!! Optimum on bound for: %s\n', strjoin(cfg.names(onB), ', ')); end

inf_ = L(L(:,2)>=4 & ~isnan(L(:,7)),:);
if ~isempty(inf_)
    e = 100*(inf_(:,6) - inf_(:,7)) / (max(L(:,6)) - min(L(:,6)));
    fprintf('Predict-then-query errors on infill (%% of observed MCSI range):\n  %s\n', ...
            sprintf('%+.1f ', e));
end

if makePlots, bb_plots(L, m, best, u, cfg); end
end
