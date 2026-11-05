function x = bb_next(mode)
% BB_NEXT  Print the next point to dial in.
%   bb_next            fixed run sheet first, then auto (alternate exploit/explore)
%   bb_next('exploit') surrogate optimum
%   bb_next('explore') most isolated point among the top-25% predicted region
% Stores the point (and its prediction) so bb_record(y) only needs the reading.
if nargin < 1, mode = 'auto'; end
cfg = bb_config();
[L, plan] = bb_load(cfg);
nDone = size(L,1);
if nDone >= cfg.budget, error('Budget of %d calls used.', cfg.budget); end

ypred = NaN;
if nDone < size(plan,1)
    x = plan(nDone+1, 3:5); phase = plan(nDone+1, 2);
    if phase == 2 || phase == 3                     % predict-then-query on hold-out/repeats too
        tr = L(L(:,2)==1,:);
        if size(tr,1) >= 10, ypred = predict_surrogate(bb_fit(tr, cfg), normalise(x, cfg)); end
    end
else
    [m, best] = bb_fit(L, cfg);
    nInf = sum(L(:,2) == 4);
    if strcmp(mode, 'auto')
        if mod(nInf, 2) == 0, mode = 'exploit'; else, mode = 'explore'; end
    end
    Udata = normalise(L(:,3:5), cfg);
    if any(strcmp(mode, {'exploit','confirm'}))
        [u, ypred] = maximise_surrogate(m);
        dnear = min(sqrt(sum((Udata - u).^2, 2)));
        if dnear < 0.02
            fprintf('Surrogate optimum already sampled (within %.3f normalised): converged.\n', dnear);
            if strcmp(mode, 'exploit')
                fprintf('-> switching to explore so the call is not wasted.\n');
                mode = 'explore';
            else
                % Confirm = local robustness check: step 5%% of range along one
                % axis from the optimum, cycling x1, x2, x3. Gives measured
                % local sensitivities for the report.
                nC = sum(L(:,2) == 5); k = mod(nC, 3) + 1;
                step = zeros(1,3); step(k) = 0.05;
                if u(k) + step(k) > 1, step = -step; end
                u = u + step; ypred = predict_surrogate(m, u);
                fprintf('-> local check: +/-5%% step in %s\n', cfg.names{k});
            end
        end
    end
    switch mode
        case {'exploit','confirm'}
            % u, ypred already set above
        case 'explore'
            Cc = rand(20000, 3);
            yc = predict_surrogate(m, Cc);
            keep = yc >= quantile_(yc, 0.75);
            Cc = Cc(keep,:); yc = yc(keep);
            dmin = min(max(sum(Cc.^2,2) + sum(Udata.^2,2)' - 2*(Cc*Udata'), 0), [], 2);
            [~, k] = max(dmin); u = Cc(k,:); ypred = yc(k);
        otherwise, error('mode must be auto | exploit | explore | confirm');
    end
    x = denormalise(u, cfg, true);
    if strcmp(mode,'confirm'), phase = 5; else, phase = 4; end
    fprintf('Model: %s\n', best.name);
end

save('bb_pending.mat', 'x', 'phase', 'ypred');
fprintf('\n>>> CALL %d/%d  [%s]\n', nDone+1, cfg.budget, phase_name(phase));
fprintf('    x1 = %g   x2 = %g   x3 = %g\n', x);
if ~isnan(ypred), fprintf('    predicted MCSI = %.4g\n', ypred); end
fprintf('    then: bb_record(<MCSI reading>)\n\n');
end

function q = quantile_(v, p)
v = sort(v(:)); q = v(max(1, round(p*numel(v))));
end

function s = phase_name(p)
names = {'LHS training','hold-out test','repeat','infill','confirm'};
s = names{p};
end
