function [best, R] = compare_models(U, y, cfg, Ute, yte, verbose)
% COMPARE_MODELS  Poly (orders in cfg.polyOrders) vs RBF (each kernel,
% width tuned by LOOCV). Ranks by LOOCV; if a hold-out set is given it
% also reports hold-out error, which is the UNBIASED check because LOOCV
% was used to pick the RBF width (selection bias).
if nargin < 4, Ute = []; yte = []; end
if nargin < 6, verbose = true; end
R = struct('name',{},'spec',{},'nrmse',{},'q2',{},'test',{},'yloo',{});

% --- polynomials
for p = cfg.polyOrders
    s = struct('type','poly','order',p,'kernel','','width',NaN,'lambda',0);
    [e, q, yl] = loocv(U, y, s);
    R(end+1) = pack(sprintf('poly order %d (%d terms)', p, size(poly_exponents(3,p),1)), s, e, q, yl); %#ok<AGROW>
end

% --- RBFs
for k = 1:numel(cfg.rbfKernels)
    ker = cfg.rbfKernels{k};
    if any(strcmp(ker, {'cubic','thinplate'})), widths = 1; else, widths = cfg.rbfWidths; end
    bestE = inf;
    for w = widths
        for lam = cfg.rbfLambda
            s = struct('type','rbf','order',NaN,'kernel',ker,'width',w,'lambda',lam);
            [e, q, yl] = loocv(U, y, s);
            if e < bestE, bestE = e; bs = s; bq = q; byl = yl; end
        end
    end
    if isinf(bestE), continue; end
    if numel(widths) > 1, nm = sprintf('RBF %s (s=%.3f)', ker, bs.width);
    else, nm = sprintf('RBF %s', ker); end
    R(end+1) = pack(nm, bs, bestE, bq, byl); %#ok<AGROW>
end

% --- hold-out
for i = 1:numel(R)
    R(i).test = NaN;
    if ~isempty(yte) && isfinite(R(i).nrmse)
        m = fit_surrogate(U, y, R(i).spec);
        R(i).test = sqrt(mean((yte - predict_surrogate(m, Ute)).^2)) / std(y);
    end
end

[~, o] = sort([R.nrmse]); R = R(o);
best = R(1);

if verbose
    fprintf('\n%-34s %10s %8s %12s\n', 'model', 'LOO NRMSE', 'LOO Q2', 'test NRMSE');
    fprintf('%s\n', repmat('-', 1, 68));
    for i = 1:numel(R)
        fprintf('%-34s %10.3f %8.3f %12.3f\n', R(i).name, R(i).nrmse, R(i).q2, R(i).test);
    end
    fprintf('Selected (lowest LOOCV): %s\n', best.name);
    if ~isempty(yte)
        [~, it] = min([R.test]);
        if it ~= 1
            fprintf('NOTE: hold-out prefers %s -> LOOCV ranking may be selection-biased; check gap.\n', R(it).name);
        end
    end
end
end

function r = pack(name, spec, e, q, yl)
r = struct('name',name,'spec',spec,'nrmse',e,'q2',q,'test',NaN,'yloo',yl);
end
