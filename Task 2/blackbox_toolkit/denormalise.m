function X = denormalise(U, cfg, roundToDial)
% Unit cube -> real units, optionally snapped to dial resolution
if nargin < 3, roundToDial = true; end
X = cfg.lb + U .* (cfg.ub - cfg.lb);
if roundToDial
    X = round(X ./ cfg.res) .* cfg.res;
    X = min(max(X, cfg.lb), cfg.ub);
end
end
