function U = normalise(X, cfg)
% Real units -> unit cube [0,1]^3
U = (X - cfg.lb) ./ (cfg.ub - cfg.lb);
end
