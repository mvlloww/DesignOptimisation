function [m, best, R] = bb_fit(L, cfg, verbose)
% Fit the best surrogate to everything logged so far (repeats averaged).
if nargin < 3, verbose = false; end
X = L(:,3:5); y = L(:,6);
[Xu, ~, j] = unique(X, 'rows');            % average repeated points
yu = accumarray(j, y, [], @mean);
U = normalise(Xu, cfg);
[best, R] = compare_models(U, yu, cfg, [], [], verbose);
m = fit_surrogate(U, yu, best.spec);
end
