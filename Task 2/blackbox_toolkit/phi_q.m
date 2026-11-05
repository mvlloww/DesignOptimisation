function phi = phi_q(U, q)
% Morris-Mitchell space-filling criterion (lower = better spread).
% phi_q = (sum d_ij^-q)^(1/q); large q -> pure maximin.
if nargin < 2, q = 50; end
D = pdist_sq(U);
d = sqrt(D(triu(true(size(D)), 1)));
dmin = min(d);
% scale by dmin for numerical stability at large q
phi = (sum((d / dmin).^(-q)))^(1/q) / dmin;
end

function D = pdist_sq(U)
s = sum(U.^2, 2);
D = max(s + s' - 2 * (U * U'), 0);
end
