function U = lhs_random(n, d)
% Random Latin hypercube in [0,1]^d, points at cell centres.
% Toolbox-free (no lhsdesign needed).
U = zeros(n, d);
for k = 1:d
    U(:, k) = (randperm(n)' - 0.5) / n;
end
end
