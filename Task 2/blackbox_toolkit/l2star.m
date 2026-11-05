function D = l2star(U)
% L2-star discrepancy (Warnock's closed form). Fast; used to screen
% thousands of candidate plans. Points must be in [0,1]^d.
[n, d] = size(U);
t1 = 3^(-d);
t2 = (2^(1 - d) / n) * sum(prod(1 - U.^2, 2));
P = ones(n, n);
for k = 1:d
    P = P .* (1 - max(U(:, k), U(:, k)'));
end
t3 = sum(P(:)) / n^2;
D = sqrt(max(t1 - t2 + t3, 0));
end
