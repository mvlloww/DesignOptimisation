function P = poly_basis(U, E)
% Centre inputs on 0.5 so the basis is better conditioned.
V = U - 0.5;
P = ones(size(U,1), size(E,1));
for j = 1:size(E,1)
    for k = 1:size(E,2)
        if E(j,k) > 0, P(:,j) = P(:,j) .* V(:,k).^E(j,k); end
    end
end
end
