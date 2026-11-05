function E = poly_exponents(d, p)
% All exponent vectors of d variables with total degree <= p.
% d = 3: p=1 -> 4 terms, p=2 -> 10 terms, p=3 -> 20 terms.
[a,b,c] = ndgrid(0:p, 0:p, 0:p);
E = [a(:) b(:) c(:)];
E = E(sum(E,2) <= p, 1:d);
[~, i] = sort(sum(E,2)); E = E(i,:);
end
