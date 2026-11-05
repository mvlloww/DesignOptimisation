function y = predict_surrogate(m, U)
% Predict in REAL output units at normalised inputs U.
switch m.spec.type
    case 'poly'
        z = poly_basis(U, m.E) * m.beta;
    case 'rbf'
        R = sqrt(max(sum(U.^2,2) + sum(m.C.^2,2)' - 2*(U*m.C'), 0));
        z = rbf_kernel(R, m.spec.kernel, m.spec.width) * m.w + [ones(size(U,1),1) U] * m.c;
end
y = m.ymu + m.ysd * z;
end
