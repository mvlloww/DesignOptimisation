function m = fit_surrogate(U, y, spec)
% FIT_SURROGATE  Fit polynomial or RBF in NORMALISED inputs U in [0,1]^3.
% y is standardised internally (z-score) for conditioning.
%
% spec.type   = 'poly' | 'rbf'
% spec.order  = polynomial total degree (poly)
% spec.kernel = 'gaussian'|'multiquadric'|'invmultiquadric'|'cubic'|'thinplate'
% spec.width  = shape parameter (gaussian/MQ/IMQ), normalised space
% spec.lambda = nugget (regularisation) for RBF
m.spec = spec;
m.ymu = mean(y); m.ysd = std(y); if m.ysd == 0, m.ysd = 1; end
z = (y - m.ymu) / m.ysd;
m.ok = true;

switch spec.type
    case 'poly'
        m.E = poly_exponents(size(U,2), spec.order);
        P = poly_basis(U, m.E);
        if size(P,2) >= size(P,1)
            m.ok = false;  % under-determined: not enough points for this order
            m.beta = zeros(size(P,2),1);
        else
            m.beta = P \ z;
        end
    case 'rbf'
        m.C = U;
        n = size(U,1);
        A = rbf_kernel(sqrt(sqd(U,U)), spec.kernel, spec.width) + spec.lambda*eye(n);
        P = [ones(n,1) U];                     % linear tail: needed for cubic/TPS,
        K = [A P; P' zeros(size(P,2))];        % and helps the others extrapolate
        if rcond(K) < 1e-14, m.ok = false; end
        sol = K \ [z; zeros(size(P,2),1)];
        m.w = sol(1:n); m.c = sol(n+1:end);
    otherwise
        error('unknown surrogate type');
end
end

function D = sqd(A,B)
D = max(sum(A.^2,2) + sum(B.^2,2)' - 2*(A*B'), 0);
end
