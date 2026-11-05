function [nrmse, q2, yhat] = loocv(U, y, spec)
% Leave-one-out cross-validation. Brute force: n <= 50 so it is instant.
% nrmse = RMSE / std(y)   (0 = perfect, 1 = no better than predicting the mean)
% q2    = 1 - PRESS / SS_tot  (predictive R^2)
n = numel(y); yhat = nan(n,1);
for i = 1:n
    k = true(n,1); k(i) = false;
    m = fit_surrogate(U(k,:), y(k), spec);
    if ~m.ok, nrmse = inf; q2 = -inf; return; end
    yhat(i) = predict_surrogate(m, U(i,:));
end
e = y - yhat;
nrmse = sqrt(mean(e.^2)) / std(y);
q2 = 1 - sum(e.^2) / sum((y - mean(y)).^2);
end
