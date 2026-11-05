function bb_plots(L, m, best, u, cfg)
% Report figures: LOOCV parity plot + 2D slices through the optimum.
X = L(:,3:5); [Xu,~,j] = unique(X,'rows'); yu = accumarray(j, L(:,6), [], @mean);

f1 = figure('Name','LOOCV parity','Color','w'); hold on; box on; grid on;
plot(yu, best.yloo, 'ko', 'MarkerFaceColor', [0.2 0.4 0.8]);
lim = [min([yu; best.yloo]) max([yu; best.yloo])];
plot(lim, lim, 'k--');
xlabel('observed MCSI'); ylabel('LOOCV predicted MCSI');
title(sprintf('%s: NRMSE %.3f, Q^2 %.3f', best.name, best.nrmse, best.q2));
try print(f1, 'fig_loocv_parity.png', '-dpng', '-r200'); catch, end

f2 = figure('Name','Surrogate slices','Color','w');
pairs = [1 2 3; 1 3 2; 2 3 1];
g = linspace(0,1,60);
xs = denormalise(u, cfg, false);
for p = 1:3
    a = pairs(p,1); b = pairs(p,2); c = pairs(p,3);
    [A,B] = meshgrid(g,g);
    Q = zeros(numel(A),3); Q(:,a) = A(:); Q(:,b) = B(:); Q(:,c) = u(c);
    Z = reshape(predict_surrogate(m, Q), size(A));
    subplot(1,3,p); hold on; box on;
    contourf(cfg.lb(a)+A*(cfg.ub(a)-cfg.lb(a)), cfg.lb(b)+B*(cfg.ub(b)-cfg.lb(b)), Z, 20, 'LineColor','none');
    colorbar;
    plot(Xu(:,a), Xu(:,b), 'w.', 'MarkerSize', 10);
    plot(xs(a), xs(b), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
    xlabel(cfg.names{a}); ylabel(cfg.names{b});
    title(sprintf('slice at %s = %.3g', cfg.names{c}, xs(c)));
end
try print(f2, 'fig_surrogate_slices.png', '-dpng', '-r200'); catch, end
end
