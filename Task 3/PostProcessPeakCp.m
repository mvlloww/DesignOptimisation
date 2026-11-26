function res = PostProcessPeakCp(res, P, L, figDir)
% POSTPROCESSPEAKCP  Check the optimised blades and draw the figures.
% Runs at the end of OptimisePeakCp, or on its own (loads results_peakCp.mat).
%
% Each finished blade is treated as BUILT: chord and twist are now fixed, so
% the angle of attack follows the flow (BEMT alphaMode 'twist').
%   1. Re-evaluate at the design lambda -> must reproduce the optimised Cp
%   2. Same blade without tip loss       -> size of the tip-loss correction
%   3. Same blade with 120 / 240 elements -> discretisation error
%   4. Cp-lambda sweep 2..12              -> how peaky the design is
% In fixed-twist mode the BEMT iteration can oscillate (AoA depends on phi),
% so these runs use under-relaxation 0.3. That only changes the path to the
% answer, not the answer: check 1 reproduces the optimised Cp.

here = fileparts(mfilename('fullpath'));
if nargin == 0
    S = load(fullfile(here, 'results_peakCp.mat')); res = S.res; P = S.P;
    L = S822Lookup(fullfile(here, 'Wind Turbine Blade Optimisation Demo'));
    figDir = fullfile(here, 'figures');
end
res = res(arrayfun(@(s) ~isempty(s.B), res));
tw  = @(tip) struct('lookup', L, 'tipLoss', tip, 'alphaMode', 'twist', 'relax', 0.3);
run = @(B, lam, c, t, N, tip) BEMT(B, [P.Rhub P.R], lam, c, t, N, P.V, P.nu, P.tol, 5000, tw(tip));

fprintf('\n B  lambda  tip(m/s)  Cp_lin  Cp_spline  Schmitz  recheck  noTipLoss  N=120   N=240   sweepPeak@lam   P(W)  CT\n');
for k = 1:numel(res)
    B = res(k).B; o2 = res(k).s2.out; lam = res(k).s2.x(1); r = o2.r;
    o = run(B, lam, o2.c, o2.twist, P.N, true);  ok = all(o.converged);
    res(k).recheckCp = o.Cp;
    o = run(B, lam, o2.c, o2.twist, P.N, false); ok = ok && all(o.converged);
    res(k).noTipCp = o.Cp;
    Ns = [120 240]; res(k).cpN = zeros(1,2);
    for j = 1:2
        rN = linspace(P.Rhub, P.R, Ns(j));
        o = run(B, lam, interp1(r, o2.c, rN), interp1(r, o2.twist, rN), Ns(j), true);
        ok = ok && all(o.converged); res(k).cpN(j) = o.Cp;
    end
    lamS = linspace(2, 12, 81); cpS = nan(size(lamS));
    for j = 1:numel(lamS)
        o = run(B, lamS(j), o2.c, o2.twist, P.N, true);
        if all(o.converged), cpS(j) = o.Cp; end
    end
    res(k).lamSweep = lamS; res(k).cpSweep = cpS; res(k).checksConverged = ok;
    [m, j] = max(cpS);
    fprintf('%2d  %5.2f   %5.1f    %.4f  %.4f    %.4f   %.4f   %.4f     %.4f  %.4f  %.4f@%.2f  %4.0f  %.2f  %s\n', ...
        B, lam, lam*P.V, res(k).s1.Cp, res(k).s2.Cp, res(k).schmitz.Cp, res(k).recheckCp, ...
        res(k).noTipCp, res(k).cpN, m, lamS(j), o2.P, o2.CT, ...
        sprintf('converged=%d sweepNaN=%d/81', ok, sum(isnan(cpS))));
end
if nargin == 0, save(fullfile(here, 'results_peakCp.mat'), 'res', 'P'); end
makeFigures(res, P, L, figDir);
end

function makeFigures(res, P, L, figDir)
cols = lines(numel(res));
names = arrayfun(@(s) sprintf('B = %d', s.B), res, 'UniformOutput', false);
save2 = @(f, n) saveBoth(f, n, figDir);
% 1. best-L/D angle vs Re
f = figure('Visible', 'off', 'Position', [0 0 560 380]);
semilogx(L.bestTable.Re, L.bestTable.alpha, 'k.', 'MarkerSize', 12); hold on;
semilogx(L.bestTable.Re, L.alphaBest(L.bestTable.Re), 'r-', 'LineWidth', 1.5);
xlabel('Reynolds number'); ylabel('\alpha at best L/D (deg)'); grid on;
legend('S822 data (search)', 'smooth fit used by optimiser', 'Location', 'northeast');
save2(f, 't3_alpha_best');
% 2. Cp vs B
f = figure('Visible', 'off', 'Position', [0 0 620 400]);
Y = [arrayfun(@(s) s.schmitz.Cp, res); arrayfun(@(s) s.s1.Cp, res); arrayfun(@(s) s.s2.Cp, res)]';
bar([res.B], Y); hold on; yline(16/27, 'k--', 'Betz limit');
xlabel('Number of blades B'); ylabel('Peak C_p (tip loss on)'); ylim([0 0.65]); grid on;
legend('Schmitz ideal (same \lambda)', 'Stage 1: linear chord', 'Stage 2: spline chord', 'Location', 'southeast');
save2(f, 't3_cp_vs_B');
% 3. chord and twist of the stage-2 blades
f = figure('Visible', 'off', 'Position', [0 0 620 560]);
subplot(2,1,1); hold on; grid on;
for k = 1:numel(res), plot(res(k).s2.out.r, 1000*res(k).s2.out.c, 'Color', cols(k,:), 'LineWidth', 1.5); end
ylabel('Chord (mm)'); legend(names);
subplot(2,1,2); hold on; grid on;
for k = 1:numel(res), plot(res(k).s2.out.r, res(k).s2.out.twist, 'Color', cols(k,:), 'LineWidth', 1.5); end
xlabel('Radius r (m)'); ylabel('Twist \gamma (deg)');
save2(f, 't3_best_blades');
% 4. Cp-lambda of each finished blade
f = figure('Visible', 'off', 'Position', [0 0 620 400]); hold on; grid on;
for k = 1:numel(res)
    plot(res(k).lamSweep, res(k).cpSweep, 'Color', cols(k,:), 'LineWidth', 1.5);
    plot(res(k).s2.x(1), res(k).s2.Cp, 'o', 'Color', cols(k,:), 'MarkerFaceColor', cols(k,:), 'HandleVisibility', 'off');
end
xline(P.lamMax, 'k--', 'tip-speed limit'); xlabel('Tip-speed ratio \lambda'); ylabel('C_p (fixed blade shape)');
ylim([0 0.5]); legend(names, 'Location', 'southeast');
save2(f, 't3_cp_lambda');
end

function saveBoth(f, n, figDir)
exportgraphics(f, fullfile(figDir, [n '.pdf']));
exportgraphics(f, fullfile(figDir, [n '.png']), 'Resolution', 200);
close(f);
end
