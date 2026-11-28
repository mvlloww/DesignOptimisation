function VisualiseCapBlades
% VISUALISECAPBLADES  Overlay the variable- and fixed-speed blades optimised
% for winter energy with the rated-power cap (results_winter_cap.mat).
%   (a) 3-D render of both blades on the same pitch axis
%   (b) planform (chord unrolled flat about c/4)
%   (c) twist
%   (d) S822 sections at root, mid-span and tip, viewed from the tip
% Same coordinates as VisualiseBlade: X = rotation, Y = downwind, Z = radius.

here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'results_winter_cap.mat')); c = S.cap; P = S.P;
F = load(fullfile(here, 'Wind Turbine Blade Optimisation Demo', 'S822_Points.mat'));
xa = F.S822.x(:); ya = F.S822.y(:);
r = linspace(P.Rhub, P.R, numel(c.cV));
colV = [0.85 0.33 0.10]; colF = [0.49 0.18 0.56];
k15 = find(c.caps == c.Prated, 1);
labV = sprintf('variable (\\lambda_{opt} = %.2f)', c.xV(1));
labF = sprintf('fixed (tip %.1f m/s)', c.tipF(k15));

f = figure('Visible', 'off', 'Position', [0 0 1150 820], 'Color', 'w');
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Blades optimised for winter energy with a %.1f kW cap: ' ...
    'variable %.0f kWh, fixed %.0f kWh'], c.Prated/1000, c.E(2,k15), c.E(4,k15)));

% (a) 3-D overlay
ax = nexttile(tl, [2 1]); hold(ax, 'on');
[XV, YV, ZV] = loft(r, c.cV, c.twV, xa, ya);
[XF, YF, ZF] = loft(r, c.cF, c.twF, xa, ya);
hF = surf(ax, XF, YF, ZF, 'FaceColor', colF, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
hV = surf(ax, XV, YV, ZV, 'FaceColor', colV, 'FaceAlpha', 0.85, 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
for j = round(linspace(1, numel(r), 6))
    plot3(ax, XF(:,j), YF(:,j), ZF(:,j), '-', 'Color', colF, 'LineWidth', 0.8);
    plot3(ax, XV(:,j), YV(:,j), ZV(:,j), '-', 'Color', colV*0.7, 'LineWidth', 0.8);
end
axis(ax, 'equal', 'tight'); view(ax, [-55 20]); camlight(ax, 'headlight'); material(ax, 'dull');
xlabel(ax, 'X rotation (m)'); ylabel(ax, 'Y wind (m)'); zlabel(ax, 'radius (m)');
legend(ax, [hV hF], {labV, labF}, 'Location', 'southoutside');
title(ax, '(a) 3-D view'); grid(ax, 'on');

% (b) planform
ax = nexttile(tl, 2, [1 2]); hold(ax, 'on');
fill(ax, 1000*[r fliplr(r)], 1000*[0.25*c.cF fliplr(-0.75*c.cF)], colF, 'FaceAlpha', 0.30, 'EdgeColor', colF);
fill(ax, 1000*[r fliplr(r)], 1000*[0.25*c.cV fliplr(-0.75*c.cV)], colV, 'FaceAlpha', 0.55, 'EdgeColor', colV);
plot(ax, 1000*r, zeros(size(r)), 'k-.');
axis(ax, 'equal'); xlim(ax, 1000*[P.Rhub P.R]); ylim(ax, [-200 80]); grid(ax, 'on'); xlabel(ax, 'radius (mm)'); ylabel(ax, 'chord (mm), LE at top');
legend(ax, {labF, labV, 'pitch axis c/4'}, 'Location', 'eastoutside', 'FontSize', 8);
title(ax, sprintf('(b) Planform: root %.0f / %.0f mm, tip %.0f / %.0f mm', ...
    1000*c.cV(1), 1000*c.cF(1), 1000*c.cV(end), 1000*c.cF(end)));

% (c) twist
ax = nexttile(tl, 5); hold(ax, 'on');
plot(ax, 1000*r, c.twV, '-', 'Color', colV, 'LineWidth', 1.8);
plot(ax, 1000*r, c.twF, '-', 'Color', colF, 'LineWidth', 1.8);
grid(ax, 'on'); xlabel(ax, 'radius (mm)'); ylabel(ax, 'twist \gamma (deg)');
legend(ax, {labV, labF}, 'FontSize', 8);
title(ax, '(c) Twist');

% (d) sections from the tip
ax = nexttile(tl, 6); hold(ax, 'on');
rr = [P.Rhub, (P.Rhub + P.R)/2, P.R]; ls = {'-', '--', ':'}; h = gobjects(1, 6);
for j = 1:3
    [~, n] = min(abs(r - rr(j)));
    h(j)   = plot(ax, 1000*XV(:,n), 1000*YV(:,n), ls{j}, 'Color', colV, 'LineWidth', 1.6);
    h(j+3) = plot(ax, 1000*XF(:,n), 1000*YF(:,n), ls{j}, 'Color', colF, 'LineWidth', 1.6);
end
plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
axis(ax, 'equal'); grid(ax, 'on'); xlabel(ax, 'X rotation (mm)'); ylabel(ax, 'Y wind (mm)');
legend(ax, h, {sprintf('variable, r = %.2f m', rr(1)), sprintf('variable, r = %.2f m', rr(2)), ...
    sprintf('variable, r = %.2f m', rr(3)), sprintf('fixed, r = %.2f m', rr(1)), ...
    sprintf('fixed, r = %.2f m', rr(2)), sprintf('fixed, r = %.2f m', rr(3))}, ...
    'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 7);
title(ax, {'(d) Sections viewed from the tip', '(+ = pitch axis)'});

figDir = fullfile(here, 'figures');
exportgraphics(f, fullfile(figDir, 't3_cap_blades_overlay.pdf'));
exportgraphics(f, fullfile(figDir, 't3_cap_blades_overlay.png'), 'Resolution', 200);
close(f);
fprintf('Saved t3_cap_blades_overlay to %s\n', figDir);
end

function [X, Y, Z] = loft(r, c, twistDeg, xa, ya)
% S822 profile at every element, scaled by chord, about c/4, rotated by twist
g  = deg2rad(twistDeg(:)');
xc = (xa - 0.25)*c(:)';  yc = ya*c(:)';
X  = -xc.*cos(g) + yc.*sin(g);
Y  =  xc.*sin(g) + yc.*cos(g);
Z  = repmat(r(:)', numel(xa), 1);
end
