function VisualiseBlade(Bshow)
% VISUALISEBLADE  Draw the optimised (stage-2) blades from results_peakCp.mat.
%
%   VisualiseBlade      -> detailed views of the B = 3 blade + all 5 rotors
%   VisualiseBlade(4)   -> detailed views of the B = 4 blade instead
%
% Each blade element gets the S822 profile, scaled by its chord and rotated
% by its twist, stacked about the quarter-chord (pitch axis). Coordinates:
%   X = direction of rotation, Y = downwind (axial), Z = radius.
% The chord is at angle twist to the rotor plane, leading edge pointing into
% the relative wind (forward and upwind); suction side faces downwind.

if nargin < 1, Bshow = 3; end
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'results_peakCp.mat')); res = S.res; P = S.P;
res = res(arrayfun(@(s) ~isempty(s.B), res));
F = load(fullfile(here, 'Wind Turbine Blade Optimisation Demo', 'S822_Points.mat'));
xa = F.S822.x(:); ya = F.S822.y(:);               % unit chord, LE at x = 0
figDir = fullfile(here, 'figures');

%% Figure 1: one blade in detail
k = find([res.B] == Bshow, 1);
o = res(k).s2.out; lam = res(k).s2.x(1);
[X, Y, Z] = loft(o.r, o.c, o.twist, xa, ya);
f = figure('Visible', 'off', 'Position', [0 0 1150 820], 'Color', 'w');
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Optimised blade, B = %d   (\\lambda = %.2f, tip speed %.0f m/s, ' ...
    'C_p = %.3f, %.0f W at %.2f m/s)'], Bshow, lam, lam*P.V, res(k).s2.Cp, o.P, P.V));

% (a) 3-D render
ax = nexttile(tl, [2 1]); hold(ax, 'on');
surf(ax, X, Y, Z, 'FaceColor', [0.80 0.84 0.90], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
for j = round(linspace(1, numel(o.r), 7))
    plot3(ax, X(:,j), Y(:,j), Z(:,j), 'k-', 'LineWidth', 0.8);
end
drawHub(ax, P.Rhub, max(o.c));
axis(ax, 'equal'); view(ax, [-35 18]); camlight(ax, 'headlight'); material(ax, 'dull');
xlabel(ax, 'X rotation (m)'); ylabel(ax, 'Y wind (m)'); zlabel(ax, 'radius (m)');
title(ax, '(a) 3-D view'); grid(ax, 'on');

% (b) sections seen looking down the span from the tip
ax = nexttile(tl, [1 2]); hold(ax, 'on'); cm = parula(256);
rr = linspace(P.Rhub, P.R, 9);
for j = 1:numel(rr)
    [~, n] = min(abs(o.r - rr(j)));
    col = cm(round(1 + 255*(o.r(n) - P.Rhub)/(P.R - P.Rhub)), :);
    plot(ax, 1000*X(:,n), 1000*Y(:,n), '-', 'Color', col, 'LineWidth', 1.6);
end
plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
quiver(ax, -150, -55, 0, 25, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 2);
text(ax, -145, -45, 'wind', 'FontSize', 9);
quiver(ax, -150, -55, 30, 0, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 2);
text(ax, -115, -55, 'rotation', 'FontSize', 9);
colormap(ax, cm); clim(ax, [P.Rhub P.R]); cb = colorbar(ax); cb.Label.String = 'radius (m)';
axis(ax, 'equal'); grid(ax, 'on'); xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)');
title(ax, '(b) Sections from root to tip, viewed from the tip (+ = pitch axis, c/4)');

% (c) planform and twist
ax = nexttile(tl, [1 2]); hold(ax, 'on');
le = 0.25*o.c; te = -0.75*o.c;                   % chord drawn flat about c/4
fill(ax, 1000*[o.r fliplr(o.r)], 1000*[le fliplr(te)], [0.80 0.84 0.90], 'EdgeColor', 'k');
plot(ax, 1000*o.r, zeros(size(o.r)), 'k-.');
yyaxis(ax, 'left'); ylabel(ax, 'chord (mm), LE at top');
yyaxis(ax, 'right'); plot(ax, 1000*o.r, o.twist, 'r-', 'LineWidth', 1.6);
ylabel(ax, 'twist \gamma (deg)'); ax.YAxis(2).Color = 'r';
xlabel(ax, 'radius (mm)'); grid(ax, 'on');
idx = [1 round(numel(o.r)/2) numel(o.r)];
ha = {'left', 'left', 'right'};
for j = 1:3
    n = idx(j);
    text(ax, 1000*o.r(n), o.twist(n), sprintf('  c = %.0f mm, \\gamma = %.1f°  ', 1000*o.c(n), o.twist(n)), ...
        'Color', 'r', 'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', ha{j});
end
title(ax, '(c) Planform (chord unrolled flat) and twist distribution');
saveBoth(f, sprintf('t3_blade_geometry_B%d', Bshow), figDir);

%% Figure 2: all rotors, seen from upwind (dotted circle = D 2.2 m)
f = figure('Visible', 'off', 'Position', [0 0 1400 330], 'Color', 'w');
tl = tiledlayout(f, 1, numel(res), 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(res)
    o = res(k).s2.out; B = res(k).B;
    [X, Y, Z] = loft(o.r, o.c, o.twist, xa, ya);
    ax = nexttile(tl, k); hold(ax, 'on');
    for b = 0:B-1
        t = 2*pi*b/B;                                 % rotate blade about the Y axis
        surf(ax, X*cos(t) + Z*sin(t), Y, -X*sin(t) + Z*cos(t), 'FaceColor', [0.80 0.84 0.90], ...
            'EdgeColor', 'none', 'FaceLighting', 'gouraud');
    end
    drawHub(ax, P.Rhub, max(o.c));
    th = linspace(0, 2*pi, 200); plot3(ax, P.R*cos(th), zeros(size(th)), P.R*sin(th), 'k:');
    axis(ax, 'equal'); view(ax, [0 0]); camlight(ax, 'headlight'); material(ax, 'dull');
    axis(ax, 'off'); xlim(ax, [-1.15 1.15]); zlim(ax, [-1.15 1.15]);
    title(ax, {sprintf('B = %d', B), sprintf('\\lambda = %.2f,  C_p = %.3f', res(k).s2.x(1), res(k).s2.Cp)});
end
saveBoth(f, 't3_rotors_all', figDir);
fprintf('Saved t3_blade_geometry_B%d and t3_rotors_all to %s\n', Bshow, figDir);
end

function [X, Y, Z] = loft(r, c, twistDeg, xa, ya)
% Place the S822 profile at every element: scale by chord, shift so the
% pitch axis is the quarter chord, rotate by the twist.
g  = deg2rad(twistDeg(:)');
xc = (xa - 0.25)*c(:)';  yc = ya*c(:)';          % [nPoints x nElements]
X  = -xc.*cos(g) + yc.*sin(g);                   % LE points in +X (rotation)
Y  =  xc.*sin(g) + yc.*cos(g);                   % LE tilts upwind (-Y)
Z  = repmat(r(:)', numel(xa), 1);
end

function drawHub(ax, Rh, w)
[xs, ys, zs] = sphere(30);
surf(ax, Rh*xs, 0.6*w*ys + 0.01, Rh*zs, 'FaceColor', [0.45 0.45 0.5], 'EdgeColor', 'none');
end

function saveBoth(f, n, figDir)
exportgraphics(f, fullfile(figDir, [n '.pdf']));
exportgraphics(f, fullfile(figDir, [n '.png']), 'Resolution', 200);
close(f);
end
