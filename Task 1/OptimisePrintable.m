clc, clear all, close all

%Lightest PRINTABLE designs, chosen directly in whole 0.4 mm perimeters:
% 1) printable hollow box spar (no skin)
% 2) skin + box hybrid: the skin gives the aerodynamic shape anyway, so it
%    is counted structurally, and a box spar is added inside it where needed
% 3) skin only (reference: should match the stepped skin of OptimiseSkin.m)
%
%Single-objective problem statement (notes p.17-21):
%   minimise    V = sum_i q_i A_i                  volume (trapezoidal rule)
%   by varying  the section at every node i:       skin perimeters, box width,
%               box wall perimeters (discrete, notes p.11), box centre xs
%   subject to  sigma_i/sigAllow - 1 <= 0          stress at every node
%               sections only stay the same or shrink from root to tip (printing)
%Tip deflection is reported as a check, not optimised (notes p.79).
%
%Method: the design variables are discrete (whole perimeters, a table of
%widths), so gradient-based fmincon cannot be used directly. Instead:
% - every option in the table is evaluated at every node (brute force over
%   a small discrete set, notes p.26), which is cheap because volume and
%   stress are just sums/checks per node
% - the "only shrink" rule links each node to the next one, so the nodes
%   are chosen together by dynamic programming: working from tip to root,
%   each node keeps the lightest way to finish the blade outboard of it.
%   This gives the exact lightest printable design in one pass.
% - the box centre xs is the same at every node (straight webs) and is
%   swept over a list, keeping the lightest (brute force, notes p.26).

N = 100; %the beam solver discretises beam into 100 nodes
c = 0.15; %define the chord length - KEEP FIXED
R = 1.5; %blade length (m)
r = linspace(0,R,N)'; %radial position of every node, root (r = 0) to tip (r = R)
E = 2e9; %Young's modulus of ABS (Pa)
%Deflections are in mm, the units the supplied beamBending returns (as in the original
%RunBeamBending.m). The unit-load integral below uses the same formula, so it gives mm too.

SF = 1.5; %safety factor on yield
sigAllow = 30e6/SF; %allowable bending stress (Pa)
tLine = 0.4; %printed perimeter (extrusion line) width (mm)
nSkin = 2:7; %skin perimeters: 0.8-2.8 mm (mid-surface skin model within 1.5%, as OptimiseSkin.m)
nBox = 2:25; %box wall perimeters: 0.8-10 mm (as CompareConcepts.m)
wList = 10:2:150; %box widths (mm); 10 mm minimum keeps a real closed box rather than a sliver
xsList = 0.24:0.02:0.44; %box centre positions swept (x/c)

%the blade is statically determinate, so the bending moment M(r) does not
%depend on I - ONE solver call gives M for every design
[~, M] = beamBending(N,1e-8*ones(1,N));
absM = abs(M(:));
q = ([diff(r); 0] + [0; diff(r)])/2; %trapezoidal weights, so sums equal trapz exactly
G = q.*absM.*(R - r)/E; %unit-load deflection weight of each node: delta = sum(G./I)
Sreq = absM/sigAllow; %required section modulus I/ymax at each node

kinds = {'box','hybrid','skin'};
names = {'printable box','skin + box hybrid','skin only'};

%% lightest printable design of each concept (sweep the box centre)
tic
for k = 1:numel(kinds)
    if strcmp(kinds{k},'skin'), xsK = NaN; else, xsK = xsList; end
    F(k).Vxs = zeros(size(xsK)); F(k).V = inf; %#ok<SAGROW>
    for j = 1:numel(xsK)
        O = buildOptions(kinds{k}, xsK(j), c, tLine, nSkin, nBox, wList);
        [idx, V, D] = lightestDesign(O, G, Sreq, q);
        F(k).Vxs(j) = V;
        if V < F(k).V
            F(k).V = V; F(k).D = D; F(k).xs = xsK(j); F(k).idx = idx; F(k).opt = O;
        end
    end
    F(k).xsList = xsK;
end
fprintf('Lightest printable designs for %d concepts in %.1f s (1 beamBending call)\n', numel(kinds), toc)

%% results
fprintf('\nLightest stress-feasible printable designs:\n')
fprintf(' %-18s  volume (cm^3)  tip defl (mm)  xs/c   root: skin  box w  box wall   mid (r=0.75): skin  box w  box wall\n', '')
m = round(N/2);
for k = 1:numel(kinds)
    s = describe(F(k).opt, F(k).idx);
    fprintf(' %-18s     %6.1f        %6.3f      %5.2f   %4.1f mm  %3.0f mm  %4.1f mm       %4.1f mm  %3.0f mm  %4.1f mm\n', ...
            names{k}, F(k).V, F(k).D, F(k).xs, s.ts(1), s.w(1), s.tb(1), s.ts(m), s.w(m), s.tb(m))
end

%check against the full solver (deflection) and every node (stress)
fprintf('\nCheck against beamBending (and nodal stress):\n')
for k = 1:numel(kinds)
    O = F(k).opt; idx = F(k).idx;
    delta = beamBending(N,O.I(idx)');
    fprintf('  %-18s fast %.4f mm, solver %.4f mm, max stress %.2f MPa (allow %.0f)\n', ...
            names{k}, F(k).D, abs(delta(end)), max(absM./O.S(idx))/1e6, sigAllow/1e6)
end

%printability: how often the section changes along the span
fprintf('\nSection changes along the span:\n')
for k = 1:numel(kinds)
    s = describe(F(k).opt, F(k).idx);
    fprintf('  %-18s skin %d, box width %d, box wall %d changes\n', names{k}, ...
            nnz(diff(s.ts)), nnz(diff(s.w)), nnz(diff(s.tb)))
end

%% figures
cols = {[0 0.45 0.74], [0.47 0.67 0.19], [0.3 0.3 0.3]};

%volume of the lightest printable design of each concept
figure(1)
vals = [F.V];
bh = bar(categorical(names,names), vals); bh.FaceColor = 'flat';
bh.CData = vertcat(cols{:});
text(1:numel(vals), vals, compose('%.0f cm^3',vals), 'HorizontalAlignment','center','VerticalAlignment','bottom')
ylabel('volume, cm^3'), grid on
title(sprintf('Lightest printable designs in whole %.1f mm perimeters, SF = %.1f', tLine, SF))

%spanwise layout of the lightest hybrid
figure(2)
s = describe(F(2).opt, F(2).idx);
subplot(2,1,1)
stairs(r, s.ts, 'k-', 'LineWidth',1.5), hold on
stairs(r, s.tb, '-', 'Color',cols{2}, 'LineWidth',1.5)
ylabel('thickness, mm'), grid on, legend('skin','box wall')
title(sprintf('Lightest skin + box hybrid: V = %.0f cm^3, x_s/c = %.2f', F(2).V, F(2).xs))
subplot(2,1,2)
stairs(r, s.w, '-', 'Color',cols{2}, 'LineWidth',1.5)
xlabel('distance from root, m'), ylabel('box width, mm'), grid on

%cross-sections of the lightest printable box and hybrid
figure(3)
rShow = [0 0.5 1.0];
for k = 1:2
    s = describe(F(k).opt, F(k).idx);
    for j = 1:numel(rShow)
        [~, n] = min(abs(r - rShow(j)));
        subplot(2,numel(rShow),(k-1)*numel(rShow) + j)
        drawSection(s.ts(n), s.w(n), s.tb(n), F(k).xs, c)
        title(sprintf('%s, r = %.1f m', names{k}, r(n)))
    end
end

%effect of the box centre position (1-D brute-force sweep)
figure(4)
plot(F(1).xsList, F(1).Vxs, 's-', 'Color',cols{1}, 'MarkerFaceColor',cols{1}), hold on
plot(F(2).xsList, F(2).Vxs, 'o-', 'Color',cols{2}, 'MarkerFaceColor',cols{2})
xlabel('box centre x_s / c'), ylabel('lightest volume, cm^3'), grid on
legend(names(1:2), 'Location','best')
title('Box position sweep')


%full grid of printable section options for one concept (and box centre
%xs): skin perimeters x box width x box wall perimeters, in ndgrid order.
%Options whose walls would overlap (2t > w) are marked not ok.
function O = buildOptions(kind, xs, c, tLine, nSkin, nBox, wList)
    switch kind
        case 'skin',   nsL = nSkin; wL = 0;          nbL = 0;
        case 'box',    nsL = 0;     wL = wList;      nbL = nBox;
        case 'hybrid', nsL = nSkin; wL = [0 wList];  nbL = nBox; %w = 0: no box at this node
    end
    O.sz = [numel(nsL) numel(wL) numel(nbL)];
    [ns, w, nb] = ndgrid(nsL, wL, nbL);
    [s1, s2, s3] = ndgrid(1:O.sz(1), 1:O.sz(2), 1:O.sz(3));
    O.sub = [s1(:) s2(:) s3(:)]; %grid position, used for the "only shrink" rule
    ns = ns(:); w = w(:); nb = nb(:);
    O.ok = w == 0 | 2*nb*tLine <= w; %box walls cannot overlap
    ts = ns*tLine*1e-3; tb = nb*tLine*1e-3; wm = w*1e-3; %m

    [Is, As, ymax] = skinSection(ts,c);
    [Ib, Ab, yb] = boxSection(wm, tb, xs, c, ts); %box inside the skin (inset by ts)
    Ib(w == 0) = 0; Ab(w == 0) = 0;
    if strcmp(kind,'box'), ymax = yb; end
    O.I = Is + Ib; O.A = As + Ab; O.S = O.I./ymax;
    O.ts = ts*1e3; O.w = w; O.tb = tb*1e3.*(w > 0); %mm, for reporting
    O.ok = O.ok & O.I > 0;
end

%lightest printable design: minimise sum_i q_i A_i with every node meeting
%its stress limit and the section only staying the same or shrinking from
%root to tip. Dynamic programming, tip to root: g(i,s) is the lightest
%volume of nodes i..N when node i uses section s.
function [idx, V, D] = lightestDesign(O, G, Sreq, q)
    N = numel(G);
    cost = q.*O.A'; %volume of each option at each node
    cost(O.S' < Sreq | ~O.ok') = inf; %overstressed or unbuildable options are not allowed
    g = zeros(size(cost)); g(N,:) = cost(N,:);
    for i = N-1:-1:1
        h = reshape(g(i+1,:), O.sz);
        h = cummin(cummin(cummin(h,1),2),3); %lightest next section that is no larger in every dimension
        g(i,:) = cost(i,:) + h(:)';
    end
    idx = zeros(N,1);
    [gmin, idx(1)] = min(g(1,:));
    if isinf(gmin), error('no printable design carries the moment - widen the option table'), end
    for i = 2:N %follow the lightest choices from root to tip
        v = g(i,:); v(any(O.sub > O.sub(idx(i-1),:), 2)) = inf;
        [~, idx(i)] = min(v);
    end
    V = sum(q.*O.A(idx))*1e6; %cm^3
    D = sum(G./O.I(idx)); %tip deflection, mm (reported only)
end

%skin thickness, box width and box wall (mm) at every node of a design
function s = describe(O, idx)
    s.ts = O.ts(idx); s.w = O.w(idx); s.tb = O.tb(idx);
end

%NACA0012 half-thickness y_t at chordwise position x (m)
function yt = naca0012(x,c)
    xc = x/c;
    yt = c*5*0.12*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 ...
                   + 0.2843*xc.^3 - 0.1015*xc.^4);
end

%thin structural skin (same as RunBeamBending.m and OptimiseSkin.m); t = 0 gives no skin
function [I, A, ymax] = skinSection(t,c)
    nPanels = 200;
    beta = linspace(0,pi,nPanels+1);
    x = c*(1 - cos(beta))/2;
    y = naca0012(x,c);
    Lp = hypot(diff(x),diff(y));
    theta = atan2(diff(y),diff(x));
    d = (y(1:end-1) + y(2:end))/2;
    %the material sits on the skin MID-surface, t/2 inside the outer surface
    %along the panel normal (height d - t cos(theta)/2). A 2 mm skin is not
    %thin next to the 9 mm half-depth: using the outer surface overstates I by 33%.
    Iu = 2*sum(Lp.^3.*sin(theta).^2/12 + Lp.*d.^2); %I per unit thickness, outer surface (m^3)
    Ic = 2*sum(Lp.*d.*cos(theta)); %mid-surface offset corrections
    Icc = 2*sum(Lp.*cos(theta).^2)/4;
    perimeter = 2*sum(Lp); %full aerofoil perimeter (m)
    I = Iu*t - Ic*t.^2 + Icc*t.^3; %within 0.6% of exact for t <= 2 mm
    A = perimeter*t - pi*t.^2; %inward offset band (Steiner), within 3% (conservative) for t <= 2 mm
    ymax = max(y)*ones(size(t));
end

%hollow box spar (as RunBeamBending.m) whose outer faces sit an inset below
%the aerofoil surface: inset = skin thickness for the hybrid, 0 for a bare box
function [I, A, ymax] = boxSection(w,t,xs,c,inset)
    nStrips = 50;
    w = w(:); t = t(:); inset = inset(:);
    x1 = max(xs*c - w/2, 0);
    x2 = min(xs*c + w/2, c);
    [Io, Ao, ymax] = stripSum(x1, x2, inset, c, nStrips);
    [Ii, Ai] = stripSum(x1 + t, max(x2 - t, x1 + t), inset + t, c, nStrips);
    I = Io - Ii; A = Ao - Ai;
end

function [I, A, ymax] = stripSum(a,b,inset,c,nStrips)
    s = linspace(0,1,nStrips+1);
    edges = a + (b - a).*s;
    bw = diff(edges,1,2);
    mid = (edges(:,1:end-1) + edges(:,2:end))/2;
    d = max(2*(naca0012(mid,c) - inset), 0);
    I = sum(bw.*d.^3/12, 2);
    A = sum(bw.*d, 2);
    ymax = max(naca0012(edges,c), [], 2);
end

%draw one section: skin of thickness ts, box of width w and wall tb (mm)
function drawSection(ts, w, tb, xs, c)
    warning('off','MATLAB:polyshape:repairedBySimplify')
    xo = linspace(0,c,300); yo = naca0012(xo,c);
    outline = polyshape([xo fliplr(xo)], [yo -fliplr(yo)]);
    sec = polyshape();
    if ts > 0, sec = subtract(outline, polybuffer(outline,-ts*1e-3)); end
    if w > 0
        t = tb*1e-3; i0 = ts*1e-3;
        xb = linspace(max(xs*c - w/2e3,0), min(xs*c + w/2e3,c), 100);
        yb = max(naca0012(xb,c) - i0, 0);
        box = polyshape([xb fliplr(xb)], [yb -fliplr(yb)]);
        xi = linspace(xb(1) + t, xb(end) - t, 100); yi = max(naca0012(xi,c) - i0 - t, 0);
        if xi(end) > xi(1) && any(yi > 0)
            box = subtract(box, polyshape([xi fliplr(xi)], [yi -fliplr(yi)]));
        end
        sec = union(sec, box);
    end
    plot(scale(sec,1e3), 'FaceColor',[0.7 0.8 1], 'EdgeColor','b', 'FaceAlpha',1), hold on
    plot(xo*1e3, yo*1e3, 'k-', xo*1e3, -yo*1e3, 'k-')
    axis equal, grid on, xlim([0 c*1e3])
    xlabel('x, mm'), ylabel('y, mm')
end
