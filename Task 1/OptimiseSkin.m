clc, clear all, close all

%Lightest thin structural skin (Approach 2) that does not yield.
%Single-objective problem statement (notes p.17-21):
%   minimise    V(x)                            skin volume
%   by varying  x                               skin thickness along the span
%   subject to  sigma_max(x)/sigAllow - 1 <= 0  bending stress, safety factor 1.5
%               tMin <= t <= tMax                print limits
%Tip deflection is NOT an objective: every design deflects about 1 mm, far
%inside any sensible limit, so it is only reported as a check (notes p.79:
%turn extra objectives into constraints or checks).
%
%Two parameterisations of the thickness:
% 1) linear taper, x = [tRoot tTip]: 2 continuous variables, solved with
%    fmincon from several starts, then rounded UP to whole 0.4 mm perimeters
% 2) stepped skin in whole perimeters (2.8, 2.4, ... 0.8 mm), printable by
%    construction. The variables are where each step ends. Making a step
%    longer only adds material, so the lightest design makes every step as
%    short as the stress allows: the optimum lies on the constraint boundary
%    (notes p.18-19, Fig. 11) and is found exactly, with no search.

N = 100; %the beam solver discretises beam into 100 nodes
c = 0.15; %define the chord length - KEEP FIXED
R = 1.5; %blade length (m)
r = linspace(0,R,N); %radial position of every node, root (r = 0) to tip (r = R)
E = 2e9; %Young's modulus of ABS (Pa)
%Deflections are in mm, the units the supplied beamBending returns (as in the original
%RunBeamBending.m). The unit-load integral below uses the same formula, so it gives mm too.

SF = 1.5; %safety factor on yield
sigAllow = 30e6/SF; %allowable bending stress (Pa)
tMin = 0.8; %minimum skin thickness (mm) - two perimeters of 0.4 mm nozzle extrusion for structure
tMax = 2.8; %maximum skin thickness (mm) - 7 perimeters; the mid-surface skin model is within 1.5% (conservative) up to here
tLine = 0.4; %printed perimeter (extrusion line) width (mm) - skin thickness must be a multiple of this
nStart = 5; %starting points for fmincon (multistart, notes p.58)

%the blade is statically determinate, so the bending moment M(r) does not
%depend on I - ONE solver call gives M for every design
[~, M] = beamBending(N,1e-8*ones(1,N));
absM = abs(M);
opts = optimoptions('fmincon','Display','off');
rng(0) %repeatable random starts

%% 1) linear taper, x = [tRoot tTip] in mm (mm so both variables are of order 1, notes p.39)
evalSkin = @(x) skinModel(taper(x,r),c,r,M,E); %x -> [volume cm^3, tip deflection mm, max stress Pa]
lb = [tMin tMin]; ub = [tMax tMax];
volume = @(x) pick(evalSkin,x,1);
nonlcon = @(x) deal(pick(evalSkin,x,3)/sigAllow - 1, []); %normalised stress constraint, <= 0 is feasible

starts = [ub; lb + (ub - lb).*rand(nStart-1,2)];
fBest = inf;
for i = 1:nStart
    [x, f, flag] = fmincon(volume,starts(i,:),[],[],[],[],lb,ub,nonlcon,opts);
    if flag > 0 && f < fBest, xT = x; fBest = f; end
end
[VT, DT, ST] = evalSkin(xT);

%printable version: round the taper UP to whole perimeters (only adds
%material, so it stays stress-feasible)
tR = roundUp(taper(xT,r),tLine);
[VR, DR, SR] = skinModel(tR,c,r,M,E);

%% 2) stepped skin, optimised directly in whole perimeters
nLev = round(tMax/tLine):-1:round(tMin/tLine); %perimeters in each step, root to tip (7 6 5 4 3 2)
tLev = nLev*tLine*1e-3; %thickness of each step (m)
K = numel(tLev);
[Ilev, ~, ym] = skinSection(tLev,c); %I and outer fibre of each step
ym = ym(1);

%stress: |M| falls from root to tip, so each step is most stressed where it
%starts. Step k+1 can only start once |M| has dropped to the moment it can carry.
Mallow = sigAllow*Ilev/ym; %largest moment each step can carry (N m)
if absM(1) > Mallow(1)
    error('a %.1f mm root skin cannot carry the root moment - raise tMax', tLev(1)*1e3)
end
bMin = zeros(1,K-1); %earliest radius where each thinner step can start
for j = 1:K-1
    if absM(1) > Mallow(j+1)
        bMin(j) = fzero(@(x) interp1(r,absM,x,'pchip') - Mallow(j+1), [0 R]);
    end
end
tS = stepProfile(bMin,tLev,r); %lightest stepped skin: every step as short as stress allows
[VS, DS, SS] = skinModel(tS,c,r,M,E);

%% results
fprintf('Lightest stress-feasible skins (allowable stress %.0f MPa)\n', sigAllow/1e6)
fprintf('  %-28s V = %6.1f cm^3, tip defl = %.3f mm, max stress = %.2f MPa, tRoot/tTip = %.2f/%.2f mm\n', ...
        'linear taper (continuous)', VT, DT, ST/1e6, xT)
fprintf('  %-28s V = %6.1f cm^3, tip defl = %.3f mm, max stress = %.2f MPa\n', ...
        'taper rounded up to 0.4 mm', VR, DR, SR/1e6)
fprintf('  %-28s V = %6.1f cm^3, tip defl = %.3f mm, max stress = %.2f MPa\n', ...
        'stepped skin (printable)', VS, DS, SS/1e6)
fprintf('  stepped skin is %.1f%% lighter than the rounded taper\n', 100*(1 - VS/VR))
stepNames = arrayfun(@(k) sprintf('%.1f->%.1f mm', tLev(k)*1e3, tLev(k+1)*1e3), 1:K-1, 'UniformOutput', false);
fprintf('  step positions: %s\n', strjoin(compose('%s at r = %.3f m', string(stepNames'), bMin'), ', '))

%check against the full solver
fprintf('\nCheck against beamBending:\n')
designs = {taper(xT,r), tR, tS}; names = {'linear taper','rounded taper','stepped skin'}; Dfast = [DT DR DS];
for k = 1:3
    delta = beamBending(N,skinSection(designs{k},c));
    fprintf('  %-14s fast %.4f mm, solver %.4f mm\n', names{k}, Dfast(k), abs(delta(end)))
end

%% figures
%design space of the linear taper: volume contours, stress-infeasible region
%and the optimum on the constraint boundary (compare notes Fig. 11 and 13)
[g1, g2] = meshgrid(linspace(tMin,tMax,60));
Vg = zeros(size(g1)); Sg = Vg;
for k = 1:numel(g1)
    [Vg(k), ~, Sg(k)] = evalSkin([g1(k) g2(k)]);
end
figure(1)
contourf(g1,g2,Vg,20,'LineColor','none'), hold on
cb = colorbar; cb.Label.String = 'volume, cm^3';
contourf(g1,g2,double(Sg > sigAllow),[0.5 0.5],'FaceColor',[1 1 1],'FaceAlpha',0.6,'LineColor','none')
contour(g1,g2,Sg/sigAllow,[1 1],'r-','LineWidth',2)
plot(xT(1),xT(2),'kp','MarkerSize',16,'MarkerFaceColor','y')
xlabel('root thickness t_{root}, mm'), ylabel('tip thickness t_{tip}, mm')
legend('volume','violates stress (shaded)','\sigma_{max} = \sigma_{allow}','optimum','Location','northwest')
title(sprintf('Linear taper: lightest feasible skin V = %.0f cm^3', VT))

%thickness and stress along the span of the three designs
figure(2)
subplot(2,1,1)
plot(r,taper(xT,r)*1e3,'k-', r,tR*1e3,'--', r,tS*1e3,'b-','LineWidth',1.5)
ylabel('skin thickness, mm'), grid on
legend(sprintf('linear taper, %.0f cm^3',VT), sprintf('rounded up, %.0f cm^3',VR), ...
       sprintf('stepped, %.0f cm^3',VS), 'Location','northeast')
subplot(2,1,2)
plot(r,absM*ym./skinSection(taper(xT,r),c)/1e6,'k-', r,absM*ym./skinSection(tR,c)/1e6,'--', ...
     r,absM*ym./skinSection(tS,c)/1e6,'b-','LineWidth',1.5), hold on
yline(sigAllow/1e6,'r--','allowable')
xlabel('distance from root, m'), ylabel('bending stress, MPa'), grid on

%linear taper of skin thickness from root to tip, x = [tRoot tTip] in mm -> t in m
function t = taper(x,r)
    t = (x(1) + (x(2) - x(1))*r/r(end))*1e-3;
end

%round thickness t (m) UP to a whole number of perimeters of width tLine (mm)
%(small tolerance so values already on a step are not pushed up a step)
function t = roundUp(t,tLine)
    t = ceil(t*1e3/tLine - 1e-9)*tLine*1e-3;
end

%thickness at every node of a stepped skin: step k (thickness tLev(k)) runs
%from b(k-1) to b(k), with b(0) = 0 at the root
function t = stepProfile(b,tLev,r)
    starts = [0 b];
    t = zeros(size(r));
    for k = 1:numel(tLev)
        t(r >= starts(k)) = tLev(k);
    end
end

%evaluate one skin thickness distribution t (m) without calling beamBending:
%volume = integral of area along the span
%tip deflection from the unit-load (moment-area) method:
%   delta_tip = integral M(r)(R - r)/(E I(r)) dr  (matches beamBending to 0.02%)
%stress = |M| y_max / I at every node
function [V, dTip, sigMax] = skinModel(t,c,r,M,E)
    [I, A, ymax] = skinSection(t,c);
    V = trapz(r,A)*1e6;
    dTip = trapz(r, abs(M).*(r(end) - r)./(E*I));
    sigMax = max(abs(M).*ymax./I);
end

%return output number n of a multi-output function (for fmincon handles)
function v = pick(f,x,n)
    [out{1:3}] = f(x);
    v = out{n};
end

%NACA0012 half-thickness y_t at chordwise position x (m)
function yt = naca0012(x,c)
    xc = x/c;
    yt = c*5*0.12*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 ...
                   + 0.2843*xc.^3 - 0.1015*xc.^4);
end

%Approach 2: thin structural skin of thickness t (same as RunBeamBending.m)
function [I, A, ymax] = skinSection(t,c)
    nPanels = 200; %panels per surface
    beta = linspace(0,pi,nPanels+1);
    x = c*(1 - cos(beta))/2; %cosine spacing: more panels round the curved leading edge
    y = naca0012(x,c);
    Lp = hypot(diff(x),diff(y)); %panel lengths
    theta = atan2(diff(y),diff(x)); %panel angles to the chord
    d = (y(1:end-1) + y(2:end))/2; %panel centroid heights above the neutral axis
    %the material sits on the skin MID-surface, t/2 inside the outer surface
    %along the panel normal (height d - t cos(theta)/2). A 2 mm skin is not
    %thin next to the 9 mm half-depth: using the outer surface overstates I by 33%.
    Iu = 2*sum(Lp.^3.*sin(theta).^2/12 + Lp.*d.^2); %I per unit thickness, outer surface (m^3)
    Ic = 2*sum(Lp.*d.*cos(theta)); %mid-surface offset corrections
    Icc = 2*sum(Lp.*cos(theta).^2)/4;
    perimeter = 2*sum(Lp); %full aerofoil perimeter (m)
    I = Iu*t - Ic*t.^2 + Icc*t.^3; %within 0.6% of exact for t <= 2 mm
    A = perimeter*t - pi*t.^2; %inward offset band (Steiner), within 3% (conservative) for t <= 2 mm
    ymax = max(y)*ones(size(t)); %furthest fibre = max half-thickness (0.06c)
end
