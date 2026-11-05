clc, clear all, close all

%Pareto front of volume vs tip deflection for a linearly tapered skin
%(Approach 2), using the epsilon-constraint method:
%   minimise V(x)  subject to  sigma_max(x) <= sigAllow,  delta_tip(x) <= eps
%swept over eps from the stiffest to the lightest feasible design.
%Each optimum is then rounded UP to a whole number of printed perimeters
%(thickness steps) to give the printable design.

N = 100; %the beam solver discretises beam into 100 nodes
c = 0.15; %define the chord length - KEEP FIXED
R = 1.5; %blade length (m)
r = linspace(0,R,N); %radial position of every node, root (r = 0) to tip (r = R)
E = 2e9; %Young's modulus of ABS (Pa)

SF = 1.5; %safety factor on yield
sigAllow = 30e6/SF; %allowable bending stress (Pa)
tMin = 0.8; %minimum skin thickness (mm) - two perimeters of 0.4 mm nozzle extrusion for structure
tMax = 2.0; %maximum skin thickness (mm) - keeps t <= ~20% of the 9 mm half-depth so the thin-wall I is valid
tLine = 0.4; %printed perimeter (extrusion line) width (mm) - skin thickness must be a multiple of this
nFront = 15; %number of points on the Pareto front

%the blade is statically determinate, so the bending moment M(r) does not
%depend on I - ONE solver call gives M for every design
[~, M] = beamBending(N,1e-8*ones(1,N));
evalSkin = @(x) skinModel(taper(x,r),c,r,M,E); %x = [tRoot tTip] in mm -> [volume cm^3, tip deflection m, max stress Pa]
evalPrint = @(x) skinModel(stepped(taper(x,r),tLine),c,r,M,E); %same, after rounding up to whole perimeters

opts = optimoptions('fmincon','Display','off');
lb = [tMin tMin]; ub = [tMax tMax];
volume = @(x) pick(evalSkin,x,1);
stressCon = @(x) pick(evalSkin,x,3)/sigAllow - 1; %normalised: <= 0 is feasible
deflCon = @(x,eps) pick(evalSkin,x,2)/eps - 1;

tic
%anchor 1: lightest design that meets the stress constraint
xLight = fmincon(volume,[tMax tMax],[],[],[],[],lb,ub,@(x) deal(stressCon(x),[]),opts);
[~, dLight] = evalSkin(xLight);
%anchor 2: stiffest design = thickest skin everywhere
[~, dStiff] = evalSkin(ub);

%epsilon-constraint sweep: minimise volume for a range of deflection limits
epsList = linspace(dStiff,dLight,nFront);
X = zeros(nFront,2); V = zeros(1,nFront); D = V; S = V;
x0 = ub; %start from the stiffest design and warm-start each solve from the last
for k = 1:nFront
    nonlcon = @(x) deal([stressCon(x); deflCon(x,epsList(k))], []);
    X(k,:) = fmincon(volume,x0,[],[],[],[],lb,ub,nonlcon,opts);
    [V(k), D(k), S(k)] = evalSkin(X(k,:));
    x0 = X(k,:);
end

%printable designs: round each continuous optimum UP to whole perimeters.
%Rounding up only adds material, so stress stays within the limit and
%deflection can only decrease.
Vp = zeros(1,nFront); Dp = Vp; Sp = Vp;
for k = 1:nFront
    [Vp(k), Dp(k), Sp(k)] = evalPrint(X(k,:));
end
fprintf('Pareto front: %d points in %.2f s (1 beamBending call)\n', nFront, toc)

%print the front
fprintf('\n          continuous optimum                    |  printable (%.1f mm steps)\n', tLine)
fprintf(' tip defl (m)  volume (cm^3)  tRoot  tTip (mm)  |  tip defl (m)  volume (cm^3)  max stress (MPa)\n')
fprintf('   %6.3f        %7.1f      %5.2f  %5.2f     |    %6.3f        %7.1f          %5.2f\n', [D; V; X'; Dp; Vp; Sp/1e6])

%check the fast deflection against the full solver at both ends and the middle
fprintf('\nCheck against beamBending:\n')
for k = [1 ceil(nFront/2) nFront]
    I = skinSection(stepped(taper(X(k,:),r),tLine),c);
    delta = beamBending(N,I);
    fprintf('  printable point %2d: fast %.4f m, solver %.4f m\n', k, Dp(k), abs(delta(end)))
end

%brute-force grid of designs to show the whole design space behind the front
[g1, g2] = meshgrid(linspace(tMin,tMax,30));
Vg = zeros(size(g1)); Dg = Vg; Sg = Vg;
for k = 1:numel(g1)
    [Vg(k), Dg(k), Sg(k)] = evalSkin([g1(k) g2(k)]);
end
ok = Sg <= sigAllow;

%plot the Pareto front
figure(1)
scatter(Dg(ok),Vg(ok),12,[0.7 0.7 0.7],'filled'), hold on
scatter(Dg(~ok),Vg(~ok),12,[1 0.7 0.7],'filled')
plot(D,V,'ko-','LineWidth',1.5,'MarkerFaceColor','k')
plot(Dp,Vp,'bs-','LineWidth',1.5,'MarkerFaceColor','b')
xline(0.15*R,'b--','15% of span')
xlabel('tip deflection, m'), ylabel('volume, cm^3'), grid on
legend('feasible designs','violates stress','Pareto front (continuous)', ...
       sprintf('printable (%.1f mm steps)',tLine),'Location','northeast')
title(sprintf('Linearly tapered skin, SF = %.1f, t = %.1f-%.1f mm', SF, tMin, tMax))

%plot the thickness of each Pareto-optimal design
figure(2)
plot(D,X(:,1),'ko-',D,X(:,2),'bs-','MarkerFaceColor','auto')
xlabel('tip deflection, m'), ylabel('skin thickness, mm'), grid on
legend('t_{root}','t_{tip}')

%plot the continuous and printable thickness of the lightest design
figure(3)
plot(r,taper(X(end,:),r)*1e3,'k-', r,stepped(taper(X(end,:),r),tLine)*1e3,'b-','LineWidth',1.5)
xlabel('distance, m'), ylabel('skin thickness, mm'), grid on
legend('continuous optimum',sprintf('printable (%.1f mm steps)',tLine))
title(sprintf('Lightest design: V = %.0f cm^3 continuous, %.0f cm^3 printable', V(end), Vp(end)))


%linear taper of skin thickness from root to tip, x = [tRoot tTip] in mm -> t in m
function t = taper(x,r)
    t = (x(1) + (x(2) - x(1))*r/r(end))*1e-3;
end

%round thickness t (m) UP to a whole number of perimeters of width tLine (mm)
%(small tolerance so values already on a step are not pushed up a step)
function t = stepped(t,tLine)
    t = ceil(t*1e3/tLine - 1e-9)*tLine*1e-3;
end

%evaluate one skin thickness distribution t (m) without calling beamBending:
%volume = integral of area along the span
%tip deflection from the unit-load (moment-area) method:
%   delta_tip = integral M(r)(R - r)/(E I(r)) dr  (matches beamBending to 0.01%)
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
    Iunit = 2*sum(Lp.^3.*sin(theta).^2/12 + Lp.*d.^2); %I per unit skin thickness (m^3)
    perimeter = 2*sum(Lp); %full aerofoil perimeter (m)
    I = Iunit*t;
    A = perimeter*t;
    ymax = max(y)*ones(size(t)); %furthest fibre = max half-thickness (0.06c)
end
