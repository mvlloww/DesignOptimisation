clc, clear all, close all

%Pareto front of volume vs tip deflection for the thin structural skin
%(Approach 2), using the epsilon-constraint method:
%   minimise V(x)  subject to  sigma_max(x) <= sigAllow,  delta_tip(x) <= eps
%swept over eps from the stiffest to the lightest feasible design.
%
%Two design families are optimised:
% 1) linearly tapered skin, x = [tRoot tTip] - the continuous reference
% 2) STEPPED skin, the printable design: the wall is a whole number of
%    printed perimeters, dropping one perimeter at a time from root to tip.
%    The thickness of each step is fixed (2.0, 1.6, 1.2, 0.8 mm) and the
%    design variables are the radial positions where each step ends.
%    This optimises the printable design directly, rather than rounding a
%    continuous optimum up to whole perimeters afterwards.

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
nFront = 15; %number of points on each Pareto front

%the blade is statically determinate, so the bending moment M(r) does not
%depend on I - ONE solver call gives M for every design
[~, M] = beamBending(N,1e-8*ones(1,N));
opts = optimoptions('fmincon','Display','off');

tic
%% 1) linearly tapered skin (continuous reference)
evalSkin = @(x) skinModel(taper(x,r),c,r,M,E); %x = [tRoot tTip] in mm -> [volume cm^3, tip deflection m, max stress Pa]
lb = [tMin tMin]; ub = [tMax tMax];
volume = @(x) pick(evalSkin,x,1);
stressCon = @(x) pick(evalSkin,x,3)/sigAllow - 1; %normalised: <= 0 is feasible
deflCon = @(x,eps) pick(evalSkin,x,2)/eps - 1;

%anchor 1: lightest design that meets the stress constraint
xLight = fmincon(volume,[tMax tMax],[],[],[],[],lb,ub,@(x) deal(stressCon(x),[]),opts);
[~, dLight] = evalSkin(xLight);
%anchor 2: stiffest design = thickest skin everywhere
[~, dStiff] = evalSkin(ub);

%epsilon-constraint sweep: minimise volume for a range of deflection limits
epsList = linspace(dStiff,dLight,nFront);
X = zeros(nFront,2); V = zeros(1,nFront); D = V;
x0 = ub; %start from the stiffest design and warm-start each solve from the last
for k = 1:nFront
    nonlcon = @(x) deal([stressCon(x); deflCon(x,epsList(k))], []);
    X(k,:) = fmincon(volume,x0,[],[],[],[],lb,ub,nonlcon,opts);
    [V(k), D(k)] = evalSkin(X(k,:));
    x0 = X(k,:);
end

%for comparison only: the old approach of rounding each taper optimum UP to
%whole perimeters after optimising
Vr = zeros(1,nFront); Dr = Vr;
for k = 1:nFront
    [Vr(k), Dr(k)] = skinModel(roundUp(taper(X(k,:),r),tLine),c,r,M,E);
end

%% 2) stepped skin, optimised directly in whole perimeters
nLev = round(tMax/tLine):-1:round(tMin/tLine); %perimeters in each step, root to tip (5 4 3 2)
tLev = nLev*tLine*1e-3; %thickness of each step (m)
K = numel(tLev); %number of steps -> K-1 step positions b(1..K-1) are the design variables

[Iu, Pu, ym] = skinSection(1,c); %I, perimeter and outer fibre per unit thickness
absM = abs(M);
%Because I = Iu*t, volume and deflection of a stepped skin are exact sums over the steps:
%   V     = Pu * sum_k t_k (b_k - b_k-1)
%   delta = sum_k [F(b_k) - F(b_k-1)]/(E Iu t_k),  F(r) = cumulative integral of |M|(R - r)
F = cumtrapz(r, absM.*(R - r));
Fq = @(x) interp1(r,F,x,'pchip'); %smooth in the step positions, so fmincon sees a smooth problem
stepV = @(b) Pu*sum(tLev.*diff([0 b R]))*1e6; %cm^3
stepD = @(b) sum(diff(Fq([0 b R]))./(E*Iu*tLev)); %m

%stress: |M| falls from root to tip, so each step is most stressed where it
%starts. Step k+1 can only start once |M| has dropped to the moment it can
%carry -> the stress constraint becomes a simple lower bound on each b.
Mallow = sigAllow*Iu*tLev/ym; %largest moment each step can carry (N m)
if absM(1) > Mallow(1)
    error('a %.1f mm root skin cannot carry the root moment - raise tMax', tLev(1)*1e3)
end
bMin = zeros(1,K-1);
for j = 1:K-1
    if absM(1) > Mallow(j+1)
        bMin(j) = fzero(@(x) interp1(r,absM,x,'pchip') - Mallow(j+1), [0 R]);
    end
end

%steps must stay in order root to tip: b(j) <= b(j+1)
Aord = -diff(eye(K-1)); bOrd = zeros(K-2,1);
lbS = bMin; ubS = R*ones(1,K-1);

%anchors: making any step longer only adds material, so the lightest
%stress-feasible design has every step as short as stress allows (b = bMin);
%the stiffest is the full 2.0 mm skin (b = R)
dLightS = stepD(lbS);
dStiffS = stepD(ubS);

epsS = linspace(dStiffS,dLightS,nFront);
B = zeros(nFront,K-1); VS = zeros(1,nFront); DS = VS;
b0 = ubS; %start from the stiffest design and warm-start each solve from the last
for k = 1:nFront
    nonlcon = @(b) deal(stepD(b)/epsS(k) - 1, []);
    B(k,:) = fmincon(stepV,b0,Aord,bOrd,[],[],lbS,ubS,nonlcon,opts);
    VS(k) = stepV(B(k,:)); DS(k) = stepD(B(k,:));
    b0 = B(k,:);
end
fprintf('Both Pareto fronts: %d points each in %.2f s (1 beamBending call)\n', nFront, toc)

%% results
fprintf('\nLinear taper (continuous reference)\n')
fprintf(' tip defl (m)  volume (cm^3)  tRoot  tTip (mm)  |  rounded up: tip defl (m)  volume (cm^3)\n')
fprintf('   %6.3f        %7.1f      %5.2f  %5.2f     |              %6.3f        %7.1f\n', [D; V; X'; Dr; Vr])

stepNames = arrayfun(@(k) sprintf('%.1f->%.1f', tLev(k)*1e3, tLev(k+1)*1e3), 1:K-1, 'UniformOutput', false);
fprintf('\nStepped skin, optimised in whole %.1f mm perimeters\n', tLine)
fprintf(' tip defl (m)  volume (cm^3)  step positions r (m): %s\n', strjoin(stepNames,'  '))
fprintf(['   %6.3f        %7.1f      ' repmat('        %5.3f',1,K-1) '\n'], [DS; VS; B'])

fprintf('\nLightest printable design: rounded-up taper %.1f cm^3, optimised steps %.1f cm^3 (%.1f%% less)\n', ...
        Vr(end), VS(end), 100*(1 - VS(end)/Vr(end)))

%check the stepped designs against the full solver on the 100-node grid
%(steps fall between nodes, so a small discretisation difference is expected)
fprintf('\nCheck stepped designs against beamBending:\n')
for k = [1 ceil(nFront/2) nFront]
    t = stepProfile(B(k,:),tLev,r);
    delta = beamBending(N,Iu*t);
    sig = max(absM*ym./(Iu*t));
    fprintf('  point %2d: exact %.4f m, solver %.4f m, max nodal stress %.2f MPa (allow %.0f)\n', ...
            k, DS(k), abs(delta(end)), sig/1e6, sigAllow/1e6)
end

%brute-force grid of linear designs to show the design space behind the fronts
[g1, g2] = meshgrid(linspace(tMin,tMax,30));
Vg = zeros(size(g1)); Dg = Vg; Sg = Vg;
for k = 1:numel(g1)
    [Vg(k), Dg(k), Sg(k)] = evalSkin([g1(k) g2(k)]);
end
ok = Sg <= sigAllow;

%plot the Pareto fronts
figure(1)
scatter(Dg(ok),Vg(ok),12,[0.7 0.7 0.7],'filled'), hold on
scatter(Dg(~ok),Vg(~ok),12,[1 0.7 0.7],'filled')
plot(D,V,'ko-','LineWidth',1.5,'MarkerFaceColor','k')
plot(Dr,Vr,'--','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(DS,VS,'bs-','LineWidth',1.5,'MarkerFaceColor','b')
xline(0.15*R,'b--','15% of span')
xlabel('tip deflection, m'), ylabel('volume, cm^3'), grid on
legend('linear taper: feasible','linear taper: violates stress','linear taper Pareto front (continuous)', ...
       'linear taper rounded up to perimeters', sprintf('stepped skin, optimised in %.1f mm perimeters',tLine), ...
       'Location','northeast')
title(sprintf('Thin skin, SF = %.1f, t = %.1f-%.1f mm', SF, tMin, tMax))

%plot where each step ends for every Pareto-optimal stepped design
figure(2)
plot(DS,B,'s-','MarkerFaceColor','auto'), hold on
yline(bMin,':','Color',[0.5 0.5 0.5])
xlabel('tip deflection, m'), ylabel('step position r, m'), grid on
legend(strcat(stepNames,' mm'),'Location','northeast')
title('Stepped skin: where each perimeter is dropped (dotted = stress limit)')

%plot the thickness of the lightest designs
figure(3)
plot(r,taper(X(end,:),r)*1e3,'k-', r,roundUp(taper(X(end,:),r),tLine)*1e3,'--', ...
     r,stepProfile(B(end,:),tLev,r)*1e3,'b-','LineWidth',1.5)
xlabel('distance, m'), ylabel('skin thickness, mm'), grid on
legend('linear taper (continuous)','linear taper rounded up','stepped (optimised)')
title(sprintf('Lightest designs: %.0f cm^3 rounded-up taper vs %.0f cm^3 optimised steps', Vr(end), VS(end)))


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
    Iunit = 2*sum(Lp.^3.*sin(theta).^2/12 + Lp.*d.^2); %I per unit skin thickness (m^3)
    perimeter = 2*sum(Lp); %full aerofoil perimeter (m)
    I = Iunit*t;
    A = perimeter*t;
    ymax = max(y)*ones(size(t)); %furthest fibre = max half-thickness (0.06c)
end
