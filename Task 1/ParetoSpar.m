clc, clear all, close all

%Fair comparison of structural concepts: every concept is optimised with
%the same method, loads, safety factor and print limits, so the choice of
%concept is based on optimised designs rather than hand-picked ones.
%
%   solid spar  (Approach 1)        x = [xs wRoot wTip]
%   hollow box  (Approach 1, guide) x = [xs wRoot wTip tRoot tTip]
%   thin skin   (Approach 2)        x = [tRoot tTip]
%
%xs = chordwise centre of the spar (x/c), w = spar width (mm), t = wall
%thickness (mm); widths and thicknesses taper linearly from root to tip.
%For each concept an epsilon-constraint Pareto front is found:
%   minimise V(x)  subject to  sigma_max(x) <= sigAllow,  delta_tip(x) <= eps
%Spar position and width make the problem non-convex, so every
%optimisation is run from several starting points (multistart).

N = 100; %the beam solver discretises beam into 100 nodes
c = 0.15; %define the chord length - KEEP FIXED
R = 1.5; %blade length (m)
r = linspace(0,R,N); %radial position of every node, root (r = 0) to tip (r = R)
E = 2e9; %Young's modulus of ABS (Pa)

SF = 1.5; %safety factor on yield
sigAllow = 30e6/SF; %allowable bending stress (Pa)
tMin = 0.8; %minimum wall thickness (mm) - two perimeters of 0.4 mm nozzle extrusion
tMax = 2.0; %maximum SKIN thickness (mm) - thin-wall validity of the skin model
tBoxMax = 10; %maximum box wall thickness (mm) - the strip model is exact for any t
nFront = 15; %number of points on each Pareto front
nStartAnchor = 8; %multistart points for the two anchor designs
nStartSweep = 3; %extra random starts for every point on the front (plus a warm start)

%the blade is statically determinate, so the bending moment M(r) does not
%depend on I - ONE solver call gives M for every design of every concept
[~, M] = beamBending(N,1e-8*ones(1,N));
opts = optimoptions('fmincon','Display','off');
rng(0) %repeatable random starts

%% concept definitions
C(1).name = 'solid spar'; C(1).key = 'spar';
C(1).lb = [0.05 tMin tMin]; C(1).ub = [0.95 c*1e3 c*1e3];
C(1).A = []; C(1).b = [];

C(2).name = 'hollow box'; C(2).key = 'box';
C(2).lb = [0.05 2*tMin 2*tMin tMin tMin]; C(2).ub = [0.95 c*1e3 c*1e3 tBoxMax tBoxMax];
C(2).A = [0 -1 0 2 0; 0 0 -1 0 2]; C(2).b = [0; 0]; %walls cannot overlap: 2t <= w at root and tip

C(3).name = 'thin skin'; C(3).key = 'skin';
C(3).lb = [tMin tMin]; C(3).ub = [tMax tMax];
C(3).A = []; C(3).b = [];

%% Pareto front of every concept
tic
for k = 1:numel(C)
    ev = @(x) evalConcept(C(k).key,x,r,c,M,E); %x -> [volume cm^3, tip deflection m, max stress Pa]
    volume = @(x) pick(ev,x,1);
    deflection = @(x) pick(ev,x,2);
    stressCon = @(x) pick(ev,x,3)/sigAllow - 1; %normalised: <= 0 is feasible
    starts = @(n) C(k).lb + (C(k).ub - C(k).lb).*rand(n,numel(C(k).lb));

    %anchor 1: lightest design that meets the stress constraint
    xLight = bestOf(volume, starts(nStartAnchor), C(k), @(x) deal(stressCon(x),[]), opts);
    %anchor 2: stiffest stress-feasible design
    xStiff = bestOf(deflection, starts(nStartAnchor), C(k), @(x) deal(stressCon(x),[]), opts);
    [~, dLight] = ev(xLight);
    [~, dStiff] = ev(xStiff);

    %epsilon-constraint sweep, stiff to light, warm start plus random starts
    epsList = linspace(dStiff,dLight,nFront);
    X = zeros(nFront,numel(C(k).lb)); V = zeros(1,nFront); D = V; S = V;
    x0 = xStiff;
    for j = 1:nFront
        nonlcon = @(x) deal([stressCon(x); deflection(x)/epsList(j) - 1], []);
        X(j,:) = bestOf(volume, [x0; starts(nStartSweep)], C(k), nonlcon, opts);
        [V(j), D(j), S(j)] = ev(X(j,:));
        x0 = X(j,:);
    end
    C(k).X = X; C(k).V = V; C(k).D = D; C(k).S = S;
end
fprintf('Pareto fronts for %d concepts in %.1f s (1 beamBending call)\n', numel(C), toc)

%% results
fprintf('\n%-11s  %s\n', 'concept', 'lightest stress-feasible design')
for k = 1:numel(C)
    fprintf('%-11s  V = %6.1f cm^3, tip defl = %.3f m (%.0f%% of span), max stress = %.2f MPa, x = [%s]\n', ...
            C(k).name, C(k).V(end), C(k).D(end), 100*C(k).D(end)/R, C(k).S(end)/1e6, ...
            strjoin(compose('%.3g',C(k).X(end,:)),' '))
end
[Vs, Ds] = evalConcept('spar',[0.5 c*1e3 c*1e3],r,c,M,E); %fully filled blade for reference
fprintf('%-11s  V = %6.1f cm^3, tip defl = %.3f m (%.0f%% of span)\n', 'fully solid', Vs, Ds, 100*Ds/R)

fprintf('\nHollow box Pareto front\n')
fprintf(' tip defl (m)  volume (cm^3)   xs     wRoot  wTip (mm)   tRoot  tTip (mm)\n')
fprintf('   %6.3f        %7.1f      %5.3f   %6.1f %6.1f      %5.2f  %5.2f\n', [C(2).D; C(2).V; C(2).X'])

%check the lightest design of every concept against the full solver
fprintf('\nCheck lightest designs against beamBending:\n')
for k = 1:numel(C)
    I = conceptSection(C(k).key,C(k).X(end,:),r,c);
    delta = beamBending(N,I);
    fprintf('  %-11s fast %.4f m, solver %.4f m\n', C(k).name, C(k).D(end), abs(delta(end)))
end

%% figures
cols = {[0.85 0.33 0.10], [0 0.45 0.74], [0.3 0.3 0.3]};
marks = {'o-','s-','^-'};

%Pareto fronts of all concepts
figure(1)
for k = 1:numel(C)
    semilogy(C(k).D, C(k).V, marks{k}, 'Color',cols{k}, 'MarkerFaceColor',cols{k}, 'LineWidth',1.5), hold on
end
semilogy(Ds, Vs, 'kp', 'MarkerSize',14, 'MarkerFaceColor','y')
xlabel('tip deflection, m'), ylabel('volume, cm^3 (log scale)'), grid on
legend([{C.name} {'fully solid blade'}], 'Location','northeast')
title(sprintf('Optimised concepts, linear tapers, SF = %.1f', SF))

%cross-sections of the lightest design of each concept (root filled, tip dashed)
figure(2)
for k = 1:numel(C)
    subplot(numel(C),1,k)
    drawSection(C(k).key, C(k).X(end,:), c)
    title(sprintf('%s, lightest: V = %.0f cm^3 (root filled, tip dashed)', C(k).name, C(k).V(end)))
end

%how the optimal box changes along its front
figure(3)
subplot(1,3,1), plot(C(2).D, C(2).X(:,1), 's-'), grid on
xlabel('tip deflection, m'), ylabel('box centre x_s / c')
subplot(1,3,2), plot(C(2).D, C(2).X(:,2:3), 's-'), grid on
xlabel('tip deflection, m'), ylabel('box width, mm'), legend('root','tip')
subplot(1,3,3), plot(C(2).D, C(2).X(:,4:5), 's-'), grid on
xlabel('tip deflection, m'), ylabel('wall thickness, mm'), legend('root','tip')
sgtitle('Hollow box: Pareto-optimal design variables')


%run fmincon from every start and keep the best feasible optimum
function xBest = bestOf(fun, starts, con, nonlcon, opts)
    fBest = inf; xBest = starts(1,:);
    for i = 1:size(starts,1)
        [x, f, flag] = fmincon(fun, starts(i,:), con.A, con.b, [], [], con.lb, con.ub, nonlcon, opts);
        if flag > 0 && f < fBest
            xBest = x; fBest = f;
        end
    end
    if isinf(fBest), warning('no feasible optimum found from %d starts', size(starts,1)), end
end

%section properties of one design of one concept at every node
function [I, A, ymax] = conceptSection(key,x,r,c)
    lin = @(a,b) (a + (b - a)*r/r(end))*1e-3; %linear taper, mm -> m
    switch key
        case 'spar' %solid spar = box whose walls fill it (t = w/2)
            w = lin(x(2),x(3));
            [I, A, ymax] = boxSection(w,w/2,x(1),c);
        case 'box'
            [I, A, ymax] = boxSection(lin(x(2),x(3)),lin(x(4),x(5)),x(1),c);
        case 'skin'
            [I, A, ymax] = skinSection(lin(x(1),x(2)),c);
    end
end

%evaluate one design without calling beamBending:
%volume = integral of area along the span
%tip deflection from the unit-load method: integral M(r)(R - r)/(E I(r)) dr
%stress = |M| y_max / I at every node
function [V, dTip, sigMax] = evalConcept(key,x,r,c,M,E)
    [I, A, ymax] = conceptSection(key,x,r,c);
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

%hollow box spar (same as RunBeamBending.m): outer strips minus inner strips
function [I, A, ymax] = boxSection(w,t,xs,c)
    nStrips = 50;
    w = w(:); t = t(:);
    x1 = max(xs*c - w/2, 0);
    x2 = min(xs*c + w/2, c);
    [Io, Ao, ymax] = stripSum(x1, x2, 0*t, c, nStrips);
    [Ii, Ai] = stripSum(x1 + t, max(x2 - t, x1 + t), t, c, nStrips);
    I = (Io - Ii)'; A = (Ao - Ai)'; ymax = ymax';
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

%thin structural skin (same as RunBeamBending.m and ParetoSkin.m)
function [I, A, ymax] = skinSection(t,c)
    nPanels = 200;
    beta = linspace(0,pi,nPanels+1);
    x = c*(1 - cos(beta))/2;
    y = naca0012(x,c);
    Lp = hypot(diff(x),diff(y));
    theta = atan2(diff(y),diff(x));
    d = (y(1:end-1) + y(2:end))/2;
    Iunit = 2*sum(Lp.^3.*sin(theta).^2/12 + Lp.*d.^2);
    perimeter = 2*sum(Lp);
    I = Iunit*t;
    A = perimeter*t;
    ymax = max(y)*ones(size(t));
end

%draw the aerofoil with the section of one design: root filled, tip dashed
function drawSection(key,x,c)
    warning('off','MATLAB:polyshape:repairedBySimplify')
    xo = linspace(0,c,300); yo = naca0012(xo,c);
    outline = polyshape([xo fliplr(xo)], [yo -fliplr(yo)]);
    hold on
    for end_ = 1:2
        if strcmp(key,'skin')
            t = x(end_)*1e-3;
            sec = subtract(outline, polybuffer(outline,-t));
        else
            xs = x(1); w = x(1+end_)*1e-3;
            if strcmp(key,'box'), t = x(3+end_)*1e-3; else, t = w/2; end
            xb = linspace(max(xs*c - w/2,0), min(xs*c + w/2,c), 100); yb = naca0012(xb,c);
            sec = polyshape([xb fliplr(xb)], [yb -fliplr(yb)]);
            if xb(end) - xb(1) > 2*t
                xi = linspace(xb(1) + t, xb(end) - t, 100); yi = max(naca0012(xi,c) - t, 0);
                sec = subtract(sec, polyshape([xi fliplr(xi)], [yi -fliplr(yi)]));
            end
        end
        sec = scale(sec, 1e3); %m -> mm for plotting
        if end_ == 1
            plot(sec,'FaceColor',[0.7 0.8 1],'EdgeColor','b','FaceAlpha',1)
        else
            plot(sec,'FaceColor','none','EdgeColor','r','LineStyle','--','LineWidth',1)
        end
    end
    plot(xo*1e3, yo*1e3, 'k-', xo*1e3, -yo*1e3, 'k-')
    axis equal, grid on, xlim([0 c*1e3])
    xlabel('chordwise position x, mm'), ylabel('y, mm')
end
