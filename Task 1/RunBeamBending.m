clc, clear all, close all

N = 100; %the beam solver discretises beam into 100 nodes
c = 0.15; %define the chord length - KEEP FIXED
R = 1.5; %blade length (m)
r = linspace(0,R,N); %radial position of every node, root (r = 0) to tip (r = R)

structure = 'skin'; %structural concept to analyse: 'skin' (Approach 2) or 'spar' (Approach 1)

switch structure
    case 'skin'
        %linear taper of skin thickness from root to tip
        tRoot = 1.2e-3; %skin thickness at the root (m)
        tTip = 0.4e-3; %skin thickness at the tip (m) - keep > 0 so I never reaches zero
        t = tRoot + (tTip - tRoot)*r/R; %set the skin thickness at every node (m)
        [I, A, ymax] = skinSection(t,c); %second moment of area, cross-section area and outer fibre distance

    case 'spar'
        %linear taper of spar width from root to tip
        xs = 0.3; %chordwise position of the spar centre, x/c (0.3 = max thickness)
        wRoot = 35e-3; %spar width at the root (m)
        wTip = 5e-3; %spar width at the tip (m) - keep > 0 so I never reaches zero
        w = wRoot + (wTip - wRoot)*r/R; %set the spar width at every node (m)
        [I, A, ymax] = sparSection(w,xs,c); %second moment of area, cross-section area and outer fibre distance
end

%submit to beamBending:
% number of nodes for numerical solve, N
% vector of second moment of area values across the span, I
% (I must have N elements)
[delta, M, V, L, P]  = beamBending(N,I);
%beamBending function returns:
%beam deflection, delta
%bending moment distribution, M
%shear force distriubition, V
%blade length is returned, L
%distributed aerodynamic loading, P

%define the spanwise position of every node (in meters)
nodeLocations = linspace(0,L,N);

%volume of the structure = integral of cross-section area along the span
vol = trapz(nodeLocations,A); %m^3
sigma = abs(M).*ymax./I; %bending stress magnitude at every node (Pa) - solver returns M as negative
fprintf('Structure: %s\n', structure)
fprintf('Volume = %.1f cm^3 (mass = %.3f kg at ABS ~1050 kg/m^3)\n', vol*1e6, vol*1050)
fprintf('Tip deflection = %.3f m, max bending stress = %.2f MPa (yield 30 MPa)\n', abs(delta(end)), max(sigma)/1e6)

%plot the aerofoil and the structural cross-section
figure(1)
switch structure
    case 'skin', plotSkin(t,c)
    case 'spar', plotSpar(w,xs,c)
end

%plot the results
figure(2)
subplot(3,1,1) %3 subplots in one figure (check the docs for more info)
plot(nodeLocations,delta,'ko-') %plot the deflection (delta, returned in m) at each node location across the blade
xlabel('distance, m') %set labels
ylabel('deflection, m')

subplot(3,1,2) %plot on the second subplot
plot(nodeLocations,sigma/1e6,'ko-'), hold on %plot the bending stress (M*y / I) at every node location
yline(30,'r--','yield') %ABS yield strength
xlabel('distance, m') %set labels
ylabel('Bending Stress, MPa')

subplot(3,1,3) %plot the cross-section area along the span (area under this curve = volume)
area(nodeLocations,A*1e6,'FaceColor',[0.7 0.8 1])
xlabel('distance, m')
ylabel('area, mm^2')
title(sprintf('Volume = %.1f cm^3',vol*1e6))


%NACA0012 half-thickness y_t at chordwise position x (m)
function yt = naca0012(x,c)
    xc = x/c;
    yt = c*5*0.12*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 ...
                   + 0.2843*xc.^3 - 0.1015*xc.^4);
end

%Approach 2: thin structural skin of thickness t following the NACA0012
%surface (t << c). The upper surface is split into straight panels; each
%panel has length Lp, angle theta and centroid height d, so
%   I_panel = Lp^3 t sin^2(theta)/12 + (Lp t) d^2   (parallel axis theorem)
%The lower surface is a mirror image, so the totals are doubled.
%Area = skin thickness x perimeter - see theory guide.
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

%Approach 1: single solid spar of width w centred at x/c = xs, filling the
%full aerofoil depth (top and bottom follow the NACA0012 surface).
%The spar is split into thin vertical rectangles (strips), each of height
%2*y_t(x), and I = sum(b*d^3/12), A = sum(b*d) - see theory guide.
function [I, A, ymax] = sparSection(w,xs,c)
    nStrips = 50; %number of rectangles used to approximate the spar
    I = zeros(size(w)); A = zeros(size(w)); ymax = zeros(size(w));
    for k = 1:numel(w)
        x1 = max(xs*c - w(k)/2, 0); %spar front edge (clipped to the aerofoil)
        x2 = min(xs*c + w(k)/2, c); %spar rear edge
        edges = linspace(x1,x2,nStrips+1);
        b = diff(edges); %strip widths
        d = 2*naca0012((edges(1:end-1)+edges(2:end))/2,c); %strip heights (midpoint)
        I(k) = sum(b.*d.^3/12);
        A(k) = sum(b.*d);
        ymax(k) = max(naca0012(linspace(x1,x2,200),c)); %furthest fibre from neutral axis
    end
end

%plot the NACA0012 outline with the skin at the root (blue) and the inner
%surface of the tip skin (red dashed). Skin drawn inwards from the outer surface.
function plotSkin(t,c)
    beta = linspace(0,pi,400);
    x = c*(1 - cos(beta))/2;
    y = naca0012(x,c);
    nx = -gradient(y)./hypot(gradient(x),gradient(y)); %inward unit normal of the upper surface
    ny = -gradient(x)./hypot(gradient(x),gradient(y));
    nx(1) = 1; ny(1) = 0; %leading edge points straight back
    xo = [x fliplr(x)]; yo = [y -fliplr(y)]; %outer surface, upper then lower
    inner = @(tk) deal([x+tk*nx fliplr(x+tk*nx)], [max(y+tk*ny,0) -fliplr(max(y+tk*ny,0))]);
    [xr,yr] = inner(t(1)); [xt,yt] = inner(t(end));
    fill([xo fliplr(xr)]*1e3, [yo fliplr(yr)]*1e3, [0.7 0.8 1], 'EdgeColor','b'), hold on
    plot(xt*1e3, yt*1e3, 'r--', 'LineWidth',1.2)
    plot(xo*1e3, yo*1e3, 'k-')
    axis equal, grid on
    xlabel('chordwise position x, mm'), ylabel('y, mm')
    title(sprintf('NACA0012 structural skin: t = %.2f mm (root, blue), %.2f mm (tip, red dashed)', ...
          t(1)*1e3, t(end)*1e3))
end

%plot the NACA0012 outline with the spar at the root (blue) and tip (red dashed)
function plotSpar(w,xs,c)
    x = linspace(0,c,400);
    yt = naca0012(x,c);
    plot(x*1e3,yt*1e3,'k-',x*1e3,-yt*1e3,'k-'), hold on
    spars = [w(1) w(end)]; cols = {[0.7 0.8 1],[1 1 1]}; edgeCols = {'b','r'}; styles = {'-','--'};
    for k = 1:2
        xsp = linspace(max(xs*c-spars(k)/2,0), min(xs*c+spars(k)/2,c), 100);
        ysp = naca0012(xsp,c);
        fill([xsp fliplr(xsp)]*1e3, [ysp -fliplr(ysp)]*1e3, cols{k}, 'FaceAlpha',double(k==1), ...
             'EdgeColor',edgeCols{k},'LineStyle',styles{k},'LineWidth',1.2)
    end
    axis equal, grid on
    xlabel('chordwise position x, mm'), ylabel('y, mm')
    title(sprintf('NACA0012 with single spar: w = %.1f mm (root, blue), %.1f mm (tip, red dashed)', ...
          w(1)*1e3, w(end)*1e3))
end
