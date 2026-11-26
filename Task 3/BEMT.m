function out = BEMT(B, radiusLimits, lambda, chord, alphaIn, N, V, nu, tol, maxIter, opts)
% BEMT  Copy of the supplied BEMT_Demo with three additions, all switchable.
%
% Same physics, inputs and iteration as BEMT_Demo when opts.tipLoss = false
% and opts.alphaMode = 'fixed' (checked in ValidateBEMT.m). Additions:
%   1. Aerofoil lookup built once (S822Lookup) instead of griddata every
%      iteration -> same numbers, much faster.
%   2. Optional Prandtl tip- and hub-loss factor F in the momentum equations,
%      with Buhl's high-induction correction (on by default with tip loss),
%      because F -> 0 near the tip pushes a above ~0.4 where simple
%      momentum theory breaks down.
%   3. alphaMode 'bestLD': each section's angle of attack is its best-L/D
%      angle for its own Reynolds number, plus an offset (alphaIn) that the
%      optimiser can change.
%      alphaMode 'twist': alphaIn is the blade twist (deg), so the AoA
%      follows the flow. Used to test a finished blade away from its design
%      tip-speed ratio (Cp-lambda curve).
% Chord and alphaIn can be 2 values (root, tip -> linear, as in the demo)
% or N values (one per element, e.g. from a spline).
% Returns a struct with Cp, CT and the spanwise distributions.

if nargin < 11, opts = struct; end
if ~isfield(opts, 'tipLoss'),   opts.tipLoss   = false;   end
if ~isfield(opts, 'hubLoss'),   opts.hubLoss   = false;   end   % hub adds little power; F->0 there destabilises BEM
if ~isfield(opts, 'alphaMode'), opts.alphaMode = 'fixed'; end
if ~isfield(opts, 'highInduction'), opts.highInduction = opts.tipLoss; end
if ~isfield(opts, 'lookup'),    opts.lookup    = S822Lookup(); end
if ~isfield(opts, 'relax'),     opts.relax     = 1;       end   % <1 damps oscillating iterations
relax = opts.relax;
Lk  = opts.lookup;
rho = 1.225;                                   % as hard-coded in BEMT_Demo
bestLD = strcmp(opts.alphaMode, 'bestLD');
twistM = strcmp(opts.alphaMode, 'twist');      % alphaIn = twist angle (deg): a built blade

r  = linspace(radiusLimits(1), radiusLimits(2), N);
R  = r(end);  Rh = r(1);
omega = V*lambda/R;
c     = spread(chord, N);
alpIn = spread(alphaIn, N);                    % degrees

[Ft, Fn, aS, aPS, phiS, gamS, ReS, CLS, CDS, alS, FS] = deal(zeros(1, N));
conv = true(1, N);

for i = 1:N
    sigma = B*c(i)/(2*pi*r(i));                % local solidity
    if i == 1, a = 0.33; aP = 0; end           % same start as BEMT_Demo
    dA = inf; dAP = inf; it = 0; noLoad = false;
    while (abs(dA) > tol || abs(dAP) > tol) && it < maxIter
        phi = atan((1-a)/(lambda*(r(i)/R)*(1+aP)));     % inflow angle
        W   = V*(1-a)/sin(phi);                          % relative wind
        Re  = W*c(i)/nu;
        ReQ = min(max(Re, Lk.ReRange(1)*(1+1e-12)), Lk.ReRange(2)*(1-1e-12));   % demo saturates Re too
        if bestLD,     al = Lk.alphaBest(ReQ) + alpIn(i);
        elseif twistM, al = rad2deg(phi) - alpIn(i);      % fixed twist: AoA follows the flow
        else,          al = alpIn(i);
        end
        alQ = min(max(al, Lk.alphaRange(1)+1e-9), Lk.alphaRange(2)-1e-9);   % stay just inside the data edge
        CL = Lk.CL(alQ, ReQ);  CD = Lk.CD(alQ, ReQ);
        Cn = CL*cos(phi) + CD*sin(phi);
        Ct = CL*sin(phi) - CD*cos(phi);

        F = 1;                                            % Prandtl loss factor
        sp = max(sin(phi), 1e-6);                         % guard: phi <= 0 only in bad iterates
        if opts.tipLoss
            F = F*(2/pi)*acos(exp(-B*(R - r(i))/(2*r(i)*sp)));
        end
        if opts.hubLoss
            F = F*(2/pi)*acos(exp(-B*(r(i) - Rh)/(2*Rh*sp)));
        end
        if F < 1e-9, noLoad = true; break; end   % very tip/hub: no load carried

        aNew  = 1/(4*F*sin(phi)^2/(sigma*Cn) + 1);
        if opts.highInduction                              % Buhl (2005) correction:
            CTl = sigma*(1-a)^2*Cn/sin(phi)^2;               % momentum theory fails for a > ~0.4
            if CTl > 0.96*F
                aNew = (18*F - 20 - 3*sqrt(max(CTl*(50 - 36*F) + 12*F*(3*F - 4), 0)))/(36*F - 50);
            end
        end
        aPNew = 1/(4*F*sin(phi)*cos(phi)/(sigma*Ct) - 1);
        if opts.highInduction, aNew = min(aNew, 0.95); end   % keep 1-a > 0 (flow not reversed)
        if ~isfinite(aNew) || ~isfinite(aPNew), it = maxIter; break; end   % failed element: flag, don't spread NaN
        dA  = aNew - a;   dAP = aPNew - aP;
        a   = a + relax*dA;   aP = aP + relax*dAP;      % relax = 1: same as the demo
        it  = it + 1;
    end
    conv(i) = noLoad || it < maxIter;
    if noLoad
        Ft(i) = 0; Fn(i) = 0;
    else
        Ft(i) = 0.5*Ct*rho*c(i)*W^2;            % drives the rotor
        Fn(i) = 0.5*Cn*rho*c(i)*W^2;            % bends the blade (thrust)
    end
    aS(i) = a; aPS(i) = aP; phiS(i) = phi; gamS(i) = phi - deg2rad(al);
    ReS(i) = Re; CLS(i) = CL; CDS(i) = CD; alS(i) = al; FS(i) = F;
end

Q  = B*trapz(r, r.*Ft);                        % torque
P  = omega*Q;
T  = B*trapz(r, Fn);                           % thrust (the demo's CT line used torque)
A  = pi*R^2;
out = struct('Cp', P/(0.5*rho*A*V^3), 'CT', T/(0.5*rho*A*V^2), 'P', P, 'Q', Q, 'T', T, ...
    'r', r, 'c', c, 'alpha', alS, 'twist', rad2deg(gamS), 'phi', rad2deg(phiS), ...
    'Re', ReS, 'a', aS, 'aPrime', aPS, 'CL', CLS, 'CD', CDS, 'F', FS, ...
    'Ft', Ft, 'Fn', Fn, 'converged', conv, 'omega', omega);
end

function v = spread(x, N)
if numel(x) == 2, v = linspace(x(1), x(2), N);
elseif numel(x) == N, v = x(:)';
else, error('BEMT:size', 'Profile must have 2 or N values');
end
end
