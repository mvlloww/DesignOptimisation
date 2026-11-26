function OptimisePeakCp(quick)
% OPTIMISEPEAKCP  Option 2: blade with the highest peak Cp, for B = 1..5.
%
%   maximise  Cp(lambda, chord, AoA offsets)       (variable-speed turbine:
%                                                   lambda is a design variable)
%   at        V = 7.63 m/s  (cube-root-mean wind of a 6.15 m/s Rayleigh mean)
%   subject to
%     tip speed  lambda*V <= 65 m/s          (noise)
%     20 mm <= chord <= 250 mm, chord never grows toward the tip
%     every section inside the S822 data: 20k <= Re <= 850k, 0 <= AoA <= 12 deg
%     every BEMT element converged
%   physics: BEMT.m with Prandtl tip loss + Buhl high-induction correction ON
%
% AoA: each section starts at its best-L/D angle for its own Re
% (S822Lookup); the optimiser adds a linear root->tip offset.
% Stage 1: linear chord (root, tip)              -> 5 variables, 20 starts
% Stage 2: pchip chord through 5 control points  -> 8 variables, warm start
%          from stage 1 + 7 random starts
% Benchmark: Schmitz ideal chord at the same lambda.
% Afterwards PostProcessPeakCp.m checks the finished blades and makes figures.

here = fileparts(mfilename('fullpath'));
L = S822Lookup(fullfile(here, 'Wind Turbine Blade Optimisation Demo'));
figDir = fullfile(here, 'figures'); if ~exist(figDir, 'dir'), mkdir(figDir); end

% ---------------- fixed settings ----------------
P.Vmean  = 6.153;                         % Met Office winter median of the mean
P.V      = (6/pi)^(1/3)*P.Vmean;          % = 7.63 m/s, same mean power as Rayleigh
P.R      = 1.1;  P.Rhub = 0.15;           % D = 2.2 m limit; hub ~14% of R
P.nu     = 1.5e-5;                        % as in the demo
P.N      = 60;   P.tol = 1e-10; P.maxIt = 1000;   % N=60: within 0.4% of N=240
P.tipMax = 65;   P.lamMax = P.tipMax/P.V; % = 8.5
P.cMin   = 0.02; P.cMax = 0.25;
P.dAmax  = 5;                             % AoA offset range +-5 deg
P.Bs     = 1:5;
P.nStart1 = 20;  P.nStart2 = 8;
P.rk     = linspace(P.Rhub, P.R, 5);      % stage-2 chord control stations
opts     = struct('lookup', L, 'tipLoss', true, 'alphaMode', 'bestLD');
if nargin > 0 && quick, P.Bs = 3; P.nStart1 = 2; P.nStart2 = 2; figDir = tempdir; end   % smoke test
fo = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 300, ...
    'OptimalityTolerance', 1e-7, 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-6, 'FiniteDifferenceStepSize', 1e-6);

fprintf('Design wind %.2f m/s, lambda <= %.2f, %d elements, tip loss ON\n\n', P.V, P.lamMax, P.N);
res = struct([]);
for B = P.Bs
    tB = tic;
    rng(B);                                            % repeatable random starts
    % ---------- Stage 1: x = [lambda, c_root, c_tip, dA_root, dA_tip]
    lb1 = [2, P.cMin, P.cMin, -P.dAmax, -P.dAmax];
    ub1 = [P.lamMax, P.cMax, P.cMax, P.dAmax, P.dAmax];
    A1  = [0 -1 1 0 0]; b1 = 0;                        % c_tip <= c_root
    geo1 = @(x) deal(x(2:3), x(4:5));
    Z0 = rand(P.nStart1, 5); Z0(:,2:3) = sort(Z0(:,2:3), 2, 'descend');
    [best1, log1] = multistart(@(x) design(B, x(1), geo1, x), lb1, ub1, A1, b1, Z0);

    % ---------- Stage 2: x = [lambda, c1..c5, dA_root, dA_tip]
    lb2 = [2, P.cMin*ones(1,5), -P.dAmax, -P.dAmax];
    ub2 = [P.lamMax, P.cMax*ones(1,5), P.dAmax, P.dAmax];
    A2  = zeros(4, 8); for k = 1:4, A2(k, 1+k) = -1; A2(k, 2+k) = 1; end; b2 = zeros(4,1);
    geo2 = @(x) deal(max(pchip(P.rk, x(2:6), linspace(P.Rhub, P.R, P.N)), 1e-3), x(7:8));
    x1  = best1.x;                                       % warm start: stage-1 blade
    z1  = ([x1(1), interp1([P.Rhub P.R], x1(2:3), P.rk), x1(4:5)] - lb2)./(ub2 - lb2);
    Z0  = [z1; rand(P.nStart2-1, 8)]; Z0(2:end,2:6) = sort(Z0(2:end,2:6), 2, 'descend');
    [best2, log2] = multistart(@(x) design(B, x(1), geo2, x), lb2, ub2, A2, b2, Z0);

    % ---------- Benchmark: Schmitz ideal chord at the stage-2 lambda
    lam = best2.x(1); r = linspace(P.Rhub, P.R, P.N);
    ReD = 2e5; CLd = L.CL(L.alphaBest(ReD), ReD);
    cS  = 16*pi*r/(B*CLd).*sin(atan(P.R./(lam*r))/3).^2;
    cS  = min(max(cS, P.cMin), P.cMax);
    oS  = BEMT(B, [P.Rhub P.R], lam, cS, [0 0], P.N, P.V, P.nu, P.tol, P.maxIt, opts);

    res(B).B = B; res(B).s1 = best1; res(B).s2 = best2; res(B).log1 = log1; res(B).log2 = log2;
    res(B).schmitz = oS;
    fprintf('B=%d  Schmitz %.4f | stage1 %.4f (lam %.2f, %d/%d starts within 0.1%%) | stage2 %.4f (lam %.2f, %d/%d)\n', ...
        B, oS.Cp, best1.Cp, best1.x(1), log1.nNear, P.nStart1, best2.Cp, lam, log2.nNear, P.nStart2);
    fprintf('      time %.0f s, function evaluations %d + %d\n', toc(tB), sum(log1.nFun), sum(log2.nFun));
end
if ~(nargin > 0 && quick), save(fullfile(here, 'results_peakCp.mat'), 'res', 'P'); end
PostProcessPeakCp(res, P, L, figDir);   % checks on the finished blades + figures

    % ---------------- nested helpers ----------------
    function o = design(B, lam, geo, x)
        [ch, dA] = geo(x);
        o = BEMT(B, [P.Rhub P.R], lam, ch, dA, P.N, P.V, P.nu, P.tol, P.maxIt, opts);
    end

    function [best, lg] = multistart(evalX, lb, ub, A, b, Z0)
        % optimise in normalised z in [0,1] so every variable is order one
        toX = @(z) lb + z(:)'.*(ub - lb);
        zC = []; oC = [];
        function o = cached(z)
            if isempty(zC) || any(z(:)' ~= zC), zC = z(:)'; oC = evalX(toX(z)); end
            o = oC;
        end
        function f = obj(z)
            o = cached(z); f = -o.Cp; if ~isfinite(f), f = 1; end
        end
        function [c, ceq] = con(z)
            o = cached(z); m = real(o.F) > 1e-9 & isfinite(o.Re) & isfinite(o.alpha);
            if ~any(m), c = ones(5,1); ceq = []; return; end
            c = [(L.ReRange(1) - min(o.Re(m)))/L.ReRange(1);
                 (max(o.Re(m)) - L.ReRange(2))/L.ReRange(2);
                 (L.alphaRange(1) - min(o.alpha(m)))/12;
                 (max(o.alpha(m)) - L.alphaRange(2))/12;
                 0.5 - double(all(o.converged))];
            c(~isfinite(c)) = 1; ceq = [];
        end
        n = size(Z0, 1); cp = -inf(n,1); xs = zeros(n, numel(lb)); nf = zeros(n,1); ok = false(n,1);
        tic;
        for s = 1:n
            [z, ~, flag, outp] = fmincon(@obj, Z0(s,:), A, b, [], [], zeros(size(lb)), ones(size(lb)), @con, fo);
            [cc, ~] = con(z); ok(s) = flag > 0 && all(cc <= 1e-6);
            o = cached(z); xs(s,:) = toX(z); nf(s) = outp.funcCount;
            if ok(s), cp(s) = o.Cp; end
        end
        [~, j] = max(cp);
        best.x = xs(j,:); best.out = evalX(best.x); best.Cp = best.out.Cp;
        lg = struct('cp', cp, 'x', xs, 'nFun', nf, 'feasible', ok, 'time', toc, ...
                    'nNear', sum(cp >= best.Cp*(1 - 1e-3)));
    end
end
