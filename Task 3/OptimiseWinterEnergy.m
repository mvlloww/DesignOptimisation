function OptimiseWinterEnergy(action, worker, nWorkers, mode, Prated)
% OPTIMISEWINTERENERGY  Option 4: blade + rotor speed for maximum winter energy (B = 3).
%
%   maximise  E = sum_bins P(V) x Rayleigh probability x 2160 h     (kWh, Dec-Feb)
%   wind      Rayleigh, mean 6.15 m/s, 1 m/s bins, cut-in 3, cut-out 20 m/s
%   control   FIXED speed: omega constant, so lambda = tip/V changes with wind
%   variables x = [tip speed, Vd, c_root, c_tip, dA_root, dA_tip]
%             twist = best-L/D twist at design wind Vd (lambda_d = tip/Vd),
%             i.e. the peak-Cp construction, then frozen and run at all winds
%   subject to
%     tip speed <= 65 m/s  (42 dB(A) at a house 100 m away, MCS 020)
%     chord 20-250 mm, tapering; design point inside S822 data; BEMT
%     converged at the design point and in every wind bin
%   physics   BEMT.m, tip loss + Buhl, Viterna post-stall above 12 deg
%
% Compared against
%   - VARIABLE speed with the B = 3 peak-Cp blade (tracks best lambda until
%     the tip-speed limit, then constant omega): the upper bound
%   - the B = 3 peak-Cp blade run at its best FIXED speed (1-D search)
%
% No Parallel Computing Toolbox, so starts are shared between separate
% MATLAB processes:
%   OptimiseWinterEnergy('worker', k, n)   -> starts k, k+n, ... -> winter_part_k.mat
%   OptimiseWinterEnergy('combine')        -> best design, comparisons, checks, figures
%   OptimiseWinterEnergy('quick')          -> 1 short start, smoke test
%
% VARIABLE-speed version (mode = 'variable'): same blade construction, but
%   x(1) = lambda_opt, the design tip-speed ratio, which is also the
%   controller's target: lambda = min(lambda_opt, tipMax/V). So the blade is
%   optimised for energy WITH variable-speed control, for a fair comparison.
%   OptimiseWinterEnergy('worker', k, n, 'variable') -> winterVar_part_k.mat
%   OptimiseWinterEnergy('combineVar')               -> adds it to results_winter.mat
%
% RATED-POWER CAP (Prated, W): shaft power capped at Prated in every bin.
%   OptimiseWinterEnergy('rescore')                    -> existing 4 blades at 1, 1.5, 2 kW
%   OptimiseWinterEnergy('worker', k, n, mode, 1500)   -> winterCap1500_part_k / winterVarCap1500_part_k
%   OptimiseWinterEnergy('combineCap', [], [], [], 1500) -> results_winter_cap.mat + figure

if nargin < 1, action = 'quick'; end
if nargin < 4 || isempty(mode), mode = 'fixed'; end
if nargin < 5, Prated = inf; end
here = fileparts(mfilename('fullpath'));
[P, L, W, B, peak] = setup(here);
P.Prated = Prated;

switch action
    case {'worker', 'quick'}
        if strcmp(action, 'quick'), worker = 1; nWorkers = 1; P.nStart = 1; P.fmIt = 3; end
        [lb, ub] = bounds(P, mode);
        rng(3); Z0 = rand(P.nStart, 6); Z0(:,3:4) = sort(Z0(:,3:4), 2, 'descend');
        x1 = peak.s1.x;                                   % start 1: peak-Cp linear blade
        if strcmp(mode, 'fixed'), x0 = [P.tipMax, P.tipMax/x1(1), x1(2:5)];   % at max tip speed
        else,                     x0 = [x1(1), P.V, x1(2:5)];                % at its own design point
        end
        Z0(1,:) = (x0 - lb)./(ub - lb);
        Z0 = min(max(Z0, 0), 1);
        fo = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
            'MaxFunctionEvaluations', 1500, 'MaxIterations', P.fmIt, ...
            'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8, ...
            'ConstraintTolerance', 1e-6, 'FiniteDifferenceStepSize', 1e-5);
        mine = worker:nWorkers:P.nStart;
        runs = struct('start', {}, 'x', {}, 'E', {}, 'feasible', {}, 'nFun', {}, 'time', {}, 'flag', {});
        for s = mine
            t0 = tic;
            [x, E, ok, nf, flag] = oneStart(Z0(s,:), lb, ub, fo, P, L, W, B, mode);
            runs(end+1) = struct('start', s, 'x', x, 'E', E, 'feasible', ok, 'nFun', nf, 'time', toc(t0), 'flag', flag); %#ok<AGROW>
            fprintf('[%s] worker %d start %d: E = %.1f kWh, feasible %d, x = %s, %.0f s\n', ...
                mode, worker, s, E, ok, mat2str(x, 4), toc(t0));
            if ~strcmp(action, 'quick'), save(fullfile(here, sprintf('%s_part_%d.mat', partName(mode, P.Prated), worker)), 'runs'); end
        end
    case 'combine'
        combine(here, P, L, W, B, peak);
    case 'combineVar'
        combineVar(here, P, L, W, B, peak);
    case 'rescore'
        rescore(here, P, L, W, B, peak);
    case 'combineCap'
        combineCap(here, P, L, W, B, peak);
end
end

function n = partName(mode, Prated)
if strcmp(mode, 'fixed'), n = 'winter'; else, n = 'winterVar'; end
if nargin > 1 && isfinite(Prated), n = sprintf('%sCap%d', n, round(Prated)); end
end

function [runs, xb, nNear] = loadRuns(here, name)
runs = [];
for f = dir(fullfile(here, [name '_part_*.mat']))'
    S = load(fullfile(here, f.name)); runs = [runs, S.runs]; %#ok<AGROW>
end
[~, i] = sort([runs.start]); runs = runs(i);
Ef = [runs.E]; Ef(~[runs.feasible]) = -inf;
[Ebest, j] = max(Ef); xb = runs(j).x; nNear = sum(Ef >= Ebest*(1 - 1e-3));
end

function [E, tip] = bestFixedTip(B, c, tw, W, P, L)
% best constant tip speed for a given blade (1-D search)
f = @(t) -WinterEnergy(B, c, tw, struct('mode', 'fixed', 'tip', t), W, P, L);
tip = fminbnd(f, 20, P.tipMax, optimset('TolX', 0.01)); E = -f(tip);
end

function rescore(here, P, L, W, B, peak)
% The four uncapped designs re-scored with a rated-power cap.
% Peak-Cp blade at fixed speed: its best fixed tip speed is re-found for each cap.
S = load(fullfile(here, 'results_winter.mat')); w = S.win; v = S.winVar; pk = peak.s2;
caps = [1000 1500 2000 inf]; E = zeros(4, numel(caps)); tipP = zeros(1, numel(caps));
for k = 1:numel(caps)
    P.Prated = caps(k);
    E(1,k) = WinterEnergy(B, pk.out.c, pk.out.twist, control(pk.x(1), 'variable'), W, P, L);
    E(2,k) = WinterEnergy(B, v.cVo, v.twVo, control(v.xb, 'variable'), W, P, L);
    [E(3,k), tipP(k)] = bestFixedTip(B, pk.out.c, pk.out.twist, W, P, L);
    E(4,k) = WinterEnergy(B, w.cE, w.twE, control(w.xb, 'fixed'), W, P, L);
end
names = {'variable, peak-Cp blade', 'variable, energy-opt (uncapped)', 'fixed, peak-Cp blade (best tip)', 'fixed, energy-opt (uncapped, tip 52)'};
fprintf('\nRe-scored with rated-power cap (kWh, Dec-Feb, mean %.2f m/s):\n', W.Vmean);
fprintf('%-38s %8s %8s %8s %8s\n', 'design', '1 kW', '1.5 kW', '2 kW', 'no cap');
for i = 1:4, fprintf('%-38s %8.1f %8.1f %8.1f %8.1f\n', names{i}, E(i,:)); end
fprintf('peak-Cp blade best fixed tip speed: %.1f / %.1f / %.1f / %.1f m/s\n', tipP);
fprintf('energy-opt gain over peak-Cp: variable %s %%, fixed %s %%\n', ...
    mat2str(round(1000*(E(2,:)./E(1,:) - 1))/10), mat2str(round(1000*(E(4,:)./E(3,:) - 1))/10));
save(fullfile(here, 'results_rescore_cap.mat'), 'caps', 'E', 'tipP', 'names');
end

function combineCap(here, P, L, W, B, peak)
% Designs re-optimised WITH the rated-power cap P.Prated, then tested at other caps.
Pr = P.Prated; pk = peak.s2;
[rF, xF, nF] = loadRuns(here, partName('fixed', Pr));
[rV, xV, nV] = loadRuns(here, partName('variable', Pr));
for m = {{'FIXED', rF, nF, 'tip  '}, {'VARIABLE', rV, nV, 'lamOp'}}
    m = m{1}; runs = m{2};
    fprintf('\n%s speed, cap %.0f W: %d starts, %d feasible, %d within 0.1%% of best\n', m{1}, Pr, ...
        numel(runs), sum([runs.feasible]), m{3});
    fprintf(' start   E(kWh)  feasible  %s  Vd    c_root c_tip  dA_root dA_tip  nFun  time(s)\n', m{4});
    for r = runs
        fprintf('%4d   %7.1f   %d      %5.2f %5.2f  %5.3f  %5.3f  %6.2f  %6.2f  %5d %6.0f\n', ...
            r.start, r.E, r.feasible, r.x, r.nFun, r.time);
    end
end
[cF, twF, odF] = buildBlade(designPoint(xF, 'fixed'), P, L, B);
[cV, twV, odV] = buildBlade(designPoint(xV, 'variable'), P, L, B);
caps = [1000 1500 2000 inf]; E = zeros(4, numel(caps)); tipP = zeros(1, numel(caps)); tipF = tipP;
d = cell(1, 4);
for k = 1:numel(caps)
    P.Prated = caps(k);
    [E(1,k), dd{1}] = WinterEnergy(B, pk.out.c, pk.out.twist, control(pk.x(1), 'variable'), W, P, L);
    [E(2,k), dd{2}] = WinterEnergy(B, cV, twV, control(xV, 'variable'), W, P, L);
    [E(3,k), tipP(k)] = bestFixedTip(B, pk.out.c, pk.out.twist, W, P, L);
    [~, dd{3}] = WinterEnergy(B, pk.out.c, pk.out.twist, control(tipP(k), 'fixed'), W, P, L);
    [E(4,k), tipF(k)] = bestFixedTip(B, cF, twF, W, P, L);      % re-tune tip speed at each cap
    [~, dd{4}] = WinterEnergy(B, cF, twF, control(tipF(k), 'fixed'), W, P, L);
    if caps(k) == Pr, d = dd; end
end
P.Prated = Pr;
names = {'variable, peak-Cp blade', sprintf('variable, optimised at %.1f kW cap', Pr/1000), ...
         'fixed, peak-Cp blade (best tip)', sprintf('fixed, optimised at %.1f kW cap (best tip)', Pr/1000)};
fprintf('\nWinter energy with rated-power cap (kWh, Dec-Feb, mean %.2f m/s):\n', W.Vmean);
fprintf('%-42s %8s %8s %8s %8s\n', 'design', '1 kW', '1.5 kW', '2 kW', 'no cap');
for i = 1:4, fprintf('%-42s %8.1f %8.1f %8.1f %8.1f\n', names{i}, E(i,:)); end
fprintf('optimised gain over peak-Cp: variable %s %%, fixed %s %%\n', ...
    mat2str(round(1000*(E(2,:)./E(1,:) - 1))/10), mat2str(round(1000*(E(4,:)./E(3,:) - 1))/10));
fprintf('variable vs fixed (both optimised): %s %%\n', mat2str(round(1000*(E(2,:)./E(4,:) - 1))/10));
fprintf('best fixed tip speed, peak-Cp blade: %s m/s; optimised blade: %s m/s\n', mat2str(round(10*tipP)/10), mat2str(round(10*tipF)/10));
fprintf('optimised variable: lamOpt %.2f, Vd %.2f, chord %.0f -> %.0f mm, design Cp %.4f\n', xV(1), xV(2), 1000*cV(1), 1000*cV(end), odV.Cp);
fprintf('optimised fixed:    tip %.1f, Vd %.2f, chord %.0f -> %.0f mm, design Cp %.4f\n', xF(1), xF(2), 1000*cF(1), 1000*cF(end), odF.Cp);
fprintf('peak-Cp blade: chord %.0f -> %.0f mm, Cp %.4f\n', 1000*pk.out.c(1), 1000*pk.out.c(end), pk.Cp);
cap = struct('Prated', Pr, 'caps', caps, 'E', E, 'tipP', tipP, 'tipF', tipF, 'names', {names}, 'd', {d}, ...
    'xF', xF, 'xV', xV, 'cF', cF, 'twF', twF, 'cV', cV, 'twV', twV, 'odF', odF, 'odV', odV, 'runsF', rF, 'runsV', rV);
save(fullfile(here, 'results_winter_cap.mat'), 'cap', 'P');
capFigure(cap, pk, W, fullfile(here, 'figures'));
end

function capFigure(c, pk, W, figDir)
cols = lines(4);
f = figure('Visible', 'off', 'Position', [0 0 1100 750]);
tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on;
for k = 1:4, plot(c.d{k}.V, c.d{k}.P, 'o-', 'Color', cols(k,:), 'LineWidth', 1.5, 'MarkerSize', 4); end
yline(c.Prated, 'k:', sprintf('rated %.1f kW', c.Prated/1000));
xlabel('Wind speed (m/s)'); ylabel('Shaft power (W)'); legend(c.names, 'Location', 'southeast', 'FontSize', 8);
title(sprintf('(a) Power curves, %.1f kW cap', c.Prated/1000));
nexttile; hold on; grid on;
hb = bar(categorical({'1 kW', '1.5 kW', '2 kW', 'no cap'}, {'1 kW', '1.5 kW', '2 kW', 'no cap'}), c.E', 1);
for k = 1:4, hb(k).FaceColor = cols(k,:); end
ylabel('Winter energy (kWh)'); title(sprintf('(b) Winter energy vs rated power (blades optimised at %.1f kW)', c.Prated/1000));
nexttile; hold on; grid on;
plot(pk.out.r, 1000*pk.out.c, 'k--', 'LineWidth', 1.3);
plot(c.odV.r, 1000*c.cV, 'Color', cols(2,:), 'LineWidth', 1.5); plot(c.odF.r, 1000*c.cF, 'Color', cols(4,:), 'LineWidth', 1.5);
xlabel('r (m)'); ylabel('Chord (mm)'); legend('peak-C_p', 'optimised, variable', 'optimised, fixed'); title('(c) Chord');
nexttile; hold on; grid on;
plot(pk.out.r, pk.out.twist, 'k--', 'LineWidth', 1.3);
plot(c.odV.r, c.twV, 'Color', cols(2,:), 'LineWidth', 1.5); plot(c.odF.r, c.twF, 'Color', cols(4,:), 'LineWidth', 1.5);
xlabel('r (m)'); ylabel('Twist (deg)'); title('(d) Twist');
saveBoth(f, 't3_winter_cap', figDir);
end

function [lb, ub] = bounds(P, mode)
% fixed:    x(1) = tip speed (m/s), x(2) = design wind Vd, lambda_d = tip/Vd
% variable: x(1) = lambda_opt (design + controller target), x(2) = Vd
if strcmp(mode, 'fixed'), lb1 = 20; ub1 = P.tipMax; else, lb1 = 2; ub1 = 10; end
lb = [lb1, 4,  P.cMin, P.cMin, -P.dAmax, -P.dAmax];
ub = [ub1, 14, P.cMax, P.cMax,  P.dAmax,  P.dAmax];
end

function x = designPoint(x, mode)
% buildBlade wants [tip speed, Vd, ...]; in variable mode x(1) is lambda
if strcmp(mode, 'variable'), x(1) = x(1)*x(2); end
end

function ctrl = control(x, mode)
if strcmp(mode, 'fixed'), ctrl = struct('mode', 'fixed', 'tip', x(1));
else,                     ctrl = struct('mode', 'variable', 'lamOpt', x(1));
end
end

% ------------------------------------------------------------------------
function [P, L, W, B, peak] = setup(here)
S = load(fullfile(here, 'results_peakCp.mat'));
P = S.P;  B = 3;  peak = S.res(B);
P.tipMax = 65;  P.nStart = 12;  P.fmIt = 100;
L = S822Lookup(fullfile(here, 'Wind Turbine Blade Optimisation Demo'));
W = WindModel(P.Vmean);                                % 1 m/s bins, 3-20 m/s, 2160 h
end

function [c, twist, od] = buildBlade(x, P, L, B)
% Peak-Cp construction at the design wind: AoA = best L/D + linear offset
o = struct('lookup', L, 'tipLoss', true, 'alphaMode', 'bestLD');
od = BEMT(B, [P.Rhub P.R], x(1)/x(2), x(3:4), x(5:6), P.N, x(2), P.nu, P.tol, P.maxIt, o);
c = od.c;  twist = od.twist;
end

function [x, E, ok, nf, flag] = oneStart(z0, lb, ub, fo, P, L, W, B, mode)
toX = @(z) lb + z(:)'.*(ub - lb);
zC = []; rC = [];
    function r = evalZ(z)
        if isempty(zC) || any(z(:)' ~= zC)
            x_ = toX(z);
            [c, tw, od] = buildBlade(designPoint(x_, mode), P, L, B);
            m = real(od.F) > 1e-9 & isfinite(od.Re);
            r.od = od; r.E = -inf; r.allConv = false;
            if all(od.converged) && any(m)
                [r.E, d] = WinterEnergy(B, c, tw, control(x_, mode), W, P, L);
                r.allConv = d.allConverged;
            end
            r.m = m; zC = z(:)'; rC = r;
        end
        r = rC;
    end
    function f = obj(z)
        r = evalZ(z); f = -r.E/1000; if ~isfinite(f), f = 1; end
    end
    function [cc, ceq] = con(z)
        r = evalZ(z); od = r.od; m = r.m; ceq = [];
        if ~any(m), cc = ones(5,1); return; end
        cc = [(L.ReRange(1) - min(od.Re(m)))/L.ReRange(1);
              (max(od.Re(m)) - L.ReRange(2))/L.ReRange(2);
              (L.alphaRange(1) - min(od.alpha(m)))/12;
              (max(od.alpha(m)) - L.alphaRange(2))/12;
              0.5 - double(all(od.converged) && r.allConv)];
        cc(~isfinite(cc)) = 1;
    end
A = [0 0 -1 1 0 0]; A = A.*(ub - lb);          % c_tip <= c_root, in z units
[z, ~, flag, outp] = fmincon(@obj, z0, A, 0, [], [], zeros(1,6), ones(1,6), @con, fo);
rr = evalZ(z); [cc2, ~] = con(z);
x = toX(z); E = rr.E; ok = flag > 0 && all(cc2 <= 1e-6); nf = outp.funcCount;
end

% ------------------------------------------------------------------------
function combine(here, P, L, W, B, peak)
runs = [];
for f = dir(fullfile(here, 'winter_part_*.mat'))'
    S = load(fullfile(here, f.name)); runs = [runs, S.runs]; %#ok<AGROW>
end
[~, i] = sort([runs.start]); runs = runs(i);
Ef = [runs.E]; Ef(~[runs.feasible]) = -inf;
[Ebest, j] = max(Ef); xb = runs(j).x;
nNear = sum(Ef >= Ebest*(1 - 1e-3));
fprintf('\n%d starts, %d feasible, %d within 0.1%% of best\n', numel(runs), sum([runs.feasible]), nNear);
fprintf(' start   E(kWh)  feasible  tip    Vd    c_root c_tip  dA_root dA_tip  nFun  time(s)\n');
for r = runs
    fprintf('%4d   %7.1f   %d      %5.1f %5.2f  %5.3f  %5.3f  %6.2f  %6.2f  %5d %6.0f\n', ...
        r.start, r.E, r.feasible, r.x, r.nFun, r.time);
end

% the three designs
[cE, twE, odE] = buildBlade(xb, P, L, B);
[EE, dE] = WinterEnergy(B, cE, twE, struct('mode', 'fixed', 'tip', xb(1)), W, P, L);
pk = peak.s2; cP = pk.out.c; twP = pk.out.twist;
[EV, dV] = WinterEnergy(B, cP, twP, struct('mode', 'variable', 'lamOpt', pk.x(1)), W, P, L);
fE = @(tip) -WinterEnergy(B, cP, twP, struct('mode', 'fixed', 'tip', tip), W, P, L);
tipP = fminbnd(fE, 20, P.tipMax, optimset('TolX', 0.01));
[EP, dP] = WinterEnergy(B, cP, twP, struct('mode', 'fixed', 'tip', tipP), W, P, L);
fprintf('\nWinter energy, B = 3, Rayleigh mean %.2f m/s, %d h:\n', W.Vmean, W.hours);
fprintf('  variable speed, peak-Cp blade          %6.1f kWh  (100%%)\n', EV);
fprintf('  fixed speed,    peak-Cp blade (tip %4.1f) %6.1f kWh  (%.1f%%)\n', tipP, EP, 100*EP/EV);
fprintf('  fixed speed,    energy-optimised (tip %4.1f, Vd %.2f) %6.1f kWh  (%.1f%%)\n', xb(1), xb(2), EE, 100*EE/EV);
fprintf('  energy-optimised gain over peak-Cp blade at fixed speed: %+.1f%%\n', 100*(EE/EP - 1));
fprintf('  energy-optimised blade: peak Cp %.4f vs %.4f, design lambda %.2f\n', odE.Cp, pk.Cp, xb(1)/xb(2));

% checks and sensitivity
fx = @(tip) struct('mode', 'fixed', 'tip', tip);
Wf = WindModel(W.Vmean, 0.25);
chk.binsE  = [WinterEnergy(B, cE, twE, fx(xb(1)), Wf, P, L), WinterEnergy(B, cP, twP, fx(tipP), Wf, P, L), ...
              WinterEnergy(B, cP, twP, struct('mode', 'variable', 'lamOpt', pk.x(1)), Wf, P, L)];
chk.clampE = [WinterEnergy(B, cE, twE, fx(xb(1)), W, P, L, 'clamp'), WinterEnergy(B, cP, twP, fx(tipP), W, P, L, 'clamp')];
fprintf('\nCheck 0.25 m/s bins: %.1f / %.1f / %.1f kWh (1 m/s bins: %.1f / %.1f / %.1f)\n', chk.binsE, EE, EP, EV);
fprintf('Stall model: clamp instead of Viterna -> %.1f / %.1f kWh (Viterna %.1f / %.1f)\n', chk.clampE, EE, EP);
Vm = [5.880764 6.153036 6.294867]; sens = zeros(3,3);
for k = 1:3
    Wk = WindModel(Vm(k));
    sens(k,:) = [WinterEnergy(B, cP, twP, struct('mode', 'variable', 'lamOpt', pk.x(1)), Wk, P, L), ...
                 WinterEnergy(B, cP, twP, fx(tipP), Wk, P, L), WinterEnergy(B, cE, twE, fx(xb(1)), Wk, P, L)];
    fprintf('Mean wind %.2f m/s: variable %.1f, peak-Cp fixed %.1f, energy-opt fixed %.1f kWh\n', Vm(k), sens(k,:));
end

% Cp-lambda of both blades (fixed shape)
lamS = linspace(1, 14, 105); cpS = nan(2, numel(lamS));
o = struct('lookup', L, 'tipLoss', true, 'alphaMode', 'twist', 'relax', 0.3, 'stall', 'viterna');
for k = 1:numel(lamS)
    s = BEMT(B, [P.Rhub P.R], lamS(k), cP, twP, P.N, P.V, P.nu, P.tol, 5000, o); if all(s.converged), cpS(1,k) = s.Cp; end
    s = BEMT(B, [P.Rhub P.R], lamS(k), cE, twE, P.N, P.V, P.nu, P.tol, 5000, o); if all(s.converged), cpS(2,k) = s.Cp; end
end
win = struct('runs', runs, 'xb', xb, 'EE', EE, 'EP', EP, 'EV', EV, 'tipP', tipP, 'dE', dE, 'dP', dP, 'dV', dV, ...
    'cE', cE, 'twE', twE, 'odE', odE, 'chk', chk, 'sens', sens, 'Vm', Vm, 'lamS', lamS, 'cpS', cpS, 'W', W);
save(fullfile(here, 'results_winter.mat'), 'win', 'P');
figures(win, P, pk, fullfile(here, 'figures'));
end

function combineVar(here, P, L, W, B, peak)
% Variable-speed energy optimisation: best start, checks, added to results_winter.mat
runs = [];
for f = dir(fullfile(here, 'winterVar_part_*.mat'))'
    S = load(fullfile(here, f.name)); runs = [runs, S.runs]; %#ok<AGROW>
end
[~, i] = sort([runs.start]); runs = runs(i);
Ef = [runs.E]; Ef(~[runs.feasible]) = -inf;
[Ebest, j] = max(Ef); xb = runs(j).x;
fprintf('\nVARIABLE speed: %d starts, %d feasible, %d within 0.1%% of best\n', numel(runs), ...
    sum([runs.feasible]), sum(Ef >= Ebest*(1 - 1e-3)));
fprintf(' start   E(kWh)  feasible  lamOpt  Vd    c_root c_tip  dA_root dA_tip  nFun  time(s)\n');
for r = runs
    fprintf('%4d   %7.1f   %d      %5.2f %5.2f  %5.3f  %5.3f  %6.2f  %6.2f  %5d %6.0f\n', ...
        r.start, r.E, r.feasible, r.x, r.nFun, r.time);
end
S = load(fullfile(here, 'results_winter.mat')); w = S.win;
[cVo, twVo, odVo] = buildBlade(designPoint(xb, 'variable'), P, L, B);
[EVo, dVo] = WinterEnergy(B, cVo, twVo, control(xb, 'variable'), W, P, L);
pk = peak.s2;
fprintf('\nWinter energy, B = 3, Rayleigh mean %.2f m/s, %d h:\n', W.Vmean, W.hours);
fprintf('  variable speed, energy-optimised (lamOpt %.2f, Vd %.2f) %6.1f kWh  (100%%)\n', xb(1), xb(2), EVo);
fprintf('  variable speed, peak-Cp blade          %6.1f kWh  (%.1f%%)\n', w.EV, 100*w.EV/EVo);
fprintf('  fixed speed,    energy-optimised       %6.1f kWh  (%.1f%%)\n', w.EE, 100*w.EE/EVo);
fprintf('  fixed speed,    peak-Cp blade          %6.1f kWh  (%.1f%%)\n', w.EP, 100*w.EP/EVo);
fprintf('  variable: energy-optimised gain over peak-Cp blade: %+.1f%%\n', 100*(EVo/w.EV - 1));
fprintf('  variable vs fixed (both energy-optimised): %+.1f%%\n', 100*(EVo/w.EE - 1));
fprintf('  energy above the tip limit (V > %.1f m/s): %.1f of %.1f kWh\n', P.tipMax/xb(1), ...
    sum(dVo.Ebin(dVo.V > P.tipMax/xb(1))), EVo);
fprintf('  blade: peak Cp at design %.4f vs %.4f; chord %.0f -> %.0f mm\n', odVo.Cp, pk.Cp, 1000*cVo(1), 1000*cVo(end));
ctl = control(xb, 'variable');
chk.binsE  = WinterEnergy(B, cVo, twVo, ctl, WindModel(W.Vmean, 0.25), P, L);
chk.clampE = WinterEnergy(B, cVo, twVo, ctl, W, P, L, 'clamp');
fprintf('Check 0.25 m/s bins: %.1f kWh; clamp instead of Viterna: %.1f kWh (peak-Cp blade, variable, clamp: %.1f)\n', ...
    chk.binsE, chk.clampE, WinterEnergy(B, pk.out.c, pk.out.twist, struct('mode', 'variable', 'lamOpt', pk.x(1)), W, P, L, 'clamp'));
sens = zeros(1, 3);
for k = 1:3, sens(k) = WinterEnergy(B, cVo, twVo, ctl, WindModel(w.Vm(k)), P, L); end
fprintf('Mean wind %.2f / %.2f / %.2f m/s: %.1f / %.1f / %.1f kWh (peak-Cp variable %.1f / %.1f / %.1f)\n', w.Vm, sens, w.sens(:,1));
o = struct('lookup', L, 'tipLoss', true, 'alphaMode', 'twist', 'relax', 0.3, 'stall', 'viterna');
cpV = nan(1, numel(w.lamS));
for k = 1:numel(w.lamS)
    s = BEMT(B, [P.Rhub P.R], w.lamS(k), cVo, twVo, P.N, P.V, P.nu, P.tol, 5000, o); if all(s.converged), cpV(k) = s.Cp; end
end
winVar = struct('runs', runs, 'xb', xb, 'EVo', EVo, 'dVo', dVo, 'cVo', cVo, 'twVo', twVo, 'odVo', odVo, ...
    'chk', chk, 'sens', sens, 'cpV', cpV);
win = w; save(fullfile(here, 'results_winter.mat'), 'win', 'winVar', 'P');
figures(w, P, pk, fullfile(here, 'figures'), winVar);
end

function figures(w, P, pk, figDir, wv)
if nargin < 5, wv = []; end
cols = lines(4);
names = {'Variable speed, peak-C_p blade', sprintf('Fixed speed, peak-C_p blade (tip %.0f m/s)', w.tipP), ...
         sprintf('Fixed speed, energy-optimised blade (tip %.0f m/s)', w.xb(1))};
ds = {w.dV, w.dP, w.dE};
if ~isempty(wv)
    names{4} = sprintf('Variable speed, energy-optimised blade (\\lambda_{opt} %.2f)', wv.xb(1));
    ds{4} = wv.dVo;
end
nd = numel(ds);
% power curve + energy per bin
f = figure('Visible', 'off', 'Position', [0 0 1000 400]);
tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on;
for k = 1:nd, plot(ds{k}.V, ds{k}.P, 'o-', 'Color', cols(k,:), 'LineWidth', 1.5, 'MarkerSize', 4); end
xlabel('Wind speed (m/s)'); ylabel('Power (W)'); legend(names, 'Location', 'northwest'); title('Power curves');
nexttile; hold on; grid on;
Eb = cell2mat(cellfun(@(d) d.Ebin(:), ds, 'UniformOutput', false));
hb = bar(w.dV.V, Eb, 1); for k = 1:nd, hb(k).FaceColor = cols(k,:); end
yyaxis right; plot(w.W.Vmid, w.W.prob*100, 'k--', 'LineWidth', 1.2); ylabel('Rayleigh probability (%)');
yyaxis left; ylabel('Winter energy per 1 m/s bin (kWh)'); xlabel('Wind speed (m/s)');
Etot = cellfun(@(d) sum(d.Ebin), ds);
title(['Energy per bin: ' strjoin(arrayfun(@(e) sprintf('%.0f', e), Etot, 'UniformOutput', false), ' / ') ' kWh']);
saveBoth(f, 't3_winter_power_energy', figDir);
% Cp-lambda and blade shapes
f = figure('Visible', 'off', 'Position', [0 0 1000 400]);
tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on;
plot(w.lamS, w.cpS(1,:), 'Color', cols(2,:), 'LineWidth', 1.5); plot(w.lamS, w.cpS(2,:), 'Color', cols(3,:), 'LineWidth', 1.5);
lg = {'peak-C_p blade', 'energy-opt. (fixed speed)'};
if ~isempty(wv), plot(w.lamS, wv.cpV, 'Color', cols(4,:), 'LineWidth', 1.5); lg{3} = 'energy-opt. (variable speed)'; end
xlabel('\lambda'); ylabel('C_p'); ylim([0 0.5]); legend(lg, 'Location', 'south');
title('C_p-\lambda (fixed blade shape)');
nexttile; hold on; grid on;
plot(pk.out.r, 1000*pk.out.c, 'Color', cols(2,:), 'LineWidth', 1.5); plot(w.odE.r, 1000*w.cE, 'Color', cols(3,:), 'LineWidth', 1.5);
if ~isempty(wv), plot(wv.odVo.r, 1000*wv.cVo, 'Color', cols(4,:), 'LineWidth', 1.5); end
xlabel('r (m)'); ylabel('Chord (mm)'); title('Chord');
nexttile; hold on; grid on;
plot(pk.out.r, pk.out.twist, 'Color', cols(2,:), 'LineWidth', 1.5); plot(w.odE.r, w.twE, 'Color', cols(3,:), 'LineWidth', 1.5);
if ~isempty(wv), plot(wv.odVo.r, wv.twVo, 'Color', cols(4,:), 'LineWidth', 1.5); end
xlabel('r (m)'); ylabel('Twist (deg)'); title('Twist');
saveBoth(f, 't3_winter_blades', figDir);
end

function saveBoth(f, n, figDir)
exportgraphics(f, fullfile(figDir, [n '.pdf']));
exportgraphics(f, fullfile(figDir, [n '.png']), 'Resolution', 200);
close(f);
end
