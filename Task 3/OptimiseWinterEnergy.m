function OptimiseWinterEnergy(action, worker, nWorkers)
% OPTIMISEWINTERENERGY  Option 4: blade + fixed rotor speed for maximum winter energy (B = 3).
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

if nargin < 1, action = 'quick'; end
here = fileparts(mfilename('fullpath'));
[P, L, W, B, peak] = setup(here);

switch action
    case {'worker', 'quick'}
        if strcmp(action, 'quick'), worker = 1; nWorkers = 1; P.nStart = 1; P.fmIt = 3; end
        lb = [20, 4, P.cMin, P.cMin, -P.dAmax, -P.dAmax];
        ub = [P.tipMax, 14, P.cMax, P.cMax, P.dAmax, P.dAmax];
        rng(3); Z0 = rand(P.nStart, 6); Z0(:,3:4) = sort(Z0(:,3:4), 2, 'descend');
        x1 = peak.s1.x;                                   % start 1: peak-Cp linear blade at max tip speed
        Z0(1,:) = ([P.tipMax, P.tipMax/x1(1), x1(2:5)] - lb)./(ub - lb);
        Z0 = min(max(Z0, 0), 1);
        fo = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
            'MaxFunctionEvaluations', 1500, 'MaxIterations', P.fmIt, ...
            'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8, ...
            'ConstraintTolerance', 1e-6, 'FiniteDifferenceStepSize', 1e-5);
        mine = worker:nWorkers:P.nStart;
        runs = struct('start', {}, 'x', {}, 'E', {}, 'feasible', {}, 'nFun', {}, 'time', {}, 'flag', {});
        for s = mine
            t0 = tic;
            [x, E, ok, nf, flag] = oneStart(Z0(s,:), lb, ub, fo, P, L, W, B);
            runs(end+1) = struct('start', s, 'x', x, 'E', E, 'feasible', ok, 'nFun', nf, 'time', toc(t0), 'flag', flag); %#ok<AGROW>
            fprintf('worker %d start %d: E = %.1f kWh, feasible %d, x = %s, %.0f s\n', ...
                worker, s, E, ok, mat2str(x, 4), toc(t0));
            if ~strcmp(action, 'quick'), save(fullfile(here, sprintf('winter_part_%d.mat', worker)), 'runs'); end
        end
    case 'combine'
        combine(here, P, L, W, B, peak);
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

function [x, E, ok, nf, flag] = oneStart(z0, lb, ub, fo, P, L, W, B)
toX = @(z) lb + z(:)'.*(ub - lb);
zC = []; rC = [];
    function r = evalZ(z)
        if isempty(zC) || any(z(:)' ~= zC)
            x_ = toX(z);
            [c, tw, od] = buildBlade(x_, P, L, B);
            m = real(od.F) > 1e-9 & isfinite(od.Re);
            r.od = od; r.E = -inf; r.allConv = false;
            if all(od.converged) && any(m)
                [r.E, d] = WinterEnergy(B, c, tw, struct('mode', 'fixed', 'tip', x_(1)), W, P, L);
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

function figures(w, P, pk, figDir)
cols = lines(3);
names = {'Variable speed, peak-C_p blade', sprintf('Fixed speed, peak-C_p blade (tip %.0f m/s)', w.tipP), ...
         sprintf('Fixed speed, energy-optimised blade (tip %.0f m/s)', w.xb(1))};
ds = {w.dV, w.dP, w.dE};
% power curve + energy per bin
f = figure('Visible', 'off', 'Position', [0 0 1000 400]);
tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on;
for k = 1:3, plot(ds{k}.V, ds{k}.P, 'o-', 'Color', cols(k,:), 'LineWidth', 1.5, 'MarkerSize', 4); end
xlabel('Wind speed (m/s)'); ylabel('Power (W)'); legend(names, 'Location', 'northwest'); title('Power curves');
nexttile; hold on; grid on;
bar(w.dV.V, [w.dV.Ebin; w.dP.Ebin; w.dE.Ebin]', 1);
yyaxis right; plot(w.W.Vmid, w.W.prob*100, 'k--', 'LineWidth', 1.2); ylabel('Rayleigh probability (%)');
yyaxis left; ylabel('Winter energy per 1 m/s bin (kWh)'); xlabel('Wind speed (m/s)');
title(sprintf('Energy per bin: %.0f / %.0f / %.0f kWh', w.EV, w.EP, w.EE));
saveBoth(f, 't3_winter_power_energy', figDir);
% Cp-lambda and blade shapes
f = figure('Visible', 'off', 'Position', [0 0 1000 400]);
tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on;
plot(w.lamS, w.cpS(1,:), 'Color', cols(2,:), 'LineWidth', 1.5); plot(w.lamS, w.cpS(2,:), 'Color', cols(3,:), 'LineWidth', 1.5);
xlabel('\lambda'); ylabel('C_p'); ylim([0 0.5]); legend('peak-C_p blade', 'energy-optimised blade', 'Location', 'south');
title('C_p-\lambda (fixed blade shape)');
nexttile; hold on; grid on;
plot(pk.out.r, 1000*pk.out.c, 'Color', cols(2,:), 'LineWidth', 1.5); plot(w.odE.r, 1000*w.cE, 'Color', cols(3,:), 'LineWidth', 1.5);
xlabel('r (m)'); ylabel('Chord (mm)'); title('Chord');
nexttile; hold on; grid on;
plot(pk.out.r, pk.out.twist, 'Color', cols(2,:), 'LineWidth', 1.5); plot(w.odE.r, w.twE, 'Color', cols(3,:), 'LineWidth', 1.5);
xlabel('r (m)'); ylabel('Twist (deg)'); title('Twist');
saveBoth(f, 't3_winter_blades', figDir);
end

function saveBoth(f, n, figDir)
exportgraphics(f, fullfile(figDir, [n '.pdf']));
exportgraphics(f, fullfile(figDir, [n '.png']), 'Resolution', 200);
close(f);
end
