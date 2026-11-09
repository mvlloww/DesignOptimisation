function out = bb_session(cmd, varargin)
% BB_SESSION  Whole black-box workflow in one file.
%
%   bb_session('plan')      before the session: build + check sampling plan
%   bb_session              run / resume the session (interactive)
%   bb_session('status')    model table, decision checks, figures
%   bb_session('undo')      delete the last logged row (typo fix)
%   bb_session('figure')    report figure of the final model (fig_final_model.png)
%   bb_session('rehearse')  full dry run on a fake function (deletes log!)
%   bb_session('selftest')  check toolboxes + closed-form LOOCV vs brute force
%   bb_session('maxcheck')  sensitivity of the surrogate optimum to grid size
%                           and number of seeds in maximise (uses current log)
%   bb_session('plancompare') chosen LHS vs Halton / random / other LHS
%                           (fig_plan_compare.png; uses bb_plan.csv only)
%   bb_session('study')     nTrain sweep: full rehearsals on the easy AND hard
%                           fake boxes for several nTrain and seeds
%                           (fig_ntrain_study.png). Uses its own study_*
%                           files, so the real log is safe.
%                           bb_session('study', nTrainList, seedList, kinds)
%   bb_session('rehearse','hard')  rehearse on the hard (two-peak) fake box
%
% At the MCSI> prompt type ONE of:
%   7.42                 reading, knobs set exactly as suggested
%   101.5 35 0.81 7.42   x1 x2 x3 MCSI as ACTUALLY set (continuous knobs)
%   u = undo last   s = status   q = quit (resume later, nothing lost)
%   x / e / l / c = override this call: exploit / explore / local / confirm
%
% Everything is logged to bb_log.csv after every call. Restarting MATLAB
% and calling bb_session again resumes where you left off.
%
% Phase codes in the log: 1 LHS, 2 hold-out, 3 repeat (optional), 4 exploit,
%                         5 explore, 6 local, 7 confirm
% The box is assumed deterministic. The log always holds the settings
% actually dialled, so knob-setting error never enters the model.
%
% Toolboxes: Statistics and Machine Learning (lhsdesign, pdist, pdist2,
% fitlm, fitrgp, quantile, haltonset), Optimization (fmincon) and Global
% Optimization (MultiStart; optional, falls back to an fmincon loop).

if nargin < 1, cmd = 'run'; end
cfg = config();
switch lower(cmd)
    case 'plan',     out = make_plan(cfg, true);
    case 'run',      run_loop(cfg, [], true);
    case 'status',   L = load_log(cfg); status(L, cfg, true);
    case 'undo',     undo_last(cfg);
    case 'figure'
        L = load_log(cfg); F = fit_best(L, cfg); uS = maximise(F.m);
        final_figure(L, F, uS, cfg);
    case 'rehearse', out = rehearse(cfg, varargin{:});
    case 'selftest', out = selftest(cfg);
    case 'maxcheck', out = maxcheck(cfg);
    case 'plancompare', out = plancompare(cfg);
    case 'study',    out = study(cfg, varargin{:});
    otherwise, error('unknown command %s', cmd);
end
end

% =====================================================================
%  CONFIG
% =====================================================================
function cfg = config()
cfg.lb    = [0   0  0];
cfg.ub    = [200 60 1];
cfg.res   = [0.1   0.1 0.01];        % how precisely you can SET each knob (CHECK ON THE DAY)
                                     % x1 0.0-200.0 (£k), x2 0.0-60.0 (£k), x3 0.00-1.00
cfg.yRange    = [0 10];             % MCSI display range (0-10.0000)
cfg.yDecimals = 4;                  % MCSI display precision
cfg.names = {'x1 advert (£k)','x2 targeting (£k)','x3 skew (-)'};

cfg.budget  = 50;  cfg.nTrain = 25;  cfg.nTest = 5;  cfg.nConfirm = 3;
cfg.nRepeat = 0;     % 1 = repeat LHS point 1 after training, to measure how
                     %     reproducibly the knobs can be set (no readout on box)
cfg.nPlans  = 200;   % maximin LHS plans from lhsdesign ...
cfg.lhsIter = 50;    % ... each the best of 50 (10,000 LHS screened)
cfg.seed    = 27;

cfg.polyForms  = {'linear','quadratic','poly333'};   % fitlm model specs
cfg.kernels    = {'gaussian','multiquadric','invmultiquadric','cubic','thinplate'};
cfg.widths     = logspace(-1.3, 0.5, 20);
cfg.lambda     = 1e-10;                  % RBF nugget (deterministic box)
cfg.useGP      = true;                   % add fitrgp (Kriging) to the comparison

% decision rules
cfg.stableTol   = 0.02;  % optimum moved < 2% of each range
cfg.errTol      = 0.01;  % last exploit error < 1% of observed MCSI range
cfg.peakTol     = 0.05;  % rival peak within 5% of range = competitor
cfg.q2High      = 0.80;  % above this, exploit 2:1, else 1:1
cfg.nLocalMin   = 6;     % need >= 6 calls (plus confirm) to start local
cfg.nearOptFrac = 0.02;  % local box = region within 2% of range of max
cfg.minBoxSteps = 5;     % local box half-width >= 5 dial steps
cfg.confirmStep = 0.05;  % confirm steps of 5% of range
cfg.sliceWarn   = 0.4;   % warn if an LHS point is set > 0.4 slice from target

cfg.rehearseSetErr = [2 0.75 0.01];  % rehearsal only: knob setting error (+/- uniform)

cfg.logFile   = 'bb_log.csv';
cfg.planFile  = 'bb_plan.csv';
cfg.stateFile = 'bb_state.mat';   % cache of past surrogate optima (speed only)
end

% =====================================================================
%  SESSION LOOP
% =====================================================================
function run_loop(cfg, box, verbose)
% box = [] -> readings typed by hand; box = @f -> rehearsal (no prompts)
interactive = isempty(box);
if ~exist(cfg.planFile,'file'), make_plan(cfg, verbose); end
plan = dlmread(cfg.planFile, ',', 1, 0);
override = '';
while true
    L = load_log(cfg);
    n = size(L,1);
    if n >= cfg.budget
        if verbose
            fprintf('\nBudget of %d calls used.\n', cfg.budget);
            status(L, cfg, interactive);
        end
        return;
    end
    if n == size(plan,1) && isempty(override) && verbose
        fprintf('\n########## FIXED BLOCK DONE: DECISION POINT ##########');
        status(L, cfg, false);          % train-vs-hold-out table before infill
    end
    P = next_point(L, plan, cfg, override, verbose);
    override = '';

    if ~isempty(box)
        xa = P.x + cfg.rehearseSetErr .* (2*rand(1,3) - 1);   % imperfect knob setting
        xa = min(max(xa, cfg.lb), cfg.ub);
        y  = box(xa);
    else
        s = strtrim(input('MCSI> ', 's'));
        if isempty(s), continue; end
        switch lower(s)
            case 'q', fprintf('Stopped. Run bb_session to resume.\n'); return;
            case 'u', undo_last(cfg); continue;
            case 's', status(L, cfg, true); continue;
            case {'x','e','l','c'}
                m = struct('x','exploit','e','explore','l','local','c','confirm');
                override = m.(lower(s)); continue;
        end
        v = sscanf(strrep(s, ',', ' '), '%f')';
        if numel(v) == 1
            xa = P.x; y = v;
        elseif numel(v) == 4
            xa = v(1:3); y = v(4);
            if any(xa < cfg.lb | xa > cfg.ub)
                fprintf('  settings outside bounds: [%g %g %g]\n', xa); continue;
            end
        else
            fprintf('  type MCSI, or x1 x2 x3 MCSI (got "%s")\n', s); continue;
        end
        if y < cfg.yRange(1) || y > cfg.yRange(2)
            fprintf('  MCSI must be between %g and %g (got %g): re-enter\n', cfg.yRange, y); continue;
        end
    end
    record(L, P, xa, y, cfg, verbose);
end
end

% =====================================================================
%  CHOOSE NEXT POINT
% =====================================================================
function P = next_point(L, plan, cfg, override, verbose)
n = size(L,1);
P.ypred = NaN; P.sig = NaN; P.note = '';

% ---- fixed block -------------------------------------------------------
if n < size(plan,1) && isempty(override)
    P.x = plan(n+1,3:5); P.phase = plan(n+1,2);
    tr = L(L(:,2)==1,:);
    P.m = [];
    if P.phase ~= 1 && size(tr,1) >= 10
        F = fit_best(tr, cfg);
        P.m = F.m; P.ypred = spredict(F.m, nrm(P.x,cfg)); P.sig = F.rmse;
    end
    show(P, n, cfg, verbose);
    return;
end

% ---- adaptive ---------------------------------------------------------
F   = fit_best(L, cfg);
[uS, yS, peaks] = maximise(F.m);
cache_opt(L, uS, cfg);
D   = decision_checks(L, uS, yS, peaks, cfg);
nLeft = cfg.budget - n;
Udata = nrm(L(:,3:5), cfg);

if ~isempty(override), mode = override;
elseif nLeft <= cfg.nConfirm, mode = 'confirm';
elseif any(L(:,2)==6), mode = 'local';
elseif D.all, mode = 'local';
else
    nG = sum(L(:,2)==4 | L(:,2)==5);
    if F.q2 >= cfg.q2High && ~D.rival, period = 3; else, period = 2; end
    if mod(nG, period) == period - 1, mode = 'explore'; else, mode = 'exploit'; end
end

switch mode
    case 'exploit'
        u = uS; P.phase = 4;
        if min(pdist2(u, Udata)) < cfg.stableTol || is_sampled(dn(u,cfg), L, cfg)
            P.note = 'optimum already sampled -> explore instead';
            mode = 'explore';
        end
    case 'local'
        P.phase = 6;
        [lo, hi] = local_box(F.m, uS, yS, L, cfg);
        u = spacefill(lo, hi, Udata, L, cfg);
        P.note = sprintf('local box: x1 %g-%g | x2 %g-%g | x3 %g-%g', ...
            reshape([dn(lo,cfg); dn(hi,cfg)], 1, []));
    case 'confirm'
        P.phase = 7;
        if min(pdist2(uS, Udata)) >= cfg.stableTol && ~is_sampled(dn(uS,cfg), L, cfg)
            u = uS; P.note = 'confirm: surrogate optimum';
        else
            ub = nrm(L(best_obs(L),3:5), cfg);
            k = mod(sum(L(:,2)==7), 3) + 1;
            for sgn = [1 -1]
                u = ub; u(k) = min(max(u(k) + sgn*cfg.confirmStep, 0), 1);
                if ~is_sampled(dn(u,cfg), L, cfg), break; end
            end
            P.note = sprintf('confirm: 5%% step in %s from best observed', cfg.names{k});
        end
end
if strcmp(mode, 'explore')
    % most isolated point among the top 25% of predicted MCSI
    P.phase = 5;
    C = rand(20000,3); yc = spredict(F.m, C);
    C = C(yc >= quantile(yc, 0.75), :);
    [~, k] = max(min(pdist2(C, Udata), [], 2)); u = C(k,:);
end
P.x = dn(u, cfg);
P.m = F.m; P.ypred = spredict(F.m, nrm(P.x, cfg)); P.sig = F.rmse;
if verbose
    fprintf('\nModel %s (LOO Q2 %.3f) | checks: stable %d, accurate %d, no rival %d, calls %d -> %s\n', ...
        F.name, F.q2, D.stable, D.accurate, ~D.rival && D.explored, D.calls, upper(mode));
end
show(P, n, cfg, verbose);
end

function D = decision_checks(L, uS, yS, peaks, cfg)
n = size(L,1); yr = range_(L(:,6));
% 1 stable: optimum with vs without the latest point
uP = cache_opt(L(1:end-1,:), [], cfg);
if isempty(uP)
    Fp = fit_best(L(1:end-1,:), cfg); uP = maximise(Fp.m); cache_opt(L(1:end-1,:), uP, cfg);
end
D.stable = norm(uS - uP) < cfg.stableTol;
% 2 accurate: last exploit/local prediction error
r = find(L(:,2)==4 | L(:,2)==6, 1, 'last');
D.accurate = ~isempty(r) && ~isnan(L(r,7)) && abs(L(r,6)-L(r,7)) < cfg.errTol*yr;
% 3 no rival peak, explored enough, explore never beat the optimum
D.rival = size(peaks,1) > 1 && (yS - peaks(2,4)) < cfg.peakTol*yr;
ex = L(L(:,2)==5, 6);
D.explored = numel(ex) >= 2 && all(ex <= yS);
% 4 enough calls left
D.calls = (cfg.budget - n - cfg.nConfirm) >= cfg.nLocalMin;
D.all = D.stable && D.accurate && ~D.rival && D.explored && D.calls;
end

function u = cache_opt(L, u, cfg)
% Stores the surrogate optimum for exactly this data set, so the
% "stable?" check does not refit every model a second time each call.
% Keyed on the data itself, so undo/re-entry can never return stale values.
S.keys = {}; S.vals = {};
if exist(cfg.stateFile, 'file'), S = load(cfg.stateFile); end
key = L(:,3:6);
if isempty(u)                                   % lookup
    for i = 1:numel(S.keys)
        if isequal(S.keys{i}, key), u = S.vals{i}; return; end
    end
else                                            % store
    S.keys{end+1} = key; S.vals{end+1} = u;
    keys = S.keys; vals = S.vals;               %#ok<NASGU>
    save(cfg.stateFile, 'keys', 'vals');
end
end

function [lo, hi] = local_box(m, uS, yS, L, cfg)
% near-optimal region: predicted MCSI within nearOptFrac*range of the max
thr  = yS - cfg.nearOptFrac * range_(L(:,6));
minW = cfg.minBoxSteps * cfg.res ./ (cfg.ub - cfg.lb);
lo = uS; hi = uS; t = (0.005:0.005:1)';
for k = 1:3
    for sgn = [-1 1]
        Q = repmat(uS, numel(t), 1); Q(:,k) = uS(k) + sgn*t;
        ok = Q(:,k) >= 0 & Q(:,k) <= 1;
        ok(ok) = spredict(m, Q(ok,:)) >= thr;
        last = find(~ok, 1) - 1; if isempty(last), last = numel(t); end
        if last > 0
            if sgn < 0, lo(k) = Q(last,k); else, hi(k) = Q(last,k); end
        end
    end
end
lo = max(min(lo, uS - minW), 0);
hi = min(max(hi, uS + minW), 1);
end

function u = spacefill(lo, hi, Udata, L, cfg)
% most isolated point inside the box, distances in box-scaled units
w = hi - lo;
C = lo + lhsdesign(2000, 3) .* w;
d = min(pdist2(C ./ w, Udata ./ w), [], 2);
[~, o] = sort(d, 'descend');
for k = o'
    if ~is_sampled(dn(C(k,:),cfg), L, cfg), u = C(k,:); return; end
end
u = C(o(1),:);
end

% =====================================================================
%  LOGGING
% =====================================================================
function record(L, P, xa, y, cfg, verbose)
% Logs the ACTUAL settings xa. Prediction error is re-evaluated at xa, so
% it measures model error only, not knob-setting error.
ypred = P.ypred;
if ~isempty(P.m), ypred = spredict(P.m, nrm(xa, cfg)); end
row = [size(L,1)+1, P.phase, xa, y, ypred];
write_log([L; row], cfg);
if ~verbose, return; end
yr = range_([L(:,6); y]);
fprintf('  logged run %d: [%g %g %g] MCSI = %.6g', row(1), xa, y);
if ~isnan(ypred)
    fprintf('  (pred %.5g, error %+.1f%% of range)', ypred, 100*(y-ypred)/yr);
end
fprintf('\n');
if any(abs(xa - P.x) > 1e-9)
    fprintf('  set off target by [%+.3g %+.3g %+.3g]\n', xa - P.x);
    slice = (cfg.ub - cfg.lb) / cfg.nTrain;
    if P.phase == 1 && any(abs(xa - P.x) > cfg.sliceWarn * slice)
        fprintf('  !! more than %.0f%% of an LHS slice off target: u and re-set if you can\n', 100*cfg.sliceWarn);
    end
end
if ~isnan(ypred) && ~isnan(P.sig) && abs(y - ypred) > 3*P.sig
    fprintf('  !! reading is %.1f LOO-RMSEs from prediction: CHECK FOR A TYPO (u = undo)\n', ...
        abs(y-ypred)/P.sig);
end
if P.phase == 3
    fprintf('  repeat of run 1: settings differ by [%+.3g %+.3g %+.3g], MCSI by %+.4g (%.2f%% of range)\n', ...
        xa - L(1,3:5), y - L(1,6), 100*abs(y - L(1,6))/yr);
    fprintf('     -> this is your knob-setting reproducibility; quote it with the recommendation\n');
end
end

function undo_last(cfg)
L = load_log(cfg);
if isempty(L), fprintf('Log is empty.\n'); return; end
fprintf('Removed run %d: x = [%g %g %g], MCSI = %g\n', L(end,1), L(end,3:5), L(end,6));
write_log(L(1:end-1,:), cfg);
end

function L = load_log(cfg)
L = zeros(0,7);
if exist(cfg.logFile,'file')
    t = dlmread(cfg.logFile, ',', 1, 0);
    if ~isempty(t), L = t; end
end
end

function write_log(L, cfg)
fid = fopen(cfg.logFile, 'w');
fprintf(fid, 'run,phase,x1,x2,x3,y,ypred\n');
fprintf(fid, '%d,%d,%.6g,%.6g,%.6g,%.10g,%.10g\n', L');
fclose(fid);
end

function show(P, n, cfg, verbose)
if ~verbose, return; end
names = {'LHS','hold-out','repeat','exploit','explore','local','confirm'};
fprintf('\n>>> CALL %d/%d [%s]   x1 = %g   x2 = %g   x3 = %g\n', ...
    n+1, cfg.budget, names{P.phase}, P.x);
if ~isnan(P.ypred), fprintf('    predicted MCSI %.5g\n', P.ypred); end
if ~isempty(P.note), fprintf('    %s\n', P.note); end
end

% =====================================================================
%  STATUS AND REPORT FIGURES
% =====================================================================
function status(L, cfg, makePlots)
fprintf('\n===== %d / %d calls =====\n', size(L,1), cfg.budget);
if size(L,1) < 10, return; end
ib = best_obs(L);
fprintf('Best observed: MCSI %.6g at [%g %g %g] (run %d)\n', L(ib,6), L(ib,3:5), L(ib,1));
tr = L(L(:,2)==1,:); te = L(L(:,2)==2,:);
if ~isempty(te)
    fprintf('\n--- train on LHS (%d), test on hold-out (%d) ---\n', size(tr,1), size(te,1));
    compare(tr, cfg, te, true);
end
fprintf('\n--- all data ---\n');
F = compare(L, cfg, [], true);
[uS, yS, peaks] = maximise(F.m);
fprintf('Surrogate optimum [%g %g %g], predicted %.6g\n', dn(uS,cfg), yS);
if size(peaks,1) > 1
    fprintf('Other peaks: ');
    fprintf('[%g %g %g] %.4g;  ', [dn(peaks(2:end,1:3),cfg) peaks(2:end,4)]');
    fprintf('\n');
end
onB = uS < 0.01 | uS > 0.99;
if any(onB), fprintf('!! optimum on bound: %s\n', strjoin(cfg.names(onB), ', ')); end
r = L(:,2) >= 4 & ~isnan(L(:,7));
if any(r)
    fprintf('Predict-then-query errors (%% of range): %s\n', ...
        sprintf('%+.1f ', 100*(L(r,6)-L(r,7))/range_(L(:,6))));
end
recommend(L, F, uS, cfg);
if makePlots, figures(L, F, uS, cfg); end
end

function recommend(L, F, uS, cfg)
% Recommendation = surrogate optimum rounded to setting precision, checked
% against the best measured point, with worst case over +/-1 setting step.
xr = dn(uS, cfg); ib = best_obs(L);
yr = spredict(F.m, nrm(xr, cfg));
[a,b,c] = ndgrid(-1:1); K = [a(:) b(:) c(:)] .* cfg.res;
Xn = min(max(xr + K, cfg.lb), cfg.ub);
yw = min(spredict(F.m, nrm(Xn, cfg)));
fprintf('\n=== RECOMMENDATION ===\n');
fprintf('Set x1 = %g (£k), x2 = %g (£k), x3 = %g   predicted MCSI %.5g\n', xr, yr);
fprintf('Best measured: [%g %g %g] MCSI %.5g (run %d)\n', L(ib,3:5), L(ib,6), L(ib,1));
fprintf('Worst case if each knob is off by one setting step [%g %g %g]: %.5g (%.3f%% lower)\n', ...
        cfg.res, yw, 100*(yr - yw)/abs(yr));
end

function figures(L, F, uS, cfg)
try
    f1 = figure('Color','w'); hold on; box on; grid on;
    plot(F.y, F.yloo, 'ko', 'MarkerFaceColor', [0.2 0.4 0.8]);
    lim = [min([F.y; F.yloo]) max([F.y; F.yloo])]; plot(lim, lim, 'k--');
    xlabel('observed MCSI'); ylabel('LOOCV predicted MCSI');
    title(sprintf('%s: Q^2 = %.3f', F.name, F.q2));
    print(f1, 'fig_parity.png', '-dpng', '-r200');

    final_figure(L, F, uS, cfg);
catch err
    fprintf('(figures skipped: %s)\n', err.message);
end
end

function final_figure(L, F, uS, cfg)
% Report figure of the final surrogate, 2x2:
%  (a)-(c) 2D slices through the recommended point, one shared colour scale,
%          near-optimal region outlined, samples projected onto each slice
%  (d)     MCSI along each input through the optimum (one-at-a-time sensitivity)
blue = [205 226 251; 158 197 244; 109 167 236; 57 135 229; 37 106 191; 24 79 149; 13 54 107]/255;
cmap = interp1(linspace(0,1,7), blue, linspace(0,1,256));        % one hue, light -> dark
ser  = [42 120 214; 235 104 52; 27 175 122]/255;                  % validated categorical
ink  = [0.12 0.12 0.12]; muted = [0.45 0.45 0.45];
xs = dn(uS, cfg); us = nrm(xs, cfg); ys = spredict(F.m, us);
thr = ys - cfg.nearOptFrac * range_(L(:,6));
g = linspace(0, 1, 80); pr = [1 2 3; 1 3 2; 2 3 1];
Z = cell(3,1); AB = cell(3,1);
for p = 1:3
    a = pr(p,1); b = pr(p,2); c = pr(p,3);
    [A,B] = meshgrid(g,g); Q = zeros(numel(A),3); Q(:,a) = A(:); Q(:,b) = B(:); Q(:,c) = us(c);
    Z{p} = reshape(spredict(F.m, Q), size(A)); AB{p} = {A, B};
end
clim_ = [min(cellfun(@(z) min(z(:)), Z)), max(cellfun(@(z) max(z(:)), Z))];
fixed = L(:,2) <= 3; adapt = ~fixed;

f = figure('Color','w', 'Position', [100 100 1100 900]);
lbl = {'(a)','(b)','(c)'};
pos = [0.07 0.58 0.33 0.36; 0.57 0.58 0.33 0.36; 0.07 0.13 0.33 0.36; 0.57 0.13 0.38 0.36];
for p = 1:3
    a = pr(p,1); b = pr(p,2); c = pr(p,3);
    ax = axes('Position', pos(p,:)); hold(ax,'on'); box(ax,'on');
    Xa = cfg.lb(a) + AB{p}{1}*(cfg.ub(a)-cfg.lb(a));
    Xb = cfg.lb(b) + AB{p}{2}*(cfg.ub(b)-cfg.lb(b));
    contourf(ax, Xa, Xb, Z{p}, 16, 'LineColor','none');
    [~, h1] = contour(ax, Xa, Xb, Z{p}, [thr thr], 'LineColor', ink, 'LineWidth', 1.5);
    h2 = plot(ax, L(fixed,2+a), L(fixed,2+b), 'o', 'MarkerSize', 5, ...
         'MarkerFaceColor','w', 'MarkerEdgeColor', ink);
    h3 = plot(ax, L(adapt,2+a), L(adapt,2+b), 'o', 'MarkerSize', 5, ...
         'MarkerFaceColor', ser(2,:), 'MarkerEdgeColor','w');
    h4 = plot(ax, xs(a), xs(b), 'p', 'MarkerSize', 16, 'MarkerFaceColor', ink, 'MarkerEdgeColor','w');
    if p == 1, H = [h1 h2 h3 h4]; ax1 = ax; end
    if p == 3, ax3 = ax; end
    colormap(ax, cmap); caxis(ax, clim_);
    xlabel(ax, cfg.names{a}); ylabel(ax, cfg.names{b});
    title(ax, sprintf('%s slice at %s = %g', lbl{p}, cfg.names{c}, xs(c)), 'FontWeight','normal');
    set(ax, 'Layer','top', 'TickDir','out', 'XColor', muted, 'YColor', muted);
    cb = colorbar(ax); ylabel(cb, 'predicted MCSI');
    set(ax, 'Position', pos(p,:));               % all three slices the same size
    set(cb, 'Position', [pos(p,1)+pos(p,3)+0.012, pos(p,2), 0.012, pos(p,4)]);  % beside, not on top
end
lg = legend(ax1, H, {sprintf('near-optimal region (within %.0f%% of range of max)', 100*cfg.nearOptFrac), ...
       'fixed design (LHS + hold-out)', 'adaptive (infill, local, confirm)', 'recommended'}, ...
       'Orientation','horizontal', 'Box','off');
set(lg, 'Position', [0.05 0.015 0.9 0.03]);

% (d) one-at-a-time sensitivity through the optimum
ax = axes('Position', pos(4,:)); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
t = linspace(0,1,200)'; short = {'advert','targeting','skew'};
for k = 1:3
    Q = repmat(us, numel(t), 1); Q(:,k) = t;
    yk = spredict(F.m, Q);
    plot(ax, 100*t, yk, '-', 'Color', ser(k,:), 'LineWidth', 2);
    text(ax, 99, yk(end), [short{k} ' '], 'Color', ink, 'HorizontalAlignment','right', ...
         'VerticalAlignment','bottom', 'FontSize', 9);      % direct label at line end
end
plot(ax, [0 100], [thr thr], '--', 'Color', muted, 'LineWidth', 1);
plot(ax, 100*us, ys*[1 1 1], 'p', 'MarkerSize', 12, 'MarkerFaceColor', ink, 'MarkerEdgeColor','w');
xlabel(ax, 'input setting (% of range)'); ylabel(ax, 'predicted MCSI');
title(ax, '(d) each input varied alone, others at optimum', 'FontWeight','normal');
legend(ax, {cfg.names{:}, 'near-optimal threshold'}, 'Location','southeast', 'Box','off');
set(ax, 'TickDir','out', 'XColor', muted, 'YColor', muted, 'GridAlpha', 0.15);

print(f, 'fig_final_model.png', '-dpng', '-r300');
fprintf('Saved fig_final_model.png\n');
fprintf('Caption info: %s (LOOCV Q2 = %.3f), %d samples; recommended [%g, %g, %g], predicted MCSI %.4g;\n', ...
        F.name, F.q2, size(L,1), xs, ys);
fprintf('  slices pass through the recommended point; all samples are projected onto each slice.\n');
end

% =====================================================================
%  SAMPLING PLAN
% =====================================================================
function plan = make_plan(cfg, verbose)
rng(cfg.seed, 'twister');   % explicit generator: works even if legacy rand('seed') was used
n = cfg.nTrain; N = cfg.nPlans;
% 1. N maximin LHS plans (lhsdesign keeps the best of lhsIter for each)
C = cell(N,1); md = zeros(N,1); l2 = md;
for i = 1:N
    C{i}  = lhsdesign(n, 3, 'criterion','maximin', 'iterations',cfg.lhsIter, 'smooth','off');
    md(i) = min(pdist(C{i}));
    l2(i) = l2star(C{i});
end
% 2. rank on spread (bigger min distance) and uniformity (lower L2-star)
[~,a] = sort(-md); [~,b] = sort(l2); ra(a) = 1:N; rb(b) = 1:N;
[~, ib] = min(ra + rb);
% 3. snap to dials
Xtr = dn(C{ib}, cfg); Utr = nrm(Xtr, cfg);
if size(unique(Xtr,'rows'),1) < n, error('dial rounding made duplicates'); end
% 4. hold-out: greedy fill of the largest gaps
pool = rand(20000,3); Uall = Utr; Ute = zeros(0,3);
for j = 1:cfg.nTest
    [~,k] = max(min(pdist2(pool, Uall), [], 2));
    Ute(end+1,:) = pool(k,:); Uall(end+1,:) = pool(k,:); pool(k,:) = []; %#ok<AGROW>
end
Xte = dn(Ute, cfg); Ute = nrm(Xte, cfg);
% 5. run sheet: 25 LHS, (optional repeat of #1), 5 hold-out
X   = [Xtr; repmat(Xtr(1,:), cfg.nRepeat, 1); Xte];
phs = [ones(n,1); 3*ones(cfg.nRepeat,1); 2*ones(cfg.nTest,1)];
plan = [(1:size(X,1))' phs X];
fid = fopen(cfg.planFile,'w'); fprintf(fid,'run,phase,x1,x2,x3\n');
fprintf(fid,'%d,%d,%.6g,%.6g,%.6g\n', plan'); fclose(fid);
if exist(cfg.stateFile,'file'), delete(cfg.stateFile); end
if verbose
    fprintf('\n=== Sampling plan: %d LHS + %d hold-out + %d repeat ===\n', n, cfg.nTest, cfg.nRepeat);
    fprintf('min distance %.3f  (median of %d maximin plans %.3f)\n', md(ib), N, median(md));
    fprintf('L2-star      %.4f (median %.4f)\n', l2(ib), median(l2));
    fprintf('D* exact     %.4f train | %.4f train+test\n', dstar(Utr), dstar([Utr;Ute]));
    fprintf('LHS slice widths: %.3g / %.3g / %.3g -> set each knob within +/-%.2g / %.2g / %.2g\n', ...
        (cfg.ub-cfg.lb)/n, cfg.sliceWarn*(cfg.ub-cfg.lb)/n);
    fprintf('Run sheet: %s (%d fixed calls, %d adaptive)\n', cfg.planFile, size(X,1), cfg.budget - size(X,1));
end
end

function S = plancompare(cfg)
% Compare the chosen training plan (from bb_plan.csv) with other ways of
% placing the same number of points: random uniform, random LHS,
% correlation-minimising LHS, a single maximin LHS (no screening) and
% scrambled Halton (random start), each drawn nRep times. Plain (unscrambled) Halton is a
% single deterministic plan. Metrics: min pairwise distance (spread),
% L2-star and exact D* (uniformity), and empty 1D slices (stratification:
% for each input, how many of its n equal slices hold no point).
if ~exist(cfg.planFile,'file'), error('No run sheet: run bb_session(''plan'') first'); end
plan = dlmread(cfg.planFile, ',', 1, 0);
Uc = nrm(plan(plan(:,2)==1, 3:5), cfg);
n = size(Uc,1); nRep = 100;
rng(cfg.seed + 1, 'twister');
names = {'Random','Random LHS','Corr. LHS','Maximin LHS','Scr. Halton'};
gen = {@() rand(n,3), ...
       @() lhsdesign(n,3,'criterion','none'), ...
       @() lhsdesign(n,3,'criterion','correlation','iterations',cfg.lhsIter), ...
       @() lhsdesign(n,3,'criterion','maximin','iterations',cfg.lhsIter), ...
       @() net(scramble(haltonset(3,'Skip',randi(1e4)),'RR2'), n)};   % RR2 is a fixed
       % permutation, so a random start (skip) is what makes each draw different
nG = numel(gen); M = nan(nRep, nG, 4); P = cell(nRep, nG);
for j = 1:nG
    for i = 1:nRep
        P{i,j} = gen{j}(); M(i,j,:) = plan_metrics(P{i,j});
    end
end
mc = plan_metrics(Uc);
Uh = net(haltonset(3,'Skip',1), n); mh = plan_metrics(Uh);

mnames = {'min distance','L2-star','D* exact','empty 1D slices'};
fprintf('\n=== Sampling plan comparison: n = %d, d = 3, %d draws per method (medians) ===\n', n, nRep);
fprintf('%-24s %12s %10s %10s %16s\n', 'plan', 'min dist', 'L2-star', 'D*', 'empty slices');
for j = 1:nG
    fprintf('%-24s %12.3f %10.4f %10.4f %16.0f\n', names{j}, squeeze(median(M(:,j,:),1)));
end
fprintf('%-24s %12.3f %10.4f %10.4f %16.0f\n', 'Halton (unscrambled)', mh);
fprintf('%-24s %12.3f %10.4f %10.4f %16.0f\n', 'CHOSEN (screened)', mc);
fprintf('Chosen plan beats %.0f%% of single maximin LHS and %.0f%% of scrambled Halton on min distance,\n', ...
        100*mean(M(:,4,1) < mc(1)), 100*mean(M(:,5,1) < mc(1)));
fprintf('and %.0f%% / %.0f%% on L2-star.\n', 100*mean(M(:,4,2) > mc(2)), 100*mean(M(:,5,2) > mc(2)));

% ---- figure: top = x1-x2 projections, bottom = metric distributions
ink = [0.12 0.12 0.12]; muted = [0.45 0.45 0.45]; faint = [0.88 0.88 0.88];
blue = [42 120 214]/255; orange = [235 104 52]/255; green = [27 175 122]/255;
f = figure('Color','w', 'Position', [80 80 1300 720]);
tl = tiledlayout(f, 2, 4, 'TileSpacing','compact', 'Padding','compact');
[~, iMed] = min(abs(M(:,:,1) - median(M(:,:,1),1)), [], 1);   % representative (median min-distance) draw
show = {Uc, 'Maximin LHS, screened (chosen)', orange; ...
        P{iMed(5),5}, 'Scrambled Halton', blue; ...
        P{iMed(2),2}, 'Random LHS', blue; ...
        P{iMed(1),1}, 'Random uniform', blue};
lbl = 'abcdefgh';
g = (0:n)/n; gl = [g; g; nan(1,n+1)];
for p = 1:4
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
    plot(ax, gl(:), repmat([0 1 NaN]', n+1, 1), '-', 'Color', faint);   % slice lines (x)
    plot(ax, repmat([0 1 NaN]', n+1, 1), gl(:), '-', 'Color', faint);   % slice lines (y)
    U = show{p,1};
    plot(ax, U(:,1), U(:,2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', show{p,3}, 'MarkerEdgeColor','w');
    axis(ax, 'square'); xlim(ax, [0 1]); ylim(ax, [0 1]);
    m = plan_metrics(U);
    title(ax, sprintf('(%s) %s', lbl(p), show{p,2}), 'FontWeight','normal');
    subtitle(ax, sprintf('min d %.3f | L2* %.4f | empty slices %d', m(1), m(2), m(4)), 'Color', muted);
    xlabel(ax, 'x1 (normalised)'); ylabel(ax, 'x2 (normalised)');
    set(ax, 'TickDir','out', 'XColor', muted, 'YColor', muted, 'Layer','top');
end
better = {'higher = better','lower = better','lower = better','lower = better'};
for k = 1:4
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    x = repmat(1:nG, nRep, 1);
    boxchart(ax, x(:), reshape(M(:,:,k), [], 1), 'BoxFaceColor', blue, 'MarkerColor', muted);
    yline(ax, mh(k), '--', 'Halton', 'Color', green, 'LineWidth', 1.2, ...
          'LabelHorizontalAlignment','left', 'FontSize', 8);
    yline(ax, mc(k), '-', 'chosen', 'Color', orange, 'LineWidth', 1.8, ...
          'LabelHorizontalAlignment','right', 'FontSize', 8);
    set(ax, 'XTick', 1:nG, 'XTickLabel', names, 'TickDir','out', ...
        'XColor', muted, 'YColor', muted, 'GridAlpha', 0.15);
    xtickangle(ax, 30); xlim(ax, [0.4 nG+0.6]);
    title(ax, sprintf('(%s) %s', lbl(4+k), mnames{k}), 'FontWeight','normal');
    subtitle(ax, better{k}, 'Color', muted);
end
title(tl, sprintf('Sampling plans for n = %d points in 3 inputs (%d draws per method)', n, nRep), ...
      'FontWeight','normal', 'Color', ink);
print(f, 'fig_plan_compare.png', '-dpng', '-r300');
fprintf('Saved fig_plan_compare.png\n');
S = struct('names', {names}, 'metrics', M, 'metricNames', {mnames}, 'chosen', mc, 'halton', mh);
end

function m = plan_metrics(U)
% [min pairwise distance, L2-star, exact D*, empty 1D slices summed over inputs]
n = size(U,1);
b = min(floor(U*n), n-1);
emp = 0; for k = 1:size(U,2), emp = emp + n - numel(unique(b(:,k))); end
m = [min(pdist(U)), l2star(U), dstar(U), emp];
end

function D = l2star(U)
% L2-star discrepancy, Warnock closed form (fast proxy for D*)
[n,d] = size(U); P = ones(n);
for k = 1:d, P = P .* (1 - max(U(:,k), U(:,k)')); end
D = sqrt(max(3^-d - 2^(1-d)/n*sum(prod(1-U.^2,2)) + sum(P(:))/n^2, 0));
end

function D = dstar(U)
% exact star discrepancy: every anchored box with corners on the point grid
n = size(U,1); D = 0;
g1 = unique([U(:,1);1]); g2 = unique([U(:,2);1]); g3 = unique([U(:,3);1]);
for a = g1'
    o1 = U(:,1) < a; c1 = U(:,1) <= a;
    for b = g2'
        o2 = o1 & U(:,2) < b; c2 = c1 & U(:,2) <= b;
        for c = g3'
            v = a*b*c;
            D = max([D, sum(c2 & U(:,3) <= c)/n - v, v - sum(o2 & U(:,3) < c)/n]);
        end
    end
end
end

% =====================================================================
%  SURROGATES
% =====================================================================
function F = fit_best(L, cfg)
F = compare(L, cfg, [], false);
end

function F = compare(L, cfg, Lte, verbose)
% Merge points closer than the setting precision (they would make the
% interpolation matrices near-singular), normalise to [0,1], compare by LOOCV.
[~,~,j] = unique(round(L(:,3:5) ./ cfg.res), 'rows');
X = [accumarray(j, L(:,3), [], @mean), accumarray(j, L(:,4), [], @mean), accumarray(j, L(:,5), [], @mean)];
y = accumarray(j, L(:,6), [], @mean);
U = nrm(X, cfg); n = numel(y);
R = {};   % name, nrmse, yloo, model

% --- polynomials (fitlm, exact LOO from leverage)
for i = 1:numel(cfg.polyForms)
    [e, yl, m] = loo_poly(U, y, cfg.polyForms{i});
    if isfinite(e), R(end+1,:) = {['poly ' cfg.polyForms{i}], e, yl, m}; end %#ok<AGROW>
end
% --- RBF (own code, exact LOO from Rippa's formula, width tuned)
for k = 1:numel(cfg.kernels)
    ker = cfg.kernels{k};
    if any(strcmp(ker, {'cubic','thinplate'})), W = 1; else, W = cfg.widths; end
    bestE = inf;
    for w = W
        s = struct('kernel',ker,'width',w,'lambda',cfg.lambda);
        [e, yl] = loo_rbf(U, y, s);
        if e < bestE, bestE = e; bs = s; byl = yl; end
    end
    if isinf(bestE), continue; end
    if numel(W) > 1, nm = sprintf('RBF %s s=%.3f', ker, bs.width); else, nm = ['RBF ' ker]; end
    R(end+1,:) = {nm, bestE, byl, fit_rbf(U, y, bs)}; %#ok<AGROW>
end
% --- Kriging (fitrgp, LOO by refitting)
if cfg.useGP && exist('fitrgp', 'file')
    try
        [e, yl, m] = loo_gp(U, y);
        R(end+1,:) = {'Kriging (ARD SE)', e, yl, m};
    catch err
        if verbose, fprintf('(Kriging skipped: %s)\n', err.message); end
    end
end

[~,o] = sort(cell2mat(R(:,2))); R = R(o,:);
F.name = R{1,1}; F.nrmse = R{1,2}; F.q2 = q2(R{1,2}, n);
F.rmse = F.nrmse*std(y); F.y = y; F.yloo = R{1,3}; F.m = R{1,4};

if verbose
    hasT = ~isempty(Lte);
    fprintf('%-30s %9s %8s', 'model', 'LOO NRMSE', 'LOO Q2');
    if hasT, fprintf(' %10s', 'test NRMSE'); end; fprintf('\n');
    for i = 1:size(R,1)
        fprintf('%-30s %9.3f %8.3f', R{i,1}, R{i,2}, q2(R{i,2}, n));
        if hasT
            fprintf(' %10.3f', sqrt(mean((Lte(:,6) - spredict(R{i,4}, nrm(Lte(:,3:5),cfg))).^2))/std(y));
        end
        fprintf('\n');
    end
    fprintf('Selected: %s\n', F.name);
end
end

function [e, yl, m] = loo_poly(U, y, form)
% fitlm + exact LOO residuals e_i = r_i / (1 - h_ii)
m = []; yl = nan(size(y)); e = inf;
nTerms = struct('linear',4,'quadratic',10,'poly333',20);
if nTerms.(form) >= numel(y) - 1, return; end
ws = warning('off', 'all');
mdl = fitlm(U, y, form);
warning(ws);
h = mdl.Diagnostics.Leverage; r = mdl.Residuals.Raw;
if any(h > 1 - 1e-8), return; end
el = r ./ (1 - h);
e = sqrt(mean(el.^2)) / std(y); yl = y - el;
m = struct('kind','lm','mdl',mdl);
end

function [e, yl, m] = loo_gp(U, y)
% Kriging, ARD squared-exponential kernel, constant trend. Deterministic
% box: noise fixed near zero, so the model (nearly) interpolates.
s0  = 1e-3*std(y);
mdl = fitrgp(U, y, 'KernelFunction','ardsquaredexponential', 'BasisFunction','constant', ...
             'Standardize',false, 'Sigma',s0, 'SigmaLowerBound',s0/10, 'ConstantSigma',true);
cv  = crossval(mdl, 'Leaveout', 'on');
yl  = kfoldPredict(cv);
e   = sqrt(mean((y - yl).^2)) / std(y);
m   = struct('kind','gp','mdl',mdl);
end

function [e, yl] = loo_rbf(U, y, s)
% Exact LOOCV for RBF with linear tail (Rippa 1999): e_i = a_i / (K^-1)_ii
n = numel(y); z = (y - mean(y))/std(y);
K = rbf_system(U, s);
if rcond(K) < 1e-15, e = inf; yl = nan(n,1); return; end
Ki = inv(K); a = Ki * [z; zeros(4,1)];
ez = a(1:n) ./ diag(Ki(1:n,1:n));
e = sqrt(mean(ez.^2)); yl = y - ez*std(y);
end

function m = fit_rbf(U, y, s)
m.kind = 'rbf'; m.s = s; m.mu = mean(y); m.sd = std(y); m.C = U;
a = rbf_system(U, s) \ [(y - m.mu)/m.sd; zeros(4,1)];
m.w = a(1:end-4); m.c = a(end-3:end);
end

function K = rbf_system(U, s)
n = size(U,1); P = [ones(n,1) U];
K = [kern(pdist2(U,U), s) + s.lambda*eye(n), P; P', zeros(4)];
end

function Phi = kern(R, s)
switch s.kernel
    case 'gaussian',        Phi = exp(-(R/s.width).^2);
    case 'multiquadric',    Phi = sqrt(1 + (R/s.width).^2);
    case 'invmultiquadric', Phi = 1 ./ sqrt(1 + (R/s.width).^2);
    case 'cubic',           Phi = R.^3;
    case 'thinplate',       Phi = R.^2 .* log(R + (R==0));
end
end

function y = spredict(m, U)
% one predict for every model type (named spredict so the toolbox
% predict() for fitlm / fitrgp objects is not shadowed)
switch m.kind
    case {'lm','gp'}, y = predict(m.mdl, U);
    case 'rbf'
        y = m.mu + m.sd * (kern(pdist2(U, m.C), m.s) * m.w + [ones(size(U,1),1) U] * m.c);
end
end

function [uB, yB, peaks] = maximise(m, nGrid, nSeed, method)
% grid -> up to nSeed separated seeds -> fmincon SQP (bounded) from each.
% The local runs are done by MultiStart (Global Optimization Toolbox) with
% the grid seeds as custom start points, so the result stays deterministic
% and the grid still guarantees global coverage. method = 'loop' runs the
% same fmincon calls by hand (also the fallback without that toolbox).
% peaks = [u y] of distinct local maxima, best first.
% Defaults (21^3 grid, 6 seeds) are justified by bb_session('maxcheck').
if nargin < 2, nGrid = 21; end
if nargin < 3, nSeed = 6; end
if nargin < 4, method = 'multistart'; end
if strcmp(method, 'multistart') && ~(license('test','GADS_Toolbox') && exist('MultiStart','file'))
    method = 'loop';
end
g = linspace(0,1,nGrid); [a,b,c] = ndgrid(g,g,g); G = [a(:) b(:) c(:)];
[~,o] = sort(spredict(m, G), 'descend'); S = G(o(1),:);
for i = o(2:end)'
    if size(S,1) >= nSeed, break; end
    if min(pdist2(G(i,:), S)) > 0.25, S(end+1,:) = G(i,:); end %#ok<AGROW>
end
opts = optimoptions('fmincon', 'Algorithm','sqp', 'Display','off');
switch method
    case 'multistart'
        problem = createOptimProblem('fmincon', 'objective', @(u) -spredict(m,u), ...
                  'x0', S(1,:), 'lb', [0 0 0], 'ub', [1 1 1], 'options', opts);
        ms = MultiStart('Display','off', 'XTolerance',0.15);
        [~, ~, ~, ~, sols] = run(ms, problem, CustomStartPointSet(S));
        U = vertcat(sols.X);                     % local solutions, best first
    case 'loop'
        U = zeros(size(S));
        for i = 1:size(S,1)
            U(i,:) = fmincon(@(u) -spredict(m,u), S(i,:), [],[],[],[], [0 0 0], [1 1 1], [], opts);
        end
    otherwise, error('method must be multistart | loop');
end
U = min(max(U, 0), 1);
% keep distinct peaks only (same 0.15 rule for both methods, so they are comparable)
peaks = zeros(0,4);
for i = 1:size(U,1)
    if isempty(peaks) || min(pdist2(U(i,:), peaks(:,1:3))) > 0.15
        peaks(end+1,:) = [U(i,:) spredict(m,U(i,:))]; %#ok<AGROW>
    end
end
[~,o] = sort(peaks(:,4), 'descend'); peaks = peaks(o,:);
uB = peaks(1,1:3); yB = peaks(1,4);
end

function T = maxcheck(cfg)
% Sensitivity of the surrogate optimum to the inner-optimiser settings.
% Fits the final model to the current log once, then reruns maximise for
% every combination of grid size and seed count. Each optimum is compared
% with the finest setting (41^3 grid, 12 seeds) in knob setting steps
% (cfg.res): <= 1 step means the choice makes no practical difference.
% T rows: [nGrid nSeed x1 x2 x3 yPred dy maxSteps nPeaks time_s]
L = load_log(cfg);
if size(L,1) < 10, error('maxcheck needs at least 10 logged calls'); end
F = fit_best(L, cfg);
grids = [11 21 41]; seeds = [3 6 12];
[uR, yR] = maximise(F.m, grids(end), seeds(end));
xR = dn(uR, cfg);

fprintf('\n=== Inner optimiser sensitivity (%d calls, model %s) ===\n', size(L,1), F.name);
if license('test','GADS_Toolbox'), fprintf('Rows use MultiStart (fmincon SQP from the grid seeds).\n');
else, fprintf('Global Optimization Toolbox missing: rows use the fmincon loop.\n'); end
if strcmp(F.m.kind, 'rbf') && any(strcmp(F.m.s.kernel, {'gaussian','multiquadric','invmultiquadric'}))
    fprintf('RBF width %.3f vs grid spacing %.3f (21^3)', F.m.s.width, 1/20);
    if F.m.s.width < 2/20, fprintf('  !! narrow basis: refine grid or raise min width'); end
    fprintf('\n');
end
fprintf('Reference (41^3, 12 seeds): x = [%g %g %g], predicted MCSI %.6g\n\n', xR, yR);
fprintf('%5s %5s   %-26s %10s %10s %9s %6s %7s\n', 'grid', 'seeds', 'optimum x', ...
        'pred MCSI', 'dy', 'max steps', 'peaks', 'time s');
T = zeros(0,10);
for g = grids
    for s = seeds
        t0 = tic; [u, y, pk] = maximise(F.m, g, s); t = toc(t0);
        x = dn(u, cfg);
        steps = max(abs(x - xR) ./ cfg.res);
        T(end+1,:) = [g s x y y-yR steps size(pk,1) t]; %#ok<AGROW>
        fprintf('%4d^3 %5d   [%7.4g %7.4g %7.4g] %10.6g %+10.2e %9.0f %6d %7.2f\n', ...
                g, s, x, y, y-yR, steps, size(pk,1), t);
    end
end
t0 = tic; [u, y, pk] = maximise(F.m, 21, 6, 'loop'); t = toc(t0);
x = dn(u, cfg); steps = max(abs(x - xR) ./ cfg.res);
fprintf('%-11s   [%7.4g %7.4g %7.4g] %10.6g %+10.2e %9.0f %6d %7.2f   <- 21^3, 6 seeds, fmincon loop (no MultiStart)\n', ...
        'hand loop', x, y, y-yR, steps, size(pk,1), t);
T(end+1,:) = [21 6 x y y-yR steps size(pk,1) t];
ok = all(T(:,8) <= 1);
r = find(T(:,1) == 21 & T(:,2) == 6, 1);
fprintf('\nDefault (21^3, 6 seeds): %.0f setting step(s) from reference, dy = %+.2e (%.4f%% of observed range)\n', ...
        T(r,8), T(r,7), 100*abs(T(r,7))/range_(L(:,6)));
if ok
    fprintf('All settings agree to within 1 setting step: optimum insensitive to grid and seed count.\n');
else
    fprintf('!! Some settings disagree by more than 1 setting step: check rows above.\n');
end
end

% =====================================================================
%  REHEARSAL AND SELF-TEST
% =====================================================================
function out = rehearse(cfg, kind)
if nargin < 2, kind = 'easy'; end
if exist(cfg.logFile,'file'), delete(cfg.logFile); end
make_plan(cfg, true);
t0 = tic; run_loop(cfg, @(x) fake_box(x, kind), true); t = toc(t0);
L = load_log(cfg);
[xt, yt] = true_opt(cfg, kind);
ib = best_obs(L);
fprintf('\nTRUE optimum : [%.4g %.4g %.4g]  MCSI %.5g\n', xt, yt);
fprintf('Best observed: [%.4g %.4g %.4g]  MCSI %.5g  (%.3f%% below)\n', ...
        L(ib,3:5), L(ib,6), 100*(yt-L(ib,6))/yt);
fprintf('Phase counts (LHS hold rep exploit explore local confirm): %s\n', ...
        mat2str(accumarray(L(:,2), 1, [7 1])'));
fprintf('Compute time %.1f s\n', t);
out = L;
end

function y = fake_box(x, kind)
% Rehearsal stand-in for the real box: smooth response, shown like the real
% display (clipped to 0-10, rounded to 4 decimals).
if nargin < 2, kind = 'easy'; end
y = min(max(round(fake_true(x, kind) * 1e4) / 1e4, 0), 10);
end

function y = fake_true(x, kind)
% Smooth (unrounded) fake responses, used for the true optimum because
% fmincon cannot differentiate the rounded display values.
%  'easy': one smooth interior peak, diminishing returns, skew/targeting
%          interaction. Optimum 9.62 at about [129 26 0.48].
%  'hard': stress test with two competing strategies.
%          broad    = heavy advertising, no targeting, low skew: wide, forgiving
%                     peak, MCSI 8.58 at about [117 0 0.25] (on the x2 = 0 bound)
%          targeted = targeting spend with high skew: narrow in skew (width
%                     0.06) and its best skew moves with targeting. GLOBAL
%                     optimum 9.00 at about [117 24 0.71], 4.7% above broad.
%          Advertising needs critical mass (S-curve around 24 £k) before
%          diminishing returns; spend costs; baseline index about 0.1-2.
%          Only ~0.4% of the input space beats the broad peak, and ~54% of
%          random local searches end on it: deliberately deceptive.
if nargin < 2, kind = 'easy'; end
a = x(:,1)/200; b = x(:,2)/60; c = x(:,3);
switch kind
    case 'easy'
        y = 10*((1-exp(-3*a)).*(0.6+0.4*(1-exp(-5*b))).*(1+0.5*exp(-((c-0.35-0.3*b)/0.22).^2)) ...
             - 0.45*a.^2 - 0.25*b.^2 - 0.1*a.*b);
    case 'hard'
        reach    = (1 - exp(-4*a)) ./ (1 + exp(-(a - 0.12)/0.04));
        broad    = exp(-((c - 0.25)/0.15).^2) .* (1 - 0.5*b);
        targeted = exp(-((c - 0.65 - 0.15*b)/0.06).^2) .* (1 - exp(-6*b));
        cost     = 0.3*a.^2 + 0.25*b.^2 + 0.1*a.*b;
        y = 8.839 * (reach .* (0.55 + 0.35*broad + 0.52*targeted) - cost + 0.26);
    otherwise, error('fake box kind must be easy | hard');
end
end

function [xt, yt, yr] = true_opt(cfg, kind)
% True optimum (and best rival peak yr, -Inf if none) of fake_true by
% multistart fmincon. Rehearsal and study only.
if nargin < 2, kind = 'easy'; end
rs = rng; rng(0, 'twister');
opts = optimoptions('fmincon', 'Algorithm','sqp', 'Display','off');
pk = zeros(0,4);
for s = 1:60
    u = fmincon(@(u) -fake_true(dn(u,cfg,false), kind), rand(1,3), [],[],[],[], [0 0 0],[1 1 1],[],opts);
    v = fake_true(dn(u,cfg,false), kind);
    if isempty(pk) || min(pdist2(u, pk(:,1:3))) > 0.15, pk(end+1,:) = [u v]; %#ok<AGROW>
    else, [~,i] = min(pdist2(u, pk(:,1:3))); if v > pk(i,4), pk(i,:) = [u v]; end
    end
end
rng(rs);
[~,o] = sort(pk(:,4), 'descend'); pk = pk(o,:);
xt = dn(pk(1,1:3), cfg, false); yt = pk(1,4);
if size(pk,1) > 1, yr = pk(2,4); else, yr = -Inf; end
end

function R = study(cfg, nTr, seeds, kinds)
% How many LHS training points? Full rehearsals of the whole session on the
% fake boxes for each nTrain (budget, hold-out and confirm calls unchanged),
% repeated over seeds (different plan + knob-setting errors). Measures:
%   gap of the best observed reading below the true optimum (%)
%   gap of the RECOMMENDED setting (true value there) below the optimum (%)
%   call at which the session first went local (optimum judged found)
%   success: recommendation better than the rival peak (found the global one)
% Runs on study_* files, so the real bb_log.csv is never touched.
% R rows: [kind nTrain seed bestGap recGap firstLocal success time]
if nargin < 2 || isempty(nTr), nTr = [10 15 20 25 30 35]; end
if nargin < 3 || isempty(seeds), seeds = 1:5; end
if nargin < 4 || isempty(kinds), kinds = {'easy','hard'}; end
if ischar(kinds), kinds = {kinds}; end
c0 = cfg; c0.logFile = 'study_log.csv'; c0.planFile = 'study_plan.csv'; c0.stateFile = 'study_state.mat';
nK = numel(kinds); nI = numel(nTr); nS = numel(seeds); B = cfg.budget;
R = nan(0, 8); G = nan(nK, nI, nS, B);   % G: best-so-far gap (%) after each call
yrs = nan(1, nK); yts = nan(1, nK);      % rival and global peak value per box
for q = 1:nK
    [~, yt, yr] = true_opt(cfg, kinds{q}); yrs(q) = yr; yts(q) = yt;
    fprintf('\n=== nTrain study, %s box: %d settings x %d seeds, budget %d, optimum %.5g, rival %.5g ===\n', ...
            kinds{q}, nI, nS, B, yt, yr);
    fprintf('%7s %5s %12s %12s %12s %8s %8s\n', 'nTrain', 'seed', 'best gap %', 'rec. gap %', 'first local', 'global?', 'time s');
    for i = 1:nI
        for j = 1:nS
            c = c0; c.nTrain = nTr(i); c.seed = 100 + seeds(j);
            cleanup(c); t0 = tic;
            make_plan(c, false); run_loop(c, @(x) fake_box(x, kinds{q}), false);
            L = load_log(c);
            G(q,i,j,1:size(L,1)) = 100*(yt - cummax(L(:,6)))/yt;
            F = fit_best(L, c); xr = dn(maximise(F.m), c); yrec = fake_true(xr, kinds{q});
            fl = find(L(:,2) == 6, 1); if isempty(fl), fl = NaN; end
            R(end+1,:) = [q, nTr(i), seeds(j), G(q,i,j,end), 100*(yt - yrec)/yt, fl, yrec > yr, toc(t0)]; %#ok<AGROW>
            fprintf('%7d %5d %12.4f %12.4f %12g %8d %8.1f\n', R(end,2:8));
        end
    end
end
cleanup(c0);

fprintf('\n%-5s %7s %13s %13s %14s %15s %9s\n', 'box', 'nTrain', 'best gap med', 'rec gap med', 'rec gap worst', 'local at (med)', 'global');
for q = 1:nK
    for i = 1:nI
        r = R(R(:,1) == q & R(:,2) == nTr(i), :);
        fl = r(:,6); nNever = sum(isnan(fl)); fl(isnan(fl)) = B + 1;
        if median(fl) > B, ls = 'never'; else, ls = sprintf('%g', median(fl)); end
        if nNever > 0, ls = sprintf('%s (%d never)', ls, nNever); end
        fprintf('%-5s %7d %13.4f %13.4f %14.4f %15s %6d/%d\n', kinds{q}, nTr(i), median(r(:,4)), ...
                median(r(:,5)), max(r(:,5)), ls, sum(r(:,7)), size(r,1));
    end
end

% ---- figure: one row per fake box
ink = [0.12 0.12 0.12]; muted = [0.45 0.45 0.45]; orange = [235 104 52]/255;
cm = interp1([0 1], [158 197 244; 13 54 107]/255, linspace(0,1,nI));   % light -> dark = more training points
f = figure('Color','w', 'Position', [60 60 1500 420*nK]);
tl = tiledlayout(f, nK, 4, 'TileSpacing','compact', 'Padding','compact');
floorGap = 1e-4;                       % plot floor for the log axes (gap 0 = optimum hit)
for q = 1:nK
    rq = R(R(:,1) == q, :); L0 = upper(kinds{q}(1)); kn = [L0 kinds{q}(2:end)];
    gapR = 100*(yts(q) - yrs(q))/yts(q);   % rival peak, % below the optimum
    % (a) progress
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    h = gobjects(nI,1);
    for i = 1:nI
        g = squeeze(median(G(q,i,:,:), 3));
        h(i) = plot(ax, 1:B, max(g, floorGap), '-', 'Color', cm(i,:), 'LineWidth', 1.8);
        xline(ax, nTr(i) + cfg.nTest + 0.5, ':', 'Color', cm(i,:), 'LineWidth', 1);
    end
    set(ax, 'YScale','log'); ylim(ax, [1e-3 100]);
    if isfinite(yrs(q))
        yline(ax, gapR, '--', 'rival peak', 'Color', orange, 'LabelHorizontalAlignment','left', 'FontSize', 8);
    end
    xlabel(ax, 'call number'); ylabel(ax, 'best observed, % below optimum');
    title(ax, sprintf('%s box: (a) search progress (median)', kn), 'FontWeight','normal');
    subtitle(ax, 'dotted: end of fixed block (LHS + hold-out)', 'Color', muted);
    if q == 1, legend(ax, h, compose('nTrain = %d', nTr), 'Location','northeast', 'Box','off'); end
    % (b) recommendation gap, log axis
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for i = 1:nI
        r = rq(rq(:,2) == nTr(i), :);
        plot(ax, nTr(i) + 0.8*(rand(size(r,1),1) - 0.5), max(r(:,5), floorGap), 'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', cm(i,:), 'MarkerEdgeColor','w');
    end
    med = arrayfun(@(n) median(rq(rq(:,2)==n, 5)), nTr);
    plot(ax, nTr, max(med, floorGap), '-', 'Color', ink, 'LineWidth', 1.5);
    set(ax, 'YScale','log'); ylim(ax, [floorGap/2 20]);
    if isfinite(yrs(q))
        yline(ax, gapR, '--', 'rival peak', 'Color', orange, 'LabelHorizontalAlignment','left', 'FontSize', 8);
    end
    xline(ax, cfg.nTrain, '-', sprintf('chosen %d', cfg.nTrain), 'Color', orange, 'LineWidth', 1.5, ...
          'LabelVerticalAlignment','bottom', 'FontSize', 8);
    xlabel(ax, 'nTrain'); ylabel(ax, 'true MCSI at recommendation, % below optimum');
    title(ax, '(b) Quality of the recommendation', 'FontWeight','normal');
    subtitle(ax, sprintf('dots = seeds, line = median (0 plotted at %g)', floorGap), 'Color', muted);
    % (c) first local call, 'never' shown above the budget
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    yNever = B + 3;
    for i = 1:nI
        r = rq(rq(:,2) == nTr(i), 6); jit = nTr(i) + 0.8*(rand(size(r)) - 0.5);
        ok = ~isnan(r);
        plot(ax, jit(ok), r(ok), 'o', 'MarkerSize', 6, 'MarkerFaceColor', cm(i,:), 'MarkerEdgeColor','w');
        plot(ax, jit(~ok), yNever*ones(sum(~ok),1), 'x', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', cm(i,:));
    end
    medL = arrayfun(@(n) median(fillmissing(rq(rq(:,2)==n, 6), 'constant', yNever)), nTr);
    plot(ax, nTr, min(medL, yNever), '-', 'Color', ink, 'LineWidth', 1.5);
    plot(ax, nTr, nTr + cfg.nTest, '--', 'Color', muted);
    yline(ax, B + 0.5, '-', 'Color', muted);
    xline(ax, cfg.nTrain, '-', 'Color', orange, 'LineWidth', 1.5);
    ylim(ax, [min(nTr) + cfg.nTest - 2, yNever + 2]);
    yt_ = get(ax, 'YTick'); yt_ = [yt_(yt_ <= B) yNever];
    set(ax, 'YTick', yt_, 'YTickLabel', [compose('%g', yt_(1:end-1)) {'never'}]);
    xlabel(ax, 'nTrain'); ylabel(ax, 'first local call');
    title(ax, '(c) When the optimum was judged found', 'FontWeight','normal');
    subtitle(ax, 'x = never went local; dashed = end of fixed block', 'Color', muted);
    % (d) success rate
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    sr = arrayfun(@(n) 100*mean(rq(rq(:,2)==n, 7)), nTr);
    for i = 1:nI
        bar(ax, nTr(i), sr(i), 3.5, 'FaceColor', cm(i,:), 'EdgeColor','none');
        text(ax, nTr(i), sr(i) + 3, sprintf('%d/%d', round(sr(i)*nS/100), nS), ...
             'HorizontalAlignment','center', 'FontSize', 8, 'Color', ink);
    end
    xline(ax, cfg.nTrain, '-', 'Color', orange, 'LineWidth', 1.5);
    ylim(ax, [0 112]); xlabel(ax, 'nTrain'); ylabel(ax, 'runs ending on the global peak (%)');
    title(ax, '(d) Found the global optimum', 'FontWeight','normal');
    if isinf(yrs(q)), st = 'single peak: always 100%';
    else, st = 'recommendation better than the rival peak'; end
    subtitle(ax, st, 'Color', muted);
    for ax = findobj(f, 'Type','axes')'
        set(ax, 'TickDir','out', 'XColor', muted, 'YColor', muted, 'GridAlpha', 0.15);
    end
    set(findobj(f, 'Type','axes'), 'XTick', nTr, 'XLim', [min(nTr)-3, max(nTr)+3]);
end
% (a) panels use call number on x, not nTrain
for ax = tl.Children'
    if isa(ax, 'matlab.graphics.axis.Axes') && strcmp(ax.XLabel.String, 'call number')
        set(ax, 'XTickMode','auto'); xlim(ax, [1 B]);
    end
end
for ax = findobj(f, 'Type','axes')'
    try ax.Toolbar.Visible = 'off'; catch, end
end
title(tl, sprintf('Effect of initial sample size (budget %d, %d hold-out, %d seeds per setting)', ...
      B, cfg.nTest, nS), 'FontWeight','normal', 'Color', ink);
print(f, 'fig_ntrain_study.png', '-dpng', '-r300');
fprintf('Saved fig_ntrain_study.png\n');
end

function cleanup(c)
for fn = {c.logFile, c.planFile, c.stateFile}
    if exist(fn{1}, 'file'), delete(fn{1}); end
end
end

function worst = selftest(cfg)
need = {'lhsdesign','pdist','pdist2','fitlm','fitrgp','quantile','haltonset','fmincon','optimoptions','MultiStart'};
for i = 1:numel(need)
    if exist(need{i}, 'file'), st = 'found'; else, st = 'MISSING'; end
    fprintf('%-13s %s\n', need{i}, st);
end
rng(1, 'twister'); U = rand(30,3); y = fake_true(dn(U,cfg,false)); worst = 0;
% RBF closed-form vs brute force
for s = {struct('kernel','gaussian','width',0.4,'lambda',1e-10), ...
         struct('kernel','multiquadric','width',0.4,'lambda',1e-3), ...
         struct('kernel','thinplate','width',1,'lambda',1e-10)}
    [~, yf] = loo_rbf(U, y, s{1}); yb = zeros(30,1);
    for i = 1:30
        j = true(30,1); j(i) = false; yb(i) = spredict(fit_rbf(U(j,:), y(j), s{1}), U(i,:));
    end
    d = max(abs(yf - yb))/range_(y); worst = max(worst, d);
    fprintf('RBF %-14s LOO fast vs brute: %.1e of range\n', s{1}.kernel, d);
end
% poly leverage LOO vs brute force
[~, yf] = loo_poly(U, y, 'quadratic'); yb = zeros(30,1);
for i = 1:30
    j = true(30,1); j(i) = false; [~,~,m] = loo_poly(U(j,:), y(j), 'quadratic');
    yb(i) = spredict(m, U(i,:));
end
d = max(abs(yf - yb))/range_(y); worst = max(worst, d);
fprintf('poly quadratic     LOO fast vs brute: %.1e of range\n', d);
end

% =====================================================================
%  SMALL HELPERS
% =====================================================================
function U = nrm(X, cfg), U = (X - cfg.lb) ./ (cfg.ub - cfg.lb); end
function X = dn(U, cfg, snap)
if nargin < 3, snap = true; end
X = cfg.lb + U .* (cfg.ub - cfg.lb);
if snap, X = min(max(round(X ./ cfg.res) .* cfg.res, cfg.lb), cfg.ub); end
end
function q = q2(nrmse, n), q = 1 - nrmse^2 * n/(n-1); end
function tf = is_sampled(x, L, cfg)
% already have data within half a setting step of x
tol = 0.5 * [1 0.5 0.01]; if nargin > 2, tol = 0.5 * cfg.res; end
tf = any(all(abs(L(:,3:5) - x) <= tol + 1e-12, 2));
end
function i = best_obs(L), [~, i] = max(L(:,6)); end
function r = range_(v), r = max(v) - min(v); if r == 0, r = 1; end, end
