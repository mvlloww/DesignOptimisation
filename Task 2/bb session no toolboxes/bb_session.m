function out = bb_session(cmd)
% BB_SESSION  Whole black-box workflow in one file.
%
%   bb_session('plan')      before the session: build + check sampling plan
%   bb_session              run / resume the session (interactive)
%   bb_session('status')    model table, decision checks, figures
%   bb_session('undo')      delete the last logged row (typo fix)
%   bb_session('rehearse')  full dry run on a fake function (deletes log!)
%   bb_session('selftest')  check fast LOOCV against brute force
%
% At the MCSI> prompt type the reading, or a command:
%   u = undo last   s = status   q = quit (resume later, nothing lost)
%   x / e / l / c = override this call: exploit / explore / local / confirm
%
% Everything is logged to bb_log.csv after every call. Restarting MATLAB
% and calling bb_session again resumes where you left off.
%
% Phase codes in the log: 1 LHS, 2 hold-out, 3 repeat, 4 exploit,
%                         5 explore, 6 local, 7 confirm

if nargin < 1, cmd = 'run'; end
cfg = config();
switch lower(cmd)
    case 'plan',     out = make_plan(cfg, true);
    case 'run',      run_loop(cfg, [], true);
    case 'status',   L = load_log(cfg); status(L, cfg, true);
    case 'undo',     undo_last(cfg);
    case 'rehearse', out = rehearse(cfg);
    case 'selftest', out = selftest(cfg);
    otherwise, error('unknown command %s', cmd);
end
end

% =====================================================================
%  CONFIG
% =====================================================================
function cfg = config()
cfg.lb    = [0   0  0];
cfg.ub    = [200 60 1];
cfg.res   = [1   0.5 0.01];        % dial step, real units (CHECK ON THE DAY)
cfg.names = {'x1 advert (£k)','x2 targeting (£k)','x3 skew (-)'};

cfg.budget  = 50;  cfg.nTrain = 25;  cfg.nTest = 5;  cfg.nConfirm = 3;
cfg.nCand   = 2000; cfg.seed = 27;

cfg.polyOrders = [1 2 3];
cfg.kernels    = {'gaussian','multiquadric','invmultiquadric','cubic','thinplate'};
cfg.widths     = logspace(-1.3, 0.5, 20);
cfg.lambdaDet  = 1e-10;                 % deterministic box
cfg.lambdaNoisy= [1e-6 1e-3 1e-2 1e-1]; % used automatically if repeats differ

% decision rules
cfg.stableTol   = 0.02;  % optimum moved < 2% of each range
cfg.errTol      = 0.01;  % last exploit error < 1% of observed MCSI range
cfg.peakTol     = 0.05;  % rival peak within 5% of range = competitor
cfg.q2High      = 0.80;  % above this, exploit 2:1, else 1:1
cfg.nLocalMin   = 6;     % need >= 6 calls (plus confirm) to start local
cfg.nearOptFrac = 0.02;  % local box = region within 2% of range of max
cfg.minBoxSteps = 5;     % local box half-width >= 5 dial steps
cfg.confirmStep = 0.05;  % confirm steps of 5% of range

cfg.logFile  = 'bb_log.csv';
cfg.planFile = 'bb_plan.csv';
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
        fprintf('\nBudget of %d calls used.\n', cfg.budget);
        status(L, cfg, interactive);
        return;
    end
    if n == size(plan,1) && isempty(override) && verbose
        fprintf('\n########## FIXED BLOCK DONE: DECISION POINT ##########');
        status(L, cfg, false);          % train-vs-hold-out table before infill
    end
    P = next_point(L, plan, cfg, override, verbose);
    override = '';

    if ~isempty(box)
        y = box(P.x);
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
        y = str2double(strrep(s, ',', '.'));
        if isnan(y), fprintf('  not a number: "%s"\n', s); continue; end
    end
    record(L, P, y, cfg, verbose);
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
    if P.phase ~= 1 && size(tr,1) >= 10
        F = fit_best(tr, cfg);
        P.ypred = predict(F.m, nrm(P.x,cfg)); P.sig = F.rmse;
    end
    show(P, n, cfg, verbose);
    return;
end

% ---- adaptive ---------------------------------------------------------
F   = fit_best(L, cfg);
[uS, yS, peaks] = maximise(F.m);
D   = decision_checks(L, F, uS, yS, peaks, cfg);
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
        if min(dist(u, Udata)) < cfg.stableTol || is_sampled(dn(u,cfg), L)
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
        if min(dist(uS, Udata)) >= cfg.stableTol && ~is_sampled(dn(uS,cfg), L)
            u = uS; P.note = 'confirm: surrogate optimum';
        else
            [ib] = best_obs(L); ub = nrm(L(ib,3:5), cfg);
            k = mod(sum(L(:,2)==7), 3) + 1;
            for sgn = [1 -1]
                u = ub; u(k) = min(max(u(k) + sgn*cfg.confirmStep, 0), 1);
                if ~is_sampled(dn(u,cfg), L), break; end
            end
            P.note = sprintf('confirm: 5%% step in %s from best observed', cfg.names{k});
        end
end
if strcmp(mode, 'explore')
    P.phase = 5;
    C = rand(20000,3); yc = predict(F.m, C);
    keep = yc >= quantile_(yc, 0.75); C = C(keep,:);
    [~, k] = max(min(dist2(C, Udata), [], 2)); u = C(k,:);
end
P.x = dn(u, cfg);
P.ypred = predict(F.m, nrm(P.x, cfg)); P.sig = F.rmse;
P.model = F.name;
if verbose
    fprintf('\nModel %s (LOO Q2 %.3f) | checks: stable %d, accurate %d, no rival %d, calls %d -> %s\n', ...
        F.name, F.q2, D.stable, D.accurate, ~D.rival && D.explored, D.calls, upper(mode));
end
show(P, n, cfg, verbose);
end

function D = decision_checks(L, F, uS, yS, peaks, cfg)
n = size(L,1); yr = range_(L(:,6));
% 1 stable: optimum with vs without the latest point
Fp = fit_best(L(1:end-1,:), cfg);
uP = maximise(Fp.m);
D.stable = dist(uS, uP) < cfg.stableTol;
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

function [lo, hi] = local_box(m, uS, yS, L, cfg)
% near-optimal region: predicted MCSI within nearOptFrac*range of the max
thr = yS - cfg.nearOptFrac * range_(L(:,6));
minW = cfg.minBoxSteps * cfg.res ./ (cfg.ub - cfg.lb);
lo = uS; hi = uS;
for k = 1:3
    for t = 0.005:0.005:1
        u = uS; u(k) = uS(k) - t; if u(k) < 0 || predict(m,u) < thr, break; end; lo(k) = u(k);
    end
    for t = 0.005:0.005:1
        u = uS; u(k) = uS(k) + t; if u(k) > 1 || predict(m,u) < thr, break; end; hi(k) = u(k);
    end
end
lo = max(min(lo, uS - minW), 0);
hi = min(max(hi, uS + minW), 1);
end

function u = spacefill(lo, hi, Udata, L, cfg)
% most isolated point inside the box, measured in box-scaled units
w = hi - lo;
C = lo + rand(5000,3) .* w;
d = min(dist2(C ./ w, Udata ./ w), [], 2);   % all data counts, in box units
[~, o] = sort(d, 'descend');
for k = o'
    if ~is_sampled(dn(C(k,:),cfg), L), u = C(k,:); return; end
end
u = C(o(1),:);
end

% =====================================================================
%  LOGGING
% =====================================================================
function record(L, P, y, cfg, verbose)
row = [size(L,1)+1, P.phase, P.x, y, P.ypred];
write_log([L; row], cfg);
yr = range_([L(:,6); y]);
if verbose
    fprintf('  logged run %d: MCSI = %.6g', row(1), y);
    if ~isnan(P.ypred)
        fprintf('  (pred %.5g, error %+.1f%% of range)', P.ypred, 100*(y-P.ypred)/yr);
    end
    fprintf('\n');
end
if ~isnan(P.ypred) && ~isnan(P.sig) && abs(y - P.ypred) > 3*P.sig && verbose
    fprintf('  !! reading is %.1f LOO-RMSEs from prediction: CHECK FOR A TYPO (u = undo)\n', ...
        abs(y-P.ypred)/P.sig);
end
if P.phase == 3
    i = find(all(abs(L(:,3:5) - P.x) < 1e-9, 2), 1);
    dy = y - L(i,6);
    if verbose, fprintf('  repeat of run %d: diff = %.3g\n', L(i,1), dy); end
    if abs(dy) > 1e-9*max(1,abs(y)) && verbose
        fprintf('  !! box is NOISY: switching to regression RBFs (nugget) automatically\n');
    end
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
F = fit_best(L, cfg, true);
[uS, yS, peaks] = maximise(F.m);
fprintf('Surrogate optimum [%g %g %g], predicted %.6g\n', dn(uS,cfg), yS);
if size(peaks,1) > 1
    fprintf('Other peaks: '); fprintf('[%g %g %g] %.4g;  ', [dn(peaks(2:end,1:3),cfg) peaks(2:end,4)]'); fprintf('\n');
end
onB = uS < 0.01 | uS > 0.99;
if any(onB), fprintf('!! optimum on bound: %s\n', strjoin(cfg.names(onB), ', ')); end
r = L(:,2) >= 4 & ~isnan(L(:,7));
if any(r)
    fprintf('Predict-then-query errors (%% of range): %s\n', ...
        sprintf('%+.1f ', 100*(L(r,6)-L(r,7))/range_(L(:,6))));
end
if makePlots, figures(L, F, uS, cfg); end
end

function figures(L, F, uS, cfg)
try
    f1 = figure('Color','w'); hold on; box on; grid on;
    plot(F.y, F.yloo, 'ko', 'MarkerFaceColor', [0.2 0.4 0.8]);
    lim = [min([F.y; F.yloo]) max([F.y; F.yloo])]; plot(lim, lim, 'k--');
    xlabel('observed MCSI'); ylabel('LOOCV predicted MCSI');
    title(sprintf('%s: Q^2 = %.3f', F.name, F.q2));
    print(f1, 'fig_parity.png', '-dpng', '-r200');

    f2 = figure('Color','w'); pr = [1 2 3; 1 3 2; 2 3 1]; g = linspace(0,1,60);
    xs = dn(uS, cfg, false); X = L(:,3:5);
    for p = 1:3
        a = pr(p,1); b = pr(p,2); c = pr(p,3);
        [A,B] = meshgrid(g,g); Q = zeros(numel(A),3);
        Q(:,a) = A(:); Q(:,b) = B(:); Q(:,c) = uS(c);
        subplot(1,3,p); hold on; box on;
        contourf(cfg.lb(a)+A*(cfg.ub(a)-cfg.lb(a)), cfg.lb(b)+B*(cfg.ub(b)-cfg.lb(b)), ...
                 reshape(predict(F.m,Q), size(A)), 20, 'LineColor','none'); colorbar;
        sty = {'wo','ws','w^','r.','m.','y.','g.'};
        for ph = 1:7
            k = L(:,2)==ph; if any(k), plot(X(k,a), X(k,b), sty{ph}, 'MarkerSize', 6); end
        end
        plot(xs(a), xs(b), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
        xlabel(cfg.names{a}); ylabel(cfg.names{b});
        title(sprintf('slice %s = %.3g', cfg.names{c}, xs(c)));
    end
    print(f2, 'fig_slices.png', '-dpng', '-r200');
catch err
    fprintf('(figures skipped: %s)\n', err.message);
end
end

% =====================================================================
%  SAMPLING PLAN
% =====================================================================
function plan = make_plan(cfg, verbose)
rand('seed', cfg.seed); try rng(cfg.seed); catch, end %#ok<RAND>
n = cfg.nTrain; N = cfg.nCand;
ph = zeros(N,1); l2 = ph; C = cell(N,1);
for i = 1:N
    C{i} = lhs(n); ph(i) = phi_q(C{i}); l2(i) = l2star(C{i});
end
[~,a] = sort(ph); [~,b] = sort(l2); ra(a) = 1:N; rb(b) = 1:N;
[~, ib] = min(ra + rb);
Xtr = dn(C{ib}, cfg); Utr = nrm(Xtr, cfg);
if size(unique(Xtr,'rows'),1) < n, error('dial rounding made duplicates'); end
pool = rand(20000,3); Uall = Utr; Ute = zeros(0,3);
for j = 1:cfg.nTest
    [~,k] = max(min(dist2(pool, Uall), [], 2));
    Ute(end+1,:) = pool(k,:); Uall(end+1,:) = pool(k,:); pool(k,:) = []; %#ok<AGROW>
end
Xte = dn(Ute, cfg); Ute = nrm(Xte, cfg);
X = [Xtr; Xtr(1,:); Xte; Xtr(13,:)];
phs = [ones(n,1); 3; 2*ones(cfg.nTest,1); 3];
plan = [(1:size(X,1))' phs X];
fid = fopen(cfg.planFile,'w'); fprintf(fid,'run,phase,x1,x2,x3\n');
fprintf(fid,'%d,%d,%.6g,%.6g,%.6g\n', plan'); fclose(fid);
if verbose
    fprintf('\n=== Sampling plan: %d LHS + %d hold-out + 2 repeats ===\n', n, cfg.nTest);
    fprintf('phi_q   %.3f  (median random LHS %.3f)\n', ph(ib), median(ph));
    fprintf('L2-star %.4f (median random LHS %.4f)\n', l2(ib), median(l2));
    fprintf('D*      %.4f train | %.4f train+test\n', dstar(Utr), dstar([Utr;Ute]));
    fprintf('min distance %.3f (normalised)\n', sqrt(min(nonzeros(triu(dist2(Utr,Utr),1)))));
    fprintf('Run sheet: %s\n', cfg.planFile);
end
end

function U = lhs(n)
U = zeros(n,3); for k = 1:3, U(:,k) = (randperm(n)' - 0.5)/n; end
end

function p = phi_q(U)
d = sqrt(nonzeros(triu(dist2(U,U),1))); d0 = min(d);
p = sum((d/d0).^(-50))^(1/50) / d0;
end

function D = l2star(U)
[n,d] = size(U); P = ones(n);
for k = 1:d, P = P .* (1 - max(U(:,k), U(:,k)')); end
D = sqrt(max(3^-d - 2^(1-d)/n*sum(prod(1-U.^2,2)) + sum(P(:))/n^2, 0));
end

function D = dstar(U)
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
function F = fit_best(L, cfg, verbose)
if nargin < 3, verbose = false; end
F = compare(L, cfg, [], verbose);
end

function F = compare(L, cfg, Lte, verbose)
% average repeats, normalise, compare all models by fast LOOCV
[X,~,j] = unique(L(:,3:5), 'rows'); y = accumarray(j, L(:,6), [], @mean);
U = nrm(X, cfg);
lam = cfg.lambdaDet; if is_noisy(L), lam = cfg.lambdaNoisy; end
R = {};
for p = cfg.polyOrders
    s = struct('type','poly','order',p,'kernel','','width',NaN,'lambda',0);
    [e, yl] = loo(U, y, s);
    R(end+1,:) = {sprintf('poly order %d', p), s, e, yl}; %#ok<AGROW>
end
for k = 1:numel(cfg.kernels)
    ker = cfg.kernels{k};
    if any(strcmp(ker, {'cubic','thinplate'})), W = 1; else, W = cfg.widths; end
    bestE = inf;
    for w = W
        for lm = lam
            s = struct('type','rbf','order',NaN,'kernel',ker,'width',w,'lambda',lm);
            [e, yl] = loo(U, y, s);
            if e < bestE, bestE = e; bs = s; byl = yl; end
        end
    end
    if isinf(bestE), continue; end
    if numel(W) > 1, nm = sprintf('RBF %s s=%.3f', ker, bs.width); else, nm = ['RBF ' ker]; end
    R(end+1,:) = {nm, bs, bestE, byl}; %#ok<AGROW>
end
E = cell2mat(R(:,3)); [~,o] = sort(E); R = R(o,:);
sd = std(y);
F.name = R{1,1}; F.spec = R{1,2}; F.nrmse = R{1,3}; F.q2 = 1 - R{1,3}^2*numel(y)/(numel(y)-1);
F.rmse = F.nrmse*sd; F.y = y; F.yloo = R{1,4};
F.m = fit(U, y, F.spec);
if verbose
    hasT = ~isempty(Lte);
    fprintf('%-30s %9s %8s', 'model', 'LOO NRMSE', 'LOO Q2');
    if hasT, fprintf(' %10s', 'test NRMSE'); end; fprintf('\n');
    for i = 1:size(R,1)
        fprintf('%-30s %9.3f %8.3f', R{i,1}, R{i,3}, 1 - R{i,3}^2*numel(y)/(numel(y)-1));
        if hasT
            mt = fit(U, y, R{i,2});
            fprintf(' %10.3f', sqrt(mean((Lte(:,6) - predict(mt, nrm(Lte(:,3:5),cfg))).^2))/sd);
        end
        fprintf('\n');
    end
    fprintf('Selected: %s\n', F.name);
end
end

function m = fit(U, y, s)
m.s = s; m.mu = mean(y); m.sd = std(y); if m.sd == 0, m.sd = 1; end
z = (y - m.mu)/m.sd;
if strcmp(s.type, 'poly')
    m.E = poly_exp(s.order); m.b = poly_basis(U, m.E) \ z;
else
    [K, n] = rbf_system(U, s);
    a = K \ [z; zeros(4,1)];
    m.C = U; m.w = a(1:n); m.c = a(n+1:end);
end
end

function y = predict(m, U)
if strcmp(m.s.type, 'poly')
    z = poly_basis(U, m.E) * m.b;
else
    z = kern(sqrt(dist2(U, m.C)), m.s) * m.w + [ones(size(U,1),1) U] * m.c;
end
y = m.mu + m.sd * z;
end

function [e, yloo] = loo(U, y, s)
% Exact closed-form LOOCV. Poly: e_i = r_i/(1-h_ii). RBF (Rippa 1999):
% e_i = a_i / (K^-1)_ii. Both are exact for these linear smoothers and
% are invariant to the z-scoring because a constant term is in the model.
n = numel(y); z = (y - mean(y))/std(y);
if strcmp(s.type, 'poly')
    P = poly_basis(U, poly_exp(s.order));
    if size(P,2) >= n - 1, e = inf; yloo = nan(n,1); return; end
    A = P / (P'*P); h = sum(A .* P, 2); r = z - P*(P\z);
    if any(h > 1 - 1e-8), e = inf; yloo = nan(n,1); return; end
    ez = r ./ (1 - h);
else
    [K, ~] = rbf_system(U, s);
    if rcond(K) < 1e-15, e = inf; yloo = nan(n,1); return; end
    Ki = inv(K); a = Ki * [z; zeros(4,1)];
    ez = a(1:n) ./ diag(Ki(1:n,1:n));
end
e = sqrt(mean(ez.^2));           % = NRMSE in original units
yloo = y - ez*std(y);
end

function [K, n] = rbf_system(U, s)
n = size(U,1); P = [ones(n,1) U];
K = [kern(sqrt(dist2(U,U)), s) + s.lambda*eye(n), P; P', zeros(4)];
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

function E = poly_exp(p)
[a,b,c] = ndgrid(0:p); E = [a(:) b(:) c(:)]; E = E(sum(E,2) <= p, :);
[~,i] = sort(sum(E,2)); E = E(i,:);
end

function P = poly_basis(U, E)
V = U - 0.5; P = ones(size(U,1), size(E,1));
for j = 1:size(E,1), for k = 1:3, P(:,j) = P(:,j) .* V(:,k).^E(j,k); end, end
end

function [uB, yB, peaks] = maximise(m)
% grid -> up to 6 separated seeds -> fminsearch; peaks = [u y] sorted
g = linspace(0,1,21); [a,b,c] = ndgrid(g,g,g); G = [a(:) b(:) c(:)];
[~,o] = sort(predict(m, G), 'descend'); S = G(o(1),:);
for i = o(2:end)'
    if min(dist(G(i,:), S)) > 0.25, S(end+1,:) = G(i,:); end %#ok<AGROW>
    if size(S,1) >= 6, break; end
end
f = @(u) -predict(m, min(max(u,0),1)) + 1e3*sum(max(0,u-1).^2 + max(0,-u).^2);
opt = optimset('TolX',1e-6,'TolFun',1e-8,'MaxFunEvals',2000,'Display','off');
peaks = zeros(0,4);
for i = 1:size(S,1)
    u = min(max(fminsearch(f, S(i,:), opt), 0), 1);
    if isempty(peaks) || min(dist(u, peaks(:,1:3))) > 0.15
        peaks(end+1,:) = [u predict(m,u)]; %#ok<AGROW>
    end
end
[~,o] = sort(peaks(:,4), 'descend'); peaks = peaks(o,:);
uB = peaks(1,1:3); yB = peaks(1,4);
end

% =====================================================================
%  REHEARSAL
% =====================================================================
function out = rehearse(cfg)
if exist(cfg.logFile,'file'), delete(cfg.logFile); end
make_plan(cfg, true);
box = @(x) fake_box(x);
t0 = tic; run_loop(cfg, box, true); t = toc(t0);
L = load_log(cfg);
obj = @(u) -fake_box(dn(min(max(u,0),1), cfg, false)); fb = inf;
for s = 1:20
    [u,fv] = fminsearch(obj, rand(1,3), optimset('Display','off'));
    if fv < fb, fb = fv; ut = min(max(u,0),1); end
end
ib = best_obs(L);
fprintf('\nTRUE optimum : [%.4g %.4g %.4g]  MCSI %.5g\n', dn(ut,cfg,false), -fb);
fprintf('Best observed: [%.4g %.4g %.4g]  MCSI %.5g  (%.3f%% below)\n', L(ib,3:5), L(ib,6), 100*(-fb-L(ib,6))/(-fb));
fprintf('Phase counts (LHS hold rep exploit explore local confirm): %s\n', mat2str(accumarray(L(:,2), 1, [7 1])'));
fprintf('Compute time %.1f s\n', t);
out = L;
end

function y = fake_box(x)
a = x(:,1)/200; b = x(:,2)/60; c = x(:,3);
y = 10*((1-exp(-3*a)).*(0.6+0.4*(1-exp(-5*b))).*(1+0.5*exp(-((c-0.35-0.3*b)/0.22).^2)) ...
     - 0.45*a.^2 - 0.25*b.^2 - 0.1*a.*b);
end

function worst = selftest(cfg)
% closed-form LOOCV must equal brute-force refit-and-predict
rand('seed', 1); U = rand(30,3); y = fake_box(dn(U,cfg,false)); %#ok<RAND>
specs = {struct('type','poly','order',2,'kernel','','width',NaN,'lambda',0), ...
         struct('type','rbf','order',NaN,'kernel','gaussian','width',0.4,'lambda',1e-10), ...
         struct('type','rbf','order',NaN,'kernel','multiquadric','width',0.4,'lambda',1e-3), ...
         struct('type','rbf','order',NaN,'kernel','thinplate','width',1,'lambda',1e-10)};
worst = 0;
for k = 1:numel(specs)
    [~, yf] = loo(U, y, specs{k}); yb = zeros(30,1);
    for i = 1:30
        j = true(30,1); j(i) = false;
        yb(i) = predict(fit(U(j,:), y(j), specs{k}), U(i,:));
    end
    d = max(abs(yf - yb)) / range_(y); worst = max(worst, d);
    fprintf('%-14s %-14s max |fast - brute| = %.2e of range\n', specs{k}.type, specs{k}.kernel, d);
end
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
function D = dist2(A, B), D = max(sum(A.^2,2) + sum(B.^2,2)' - 2*(A*B'), 0); end
function d = dist(a, B), d = sqrt(sum((B - a).^2, 2)); end
function tf = is_sampled(x, L), tf = any(all(abs(L(:,3:5) - x) < 1e-9, 2)); end
function i = best_obs(L), [~, i] = max(L(:,6)); end
function r = range_(v), r = max(v) - min(v); if r == 0, r = 1; end, end
function q = quantile_(v, p), v = sort(v(:)); q = v(max(1, round(p*numel(v)))); end
function tf = is_noisy(L)
tf = false;
for r = find(L(:,2)==3)'
    i = find(all(abs(L(1:r-1,3:5) - L(r,3:5)) < 1e-9, 2), 1);
    if ~isempty(i) && abs(L(r,6) - L(i,6)) > 1e-9*max(1, abs(L(r,6))), tf = true; end
end
end
