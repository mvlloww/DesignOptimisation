function plan = make_sampling_plan(cfg, makePlots)
% MAKE_SAMPLING_PLAN  Run BEFORE the session. Builds the fixed part of the
% run sheet: maximin LHS training set (screened on phi_q AND D*), a
% space-filling hold-out set, and repeat calls.
%
% Phase codes in the run sheet:
%   1 = training (LHS), 2 = hold-out test, 3 = repeat, 4 = infill, 5 = confirm
if nargin < 1, cfg = bb_config(); end
if nargin < 2, makePlots = true; end
rand('seed', cfg.seed); randn('seed', cfg.seed);   %#ok<RAND> works in MATLAB+Octave
try rng(cfg.seed); catch, end

n = cfg.nTrain; d = 3; N = cfg.nCandidates;

% ---- 1. Screen N random LHS on both criteria --------------------------
phis = zeros(N,1); l2s = zeros(N,1); C = cell(N,1);
for i = 1:N
    C{i} = lhs_random(n, d);
    phis(i) = phi_q(C{i});
    l2s(i)  = l2star(C{i});
end
% rank-sum: good spread (maximin) AND good uniformity (low discrepancy)
[~, r1] = sort(phis); [~, r2] = sort(l2s);
rk1(r1) = 1:N; rk2(r2) = 1:N;
[~, best] = min(rk1 + rk2);
Utrain = C{best};

% ---- 2. Snap to dial resolution, re-check -------------------------------
Xtrain = denormalise(Utrain, cfg, true);
Utrain = normalise(Xtrain, cfg);
if size(unique(Xtrain, 'rows'), 1) < n
    error('Dial rounding created duplicate points: reduce cfg.res or re-seed');
end

% ---- 3. Hold-out set: greedy maximin augmentation -----------------------
% Test points sit in the gaps of the training set, so they genuinely test
% interpolation rather than sitting on top of training data.
pool = rand(20000, d);
Utest = zeros(0, d); Uall = Utrain;
for j = 1:cfg.nTest
    dmin = min(sqdist(pool, Uall), [], 2);
    [~, k] = max(dmin);
    Utest(end+1,:) = pool(k,:); %#ok<AGROW>
    Uall = [Uall; pool(k,:)];   %#ok<AGROW>
    pool(k,:) = [];
end
Xtest = denormalise(Utest, cfg, true);
Utest = normalise(Xtest, cfg);

% ---- 4. Report plan quality ---------------------------------------------
Dtrain = dstar_exact(Utrain);
Dall   = dstar_exact([Utrain; Utest]);
Drand  = median(l2s);
fprintf('\n=== Sampling plan (n_train = %d, n_test = %d) ===\n', n, cfg.nTest);
fprintf('Chosen candidate %d of %d\n', best, N);
fprintf('phi_q        : %.3f   (median random LHS %.3f)\n', phis(best), median(phis));
fprintf('L2-star      : %.4f  (median random LHS %.4f, %.0f%% lower)\n', ...
        l2star(Utrain), Drand, 100*(1 - l2star(Utrain)/Drand));
fprintf('D* exact     : train %.4f | train+test %.4f\n', Dtrain, Dall);
fprintf('min distance : %.3f (normalised)\n', sqrt(min(nonzeros(triu(sqdist(Utrain,Utrain),1)))));

% ---- 5. Assemble run sheet -----------------------------------------------
% Repeats: point 1 straight after training, point 13 at the end of the
% fixed block. Two repeats far apart in time also catch drift.
rep = [1, 13];
X = [Xtrain; Xtrain(rep(1),:); Xtest; Xtrain(rep(2),:)];
phase = [ones(n,1); 3; 2*ones(cfg.nTest,1); 3];
run = (1:size(X,1))';
plan = [run phase X];

fid = fopen(cfg.planFile, 'w');
fprintf(fid, 'run,phase,x1,x2,x3\n');
fprintf(fid, '%d,%d,%.4g,%.4g,%.4g\n', plan');
fclose(fid);
fprintf('Run sheet written to %s (%d fixed calls, %d left for infill/confirm)\n\n', ...
        cfg.planFile, size(X,1), cfg.budget - size(X,1));

if makePlots
    f = figure('Name','Sampling plan','Color','w');
    pairs = [1 2; 1 3; 2 3];
    for p = 1:3
        subplot(1,3,p); hold on; box on; grid on;
        a = pairs(p,1); b = pairs(p,2);
        plot(Xtrain(:,a), Xtrain(:,b), 'ko', 'MarkerFaceColor', [0.2 0.4 0.8]);
        plot(Xtest(:,a),  Xtest(:,b),  'rs', 'MarkerFaceColor', [0.9 0.3 0.2]);
        xlabel(cfg.names{a}); ylabel(cfg.names{b});
        xlim([cfg.lb(a) cfg.ub(a)]); ylim([cfg.lb(b) cfg.ub(b)]);
    end
    legend({'LHS training','hold-out'}, 'Location','best');
    try print(f, 'fig_sampling_plan.png', '-dpng', '-r200'); catch, end
end
end

function D = sqdist(A, B)
D = max(sum(A.^2,2) + sum(B.^2,2)' - 2*(A*B'), 0);
end
