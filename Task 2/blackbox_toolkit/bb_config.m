function cfg = bb_config()
% BB_CONFIG  Single source of truth for the black-box session.
% Edit dial resolution once you see the actual box.

cfg.lb    = [0   0  0];      % x1 advert (£k), x2 targeting (£k), x3 skew (-)
cfg.ub    = [200 60 1];
cfg.res   = [0.01 0.01 0.01];  % smallest dial step in real units (CHECK ON THE DAY)
cfg.names = {'x1 advert (£k)', 'x2 targeting (£k)', 'x3 skew (-)'};

% Call budget (total must be <= 50)
cfg.nTrain   = 25;   % initial maximin LHS
cfg.nTest    = 5;    % independent hold-out, space-filling augment of training
cfg.nRepeat  = 2;    % repeat calls to check determinism / noise
cfg.nInfill  = 15;   % adaptive infill (exploit / explore alternating)
cfg.nConfirm = 3;    % final confirmation + local robustness
cfg.budget   = 50;

% Sampling-plan search
cfg.nCandidates = 2000;  % random LHS candidates to screen
cfg.seed        = 27;

% Surrogate candidates
cfg.polyOrders  = [1 2 3];
cfg.rbfKernels  = {'gaussian','multiquadric','invmultiquadric','cubic','thinplate'};
cfg.rbfWidths   = logspace(-1.3, 0.5, 20);   % shape parameter grid, normalised space
cfg.rbfLambda   = 1e-10;                      % nugget; raise (e.g. [1e-6 1e-3 1e-2]) if repeats show noise

cfg.logFile  = 'bb_log.csv';
cfg.planFile = 'bb_plan.csv';
end
