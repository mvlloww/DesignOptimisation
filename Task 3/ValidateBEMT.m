% VALIDATEBEMT  Check that BEMT.m reproduces the supplied BEMT_Demo.m.
% With tip loss off and fixed AoA, both should give the same Cp. Also times
% both, and shows how much the Prandtl tip-loss correction changes Cp.
clear; clc;
here = fileparts(mfilename('fullpath'));
demoDir = fullfile(here, 'Wind Turbine Blade Optimisation Demo');
addpath(demoDir);
L = S822Lookup(demoDir);

% 1. Lookup: scatteredInterpolant vs griddata at random points
A  = load(fullfile(demoDir,'S822_Alpha_Values.mat')); Rv = load(fullfile(demoDir,'S822_Re_Values.mat'));
Cl = load(fullfile(demoDir,'S822_CL_Lookup.mat'));
rng(1); aq = 12*rand(200,1); rq = 2e4 + 8.3e5*rand(200,1);
g = arrayfun(@(x,y) griddata(A.S822_Alpha_Values, Rv.S822_Re_Values, Cl.S822_CL_Lookup', x, y), aq, rq);
fprintf('Lookup: max |CL(griddata) - CL(new)| = %.2e over 200 random points\n', max(abs(g - L.CL(aq, rq))));

% 2. Whole solver: BEMT_Demo vs BEMT on a few test blades
cases = { 3, [0.15 1.1], 6, [0.12 0.05], [8 5];
          2, [0.15 1.1], 8, [0.10 0.04], [7 4];
          3, [0.20 0.5], 4, [0.08 0.038], [5 2];     % the demo's own example
          4, [0.15 1.1], 5, [0.15 0.06], [9 6] };
V = 7.63; nu = 1.5e-5; N = 30; tol = 1e-8; it = 500;
fprintf('\n B  lambda   Cp(demo)   Cp(new)    diff      t_demo  t_new\n');
for k = 1:size(cases,1)
    [B, rl, lam, ch, al] = cases{k,:};
    tic; evalc('cpD = BEMT_Demo(B, rl, lam, ch, al, N, V, nu, tol, it, false);'); tD = toc;
    tic; o = BEMT(B, rl, lam, ch, al, N, V, nu, tol, it, struct('lookup', L)); tN = toc;
    fprintf('%2d  %5.1f   %.6f   %.6f   %.1e   %5.2fs  %5.3fs\n', B, lam, cpD, o.Cp, abs(cpD-o.Cp), tD, tN);
end

% 3. Effect of Prandtl tip loss (+ Buhl correction) on the first case
o0 = BEMT(3, [0.15 1.1], 6, [0.12 0.05], [8 5], N, V, nu, tol, it, struct('lookup', L));
o1 = BEMT(3, [0.15 1.1], 6, [0.12 0.05], [8 5], N, V, nu, tol, it, struct('lookup', L, 'tipLoss', true));
fprintf('\nTip loss: Cp %.4f -> %.4f (%.1f%% lower), all elements converged: %d\n', ...
    o0.Cp, o1.Cp, 100*(1 - o1.Cp/o0.Cp), all(o1.converged));
