function L = S822Lookup(dataDir)
% S822LOOKUP  Build the S822 aerofoil lookup once, so BEMT does not have to.
%
% The supplied BEMT_Demo calls griddata(...) inside its iteration loop.
% griddata builds a Delaunay triangulation of the (alpha, Re) data every
% time it is called, which is what makes the supplied solver slow.
% Here the SAME linear interpolation on the SAME triangulation is built
% once with scatteredInterpolant and reused. ValidateBEMT.m checks that the
% result matches griddata.
%
% It also finds the angle of attack with the best lift/drag ratio at each
% Reynolds number, and fits a smooth curve alphaBest(Re) through it, so the
% optimiser can start every blade section at its best-L/D angle.

if nargin < 1
    dataDir = fullfile(fileparts(mfilename('fullpath')), 'Wind Turbine Blade Optimisation Demo');
end
A  = load(fullfile(dataDir, 'S822_Alpha_Values.mat'));  alpha = A.S822_Alpha_Values(:)';
Rd = load(fullfile(dataDir, 'S822_Re_Values.mat'));     Re    = Rd.S822_Re_Values(:)';
C  = load(fullfile(dataDir, 'S822_CL_Lookup.mat'));     CL    = C.S822_CL_Lookup;   % [nAlpha x nRe]
D  = load(fullfile(dataDir, 'S822_CD_Lookup.mat'));     CD    = D.S822_CD_Lookup;

% Same point set and ordering that griddata(alpha, Re, CL', ...) uses
[X, Y] = meshgrid(alpha, Re);                 % [nRe x nAlpha]
L.CL = scatteredInterpolant(X(:), Y(:), reshape(CL', [], 1), 'linear', 'none');
L.CD = scatteredInterpolant(X(:), Y(:), reshape(CD', [], 1), 'linear', 'none');
L.alphaRange = [min(alpha) max(alpha)];
L.ReRange    = [min(Re) max(Re)];

% Best-L/D angle at each Reynolds number (fine alpha search on the lookup)
aFine  = 0:0.01:max(alpha);
ReFine = logspace(log10(L.ReRange(1)*1.0001), log10(L.ReRange(2)*0.9999), 60);   % just inside the data edge
aBest  = zeros(size(ReFine));  LDmax = aBest;
for k = 1:numel(ReFine)
    ld = L.CL(aFine, ReFine(k)*ones(size(aFine))) ./ L.CD(aFine, ReFine(k)*ones(size(aFine)));
    ld(~isfinite(ld)) = -inf;                  % lookup is NaN exactly on the data edge
    [LDmax(k), j] = max(ld);
    aBest(k) = aFine(j);
end
% Smooth fit alphaBest = p1*ln(Re) + p2 (a step-wise curve would give the
% gradient-based optimiser zero/infinite slopes)
L.pBest     = polyfit(log(ReFine), aBest, 1);
L.alphaBest = @(re) polyval(L.pBest, log(re));
L.bestTable = struct('Re', ReFine, 'alpha', aBest, 'LD', LDmax);
end
