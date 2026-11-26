function [E, d] = WinterEnergy(B, c, twist, ctrl, W, P, L, stall)
% WINTERENERGY  Winter (Dec-Feb) energy of a BUILT blade, in kWh.
%
%   E = sum over wind bins of  power(V_bin) x probability(bin) x hours
%
%   B, c, twist : blade (c in m and twist in deg, one value per element)
%   ctrl.mode   : 'variable' -> lambda = min(ctrl.lamOpt, tipMax/V)
%                 (track the best lambda until the tip-speed limit, then
%                 hold omega constant)
%                 'fixed'    -> omega fixed by ctrl.tip (tip speed, m/s),
%                 so lambda = tip/V changes with the wind
%   W           : wind model from WindModel (bin edges, probabilities, hours)
%   P, L        : settings struct and S822 lookup (as in OptimisePeakCp)
%   stall       : 'viterna' (default) or 'clamp'
% Power below zero (rotor would need driving at low wind) is set to zero:
% the generator just stops drawing power.

if nargin < 8, stall = 'viterna'; end
o.lookup = L; o.tipLoss = true; o.alphaMode = 'twist'; o.relax = 0.3; o.stall = stall;
nb = numel(W.Vmid);
[Pw, lam, Cp, conv] = deal(zeros(1, nb));
for k = 1:nb
    V = W.Vmid(k);
    if strcmp(ctrl.mode, 'variable'), lam(k) = min(ctrl.lamOpt, P.tipMax/V);
    else,                             lam(k) = ctrl.tip/V;
    end
    s = BEMT(B, [P.Rhub P.R], lam(k), c, twist, numel(c), V, P.nu, P.tol, 5000, o);
    Pw(k) = max(s.P, 0); Cp(k) = s.Cp; conv(k) = mean(s.converged);
end
Eb = Pw.*W.prob*W.hours/1000;                 % kWh from each bin
E  = sum(Eb);
d  = struct('V', W.Vmid, 'P', Pw, 'lambda', lam, 'Cp', Cp, 'Ebin', Eb, ...
            'convFrac', conv, 'allConverged', all(conv == 1));
end
