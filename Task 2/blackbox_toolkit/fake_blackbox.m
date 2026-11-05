function y = fake_blackbox(x, noise)
% Stand-in for rehearsal ONLY. Diminishing returns on spend, a spend
% penalty ("bang for buck"), and a skew effect that interacts with targeting.
% True optimum is interior; rehearse.m finds it numerically for comparison.
if nargin < 2, noise = 0; end
a = x(:,1)/200; b = x(:,2)/60; c = x(:,3);
reach   = (1 - exp(-3.0*a)) .* (0.6 + 0.4*(1 - exp(-5*b)));
skewFit = 1 + 0.5*exp(-((c - 0.35 - 0.3*b)/0.22).^2);
cost    = 0.45*a.^2 + 0.25*b.^2 + 0.1*a.*b;
y = 10*(reach .* skewFit - cost) + noise*randn(size(a));
end
