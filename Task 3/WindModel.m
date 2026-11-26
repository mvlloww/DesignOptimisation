function W = WindModel(Vmean, binWidth, Vin, Vout, hours)
% WINDMODEL  Rayleigh wind-speed distribution split into bins.
%   Only the winter MEAN wind is known (Met Office), so the spread is assumed
%   Rayleigh (Weibull with k = 2), the usual choice when only the mean is known:
%       F(V) = 1 - exp(-pi/4 (V/Vmean)^2)
%   Bins run from cut-in Vin to cut-out Vout. Each bin's probability is
%   exact (from F); power is evaluated at the bin midpoint.
if nargin < 2, binWidth = 1;  end
if nargin < 3, Vin = 3;       end
if nargin < 4, Vout = 20;     end
if nargin < 5, hours = 2160;  end              % Dec + Jan + Feb = 90 days
F = @(v) 1 - exp(-pi/4*(v/Vmean).^2);
edges  = Vin:binWidth:Vout;
W.Vmean = Vmean;  W.edges = edges;  W.hours = hours;
W.Vmid = edges(1:end-1) + binWidth/2;
W.prob = diff(F(edges));
W.pOperating = F(Vout) - F(Vin);               % fraction of time the turbine runs
end
