function [Ubest, ybest] = maximise_surrogate(m)
% Global max of a cheap surrogate: dense grid -> top seeds -> fminsearch.
% Bounds handled by a quadratic penalty (no Optimization Toolbox needed).
g = linspace(0, 1, 21);
[a,b,c] = ndgrid(g,g,g); G = [a(:) b(:) c(:)];
yg = predict_surrogate(m, G);
[~, o] = sort(yg, 'descend');
seeds = G(o(1),:);
for i = o(2:end)'                       % up to 6 well-separated seeds
    if min(sqrt(sum((seeds - G(i,:)).^2, 2))) > 0.25, seeds(end+1,:) = G(i,:); end %#ok<AGROW>
    if size(seeds,1) >= 6, break; end
end
obj = @(u) -predict_surrogate(m, min(max(u,0),1)) + 1e3*sum(max(0,u-1).^2 + max(0,-u).^2);
opts = optimset('TolX',1e-6,'TolFun',1e-8,'MaxFunEvals',2000,'Display','off');
ybest = -inf;
for s = 1:size(seeds,1)
    u = fminsearch(obj, seeds(s,:), opts);
    u = min(max(u,0),1);
    yv = predict_surrogate(m, u);
    if yv > ybest, ybest = yv; Ubest = u; end
end
end
