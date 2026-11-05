function D = dstar_exact(U)
% Exact star discrepancy D* (L-infinity) for small n in 3D.
% Evaluates every anchored box [0,u) with corners on the grid of point
% coordinates (plus 1), counting open and closed boxes. O(n^4) in 3D:
% ~1e6 operations for n = 30, so it runs in well under a second.
[n, d] = size(U);
if d ~= 3, error('dstar_exact written for d = 3'); end
g1 = unique([U(:,1); 1]); g2 = unique([U(:,2); 1]); g3 = unique([U(:,3); 1]);
D = 0;
for a = g1'
    in1o = U(:,1) < a;  in1c = U(:,1) <= a;
    for b = g2'
        in12o = in1o & U(:,2) < b;  in12c = in1c & U(:,2) <= b;
        for c = g3'
            vol = a * b * c;
            no  = sum(in12o & U(:,3) <  c);
            nc  = sum(in12c & U(:,3) <= c);
            D = max([D, nc / n - vol, vol - no / n]);
        end
    end
end
end
