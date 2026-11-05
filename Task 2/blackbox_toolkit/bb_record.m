function bb_record(y)
% BB_RECORD  Log the MCSI reading for the point printed by bb_next.
% Writes to CSV immediately, so a MATLAB crash loses nothing.
cfg = bb_config();
if ~exist('bb_pending.mat', 'file'), error('Run bb_next first.'); end
S = load('bb_pending.mat');
L = bb_load(cfg);
run = size(L,1) + 1;
row = [run S.phase S.x y S.ypred];
newFile = ~exist(cfg.logFile, 'file');
fid = fopen(cfg.logFile, 'a');
if newFile, fprintf(fid, 'run,phase,x1,x2,x3,y,ypred\n'); end
fprintf(fid, '%d,%d,%.6g,%.6g,%.6g,%.10g,%.10g\n', row);
fclose(fid);
delete('bb_pending.mat');
fprintf('Logged run %d: MCSI = %.6g', run, y);
if ~isnan(S.ypred)
    yr = max([L(:,6); y]) - min([L(:,6); y]);
    fprintf('  (predicted %.4g, error %+.1f%% of observed range)', S.ypred, 100*(y - S.ypred)/yr);
end
fprintf('\n');
if S.phase == 3                                     % repeat check
    same = all(abs(L(:,3:5) - S.x) < 1e-9, 2);
    y0 = L(find(same,1), 6);
    fprintf('Repeat of run %d: %.6g vs %.6g -> diff %.3g\n', L(find(same,1),1), y0, y, y - y0);
    if abs(y - y0) > 1e-9 * max(1, abs(y0))
        fprintf('!! Box is NOT deterministic. Set cfg.rbfLambda = [1e-6 1e-3 1e-2 1e-1]\n');
        fprintf('   (regression not interpolation) and treat infill decisions with noise in mind.\n');
    end
end
end
