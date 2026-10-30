function out = run_UP_experiment(uix, mask, B_plus, VOPs)
%RUN_UP_EXPERIMENT Train one universal RF shim and evaluate every subject.
%   OUT = RUN_UP_EXPERIMENT(UIX, MASK, B_PLUS, VOPS) trains on subject
%   indices UIX using anatomical MASK, B1+ maps B_PLUS (T/V), and local-SAR
%   matrices VOPS. OUT contains complex RF weights, flip-angle maps
%   (degrees), fractional NRMSE values, peak local SAR (W/kg), and the
%   target flip angle.

stepdur = 1e-3;
TR = 1e-3;
targFA = 5;

[nx, ny, nc, ~] = size(B_plus);
nu = length(uix);

mask_UP = reshape(mask(:,:,uix), size(mask,1), []);
B_plus_UP = reshape(permute(B_plus(:,:,:,uix), [1, 2, 4, 3]), nx, ny*nu, nc);

optmask = double(mask_UP);
optmask(optmask == 0) = NaN;

stepSize = [2,2,3].*1e-3;
dist = createVoxelGrid(optmask, stepSize) .* optmask;
dist = reshape(dist(~isnan(dist)), [], 3);

B1m_UP = double(B_plus_UP .* optmask);
s = flattenize(B1m_UP, 3);

A = createAmat(s, dist, [0, 0, 0], stepdur);

% Scale the CP reference over the combined training anatomy.
p_CP = getCPmode(16, 1, 1);
FAtemp = calcFA(p_CP, A);
scalefactor_cp = targFA ./ mean(abs(FAtemp));
p_CP = p_CP .* scalefactor_cp;
[~, ~, ~, VOPSAR_CP_2d] = eval_pulse( ...
    p_CP, A, optmask, targFA, VOPs, [], stepdur, TR);

% Constrain the universal shim to the CP mode's peak local SAR.
maxSAR = VOPSAR_CP_2d;

maxiter = 500;

p0 = p_CP;

[po0, nrmse0, VOPSAR0] = calcLcurve(double(A), targFA, p0, VOPs, maxSAR, ...
    optmask, 0, 'interior-point', stepdur, TR, maxiter);

% Rescale and evaluate the universal and CP shims subject by subject.
FA_indi = zeros(size(mask,1), size(mask,2), 1, size(mask,3));
nrmse_indi = zeros(size(mask,3), 1);
FA_CP_indi = zeros(size(mask,1), size(mask,2), 1, size(mask,3));
nrmse_CP_indi = zeros(size(mask,3), 1);
nrmse_CP_trainscale = zeros(size(mask,3), 1);
VOPSAR_CP_indi = zeros(size(mask,3), 1);
for jj = 1:size(mask,3)
    [Ai, maski] = get_Amat_individ(B_plus, mask, jj);

    FAtemp = calcFA(po0, Ai);
    scalefactor_up = targFA ./ mean(abs(FAtemp));
    po0_scaled = po0 .* scalefactor_up;

    [FA_indi(:,:,:,jj), nrmse_indi(jj), ~, ~] = ...
        eval_pulse(po0_scaled, Ai, maski, targFA, VOPs, [], stepdur, TR);
    fprintf('subject %d\n', jj)

    p_CP = getCPmode(16, 1, 1);
    FAtemp = calcFA(p_CP, Ai);
    scalefactor_cp = targFA ./ mean(abs(FAtemp));
    p_CP = p_CP .* scalefactor_cp;
    [FA_CP_indi(:,:,:,jj), nrmse_CP_indi(jj), ~, VOPSAR_CP_indi(jj)] = eval_pulse(p_CP, Ai, maski, targFA, VOPs, [], stepdur, TR);

    scalefactor_cp2 = targFA ./ mean(abs(calcFA(getCPmode(16, 1, 1), A)));
    [~, nrmse_CP_trainscale(jj), ~, ~] = eval_pulse( ...
        getCPmode(16, 1, 1).*scalefactor_cp2, ...
        Ai, maski, targFA, VOPs, [], stepdur, TR);

end

% Return both training-set and per-subject evaluation results.
out.po0 = po0;
out.nrmse0 = nrmse0;
out.VOPSAR0 = VOPSAR0;
out.FA_indi = FA_indi;
out.nrmse_indi = nrmse_indi;
out.FA_CP_indi = FA_CP_indi;
out.nrmse_CP_indi = nrmse_CP_indi;
out.VOPSAR_CP_indi = VOPSAR_CP_indi;
out.targFA = targFA;
out.nrmse_CP_trainscale = nrmse_CP_trainscale;
end
