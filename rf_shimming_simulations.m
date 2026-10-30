%% RF shimming simulations for the Duke and Ella voxel models
% Reproduces the simulated CP and SAR-constrained RF-shim comparison.
% Flip angles are in degrees and peak 10 g local SAR is in W/kg.

%% Load data

paths = setup_paths();

ld_data = load(fullfile(paths.data, "simu_data.mat"));
ld_VOP = load(fullfile(paths.data, "VOPs.mat"));
VOPs = double(ld_VOP.VOPs);

%% Select anatomical model
modelid = 1; % 1: Duke, 2: Ella
if modelid == 1
    maskpath = fullfile(paths.data, 'bcmask_Duke.mat');
else
    maskpath = fullfile(paths.data, 'bcmask_Ella.mat');
end

%% Run CP and constrained RF shims
tic

ALLOUT = {};
ALLOUT{1} = simu_shim(maskpath, ld_data.simu_data, VOPs, modelid, 0, 1); 
ALLOUT{2} = simu_shim(maskpath, ld_data.simu_data, VOPs, modelid, 0, 2);
ALLOUT{3} = simu_shim(maskpath, ld_data.simu_data, VOPs, modelid, 0, 3);
ALLOUT{4} = simu_shim(maskpath, ld_data.simu_data, VOPs, modelid, 1, 1);

toc

%% Plot L-curves and representative flip-angle maps

PALLOUT = {ALLOUT{[1, 2, 3, 4]}};

f = figure;
p = panel();
p.pack(length(PALLOUT), {.22, .22, .22, []});

for kk = 1:length(PALLOUT)
    out = PALLOUT{kk};
    p(kk, 4).select();
    hold on;
    plot(out.nrmse_CP_2D(1).*100, out.VOPSAR_CP_2d(1), 'kx', ...
        'LineWidth', 1.5*1.5, 'MarkerSize', 8);
    [xsc, ysc] = convhull_sorted(out.nrmse_constr_2D.*100, out.VOPSAR_constr_2d);
    plot(xsc, ysc, '.-', 'LineWidth', 1.5);
    xlim([0, 40]);
    ylim([0, 50]);
    set(gca, 'yscale', 'log');
    xl = xlim;

    plot(xl, out.VOPSAR_CP_2d(1).*[1, 1], '--', 'Color', '#AAAAAA');
    box on;
    xlabel('NRMSE [%]');
    ylabel('pSAR_{10g} [W/kg]');
    legend('CP', 'RF shim');

    sz = size(out.FA_CP_all_2D);
    ix = round(sz./2);

    % Compare the CP map with equal-SAR and minimum-NRMSE RF shims.
    [~, ixc1] = min(abs(out.VOPSAR_CP_2d(1) - out.VOPSAR_constr_2d));
    [~, ixc2] = min(out.nrmse_constr_2D);

    SHOWMAPS = {out.FA_CP_all_2D, out.FA_constr_all_2D(:,:,:,ixc1), out.FA_constr_all_2D(:,:,:,ixc2)};
    TIT = {{'', sprintf('%.2f %%', out.nrmse_CP_2D(1).*100)}, ...
        {'', sprintf('%.2f %%', out.nrmse_constr_2D(ixc1).*100)}, ...
        {'', sprintf('%.2f %%', out.nrmse_constr_2D(ixc2).*100)}};

    % Retain the publication's sagittal orientation and inferior crop.
    ctrafo = @(x) flipud(squeeze(x(ix(1),:,45:end,1)));

    for jj = 1:length(SHOWMAPS)
        p(kk, jj).select();
        imagescc(abs(ctrafo(SHOWMAPS{jj})), 0, 1, ctrafo(out.optmask));
        caxis([0, 10]);
        colormap(turbo);
        cb = colorbar;
        title(cb, 'FA [deg]');
        cb.Location = 'westoutside';
        if jj ~= 1
            colorbar off;
        end
        title(TIT{jj});
    end
end
p.margintop = 15;
p.marginleft = 20;
default_plot_settings;
f.OuterPosition = [163         107        1624        1223];
p.fontsize = 14;
font_to_allAxes;

%% Local RF-shimming function
function out = simu_shim(maskpath, SARdat, VOPs, modelid, do2d, maskid)
%SIMU_SHIM Compute CP and SAR-constrained RF shims for one voxel model.
%   OUT = SIMU_SHIM(MASKPATH, SARDAT, VOPS, MODELID, DO2D, MASKID) builds
%   small-tip-angle system matrices and evaluates the CP mode and the
%   Pareto-optimal magnitude-least-squares RF shims. MASKPATH contains the
%   brain/C-spine mask; SARDAT contains B1+ maps (T/V), tissue masks, and
%   voxel sizes; VOPS contains local-SAR matrices; MODELID selects the
%   anatomical model; DO2D restricts optimization to one sagittal plane;
%   and MASKID selects brain plus cord (1), brain (2), or C-spine (3).
%   OUT contains flip-angle maps (degrees), NRMSE, average power (W), masks,
%   and peak 10 g local SAR (W/kg).

% Index 126 marks the medulla-to-C-spine transition in the voxel model.
splitix = 126;

% Brain/C-spine target mask used by the optimization.
ld = load(maskpath, "bcmask");
bcmask = ld.bcmask;

if maskid == 2
    bcmask(:,:,1:splitix) = 0;
elseif maskid == 3
    bcmask(:,:,splitix+1:end) = 0;
end

% Use the filled tissue extent only to display maps outside the target ROI.
mask = SARdat{modelid}.tissueMask;
plotmask = double(imfill_multislice(mask));
plotmask(plotmask == 0) = NaN;

B_plus = SARdat{modelid}.B_plus;

%% Build small-tip-angle system matrices

stepdur = 2e-3;
TR = stepdur; % 100 % DC

optmask = double(bcmask);
optmask(optmask == 0) = NaN;
if do2d
    optmask = double(bcmask);
    optmask(optmask == 0) = NaN;
    ix2d = 61;
    optmask(1:ix2d-1,:,:) = NaN;
    optmask(ix2d+1:end,:,:) = NaN;
end

stepSize = SARdat{modelid}.stepSize.*1e-3;
dist = createVoxelGrid(optmask, stepSize) .* optmask;
dist = reshape(dist(~isnan(dist)), [], 3);

distp = createVoxelGrid(plotmask, stepSize) .* plotmask;
distp = reshape(distp(~isnan(distp)), [], 3);
sp = flattenize(B_plus.*plotmask);
Ap = createAmat(sp, distp, [0, 0, 0], stepdur);

targFA = 5;

B1m = double(B_plus .* optmask);
s = flattenize(B1m);

A = createAmat(s, dist, [0, 0, 0], stepdur);

%% Scale and evaluate the CP reference
row_phase_shifts = 0;
for jj = 1:length(row_phase_shifts)
    p_CP = getCPmode(16, 1, -1);
    p_CP(1:2:16) = p_CP(1:2:16).*exp(1j.*deg2rad(row_phase_shifts(jj)));
    FAtemp = calcFA(p_CP, A);
    scalefactor_cp = targFA ./ mean(abs(FAtemp));
    p_CP = p_CP .* scalefactor_cp;
    [~, nrmse_CP_2D(jj), power_CP_2D(jj), SAR_VOP_CP_2d(jj)] = ...
        eval_pulse(p_CP, A, optmask, targFA, VOPs, [], stepdur, TR);
    [FA_CP_all_2D(:,:,:,jj), ~] = ...
        eval_pulse(p_CP, Ap, plotmask, targFA, VOPs, [], stepdur, TR);
    fprintf('%d/%d\n', jj, length(row_phase_shifts))
end

out = struct();
out.FA_CP_all_2D = FA_CP_all_2D;
out.nrmse_CP_2D = nrmse_CP_2D;
out.power_CP_2D = power_CP_2D;
out.VOPSAR_CP_2d = SAR_VOP_CP_2d;

out.optmask = optmask;
out.plotmask = plotmask;

%% Trace the VOP-constrained L-curve

minSAR = 0.1; 
maxSAR = 2000;
NSAR = 256;
dowarmstart = 1;
maxiter = 500;

p0 = p_CP;

[~, ~, VOPSAR0] = calcLcurve(double(A), targFA, p0, VOPs, maxSAR, ...
    optmask, 0, 'interior-point', stepdur, TR, maxiter);

% Sweep both directions because warm starts can follow different local minima.
SARLIMs_L = exp(linspace(log(minSAR), log(VOPSAR0), NSAR));

[pp, err, sar] = calcLcurve(double(A), targFA, p0, VOPs, SARLIMs_L, ...
    optmask, dowarmstart, 'interior-point', stepdur, TR, maxiter);

SARLIMs_L_f = flip(exp(linspace(log(minSAR), log(VOPSAR0), NSAR)));

[pp_f, err_f, sar_f] = calcLcurve(double(A), targFA, p0, VOPs, SARLIMs_L_f, ...
    optmask, dowarmstart, 'interior-point', stepdur, TR, maxiter);

% Retain only non-dominated NRMSE/SAR solutions from both sweep directions.
x = cat(2, cat(1, err(:), err_f(:)), cat(1, sar(:), sar_f(:)));
[membership, ~] = find_pareto_frontier(x);

p_all = cat(2, pp, pp_f);
p_pareto = p_all(:,membership);

clear FA_constr_2D nrmse_constr_2D power_constr_2D VOPSAR_constr_2d FA_constr_all_2D
for jj = 1:size(p_pareto, 2)
    [FA_constr_2D(:,:,:,jj), nrmse_constr_2D(jj), power_constr_2D(jj), VOPSAR_constr_2d(jj)] = ...
        eval_pulse(p_pareto(:,jj), A, optmask, targFA, VOPs, [], stepdur, TR);
    [FA_constr_all_2D(:,:,:,jj), ~] = ...
        eval_pulse(p_pareto(:,jj), Ap, plotmask, targFA, VOPs, [], stepdur, TR);

    fprintf('%d/%d\n', jj, size(p_pareto, 2))
end

out.FA_constr_all_2D = FA_constr_all_2D;
out.nrmse_constr_2D = nrmse_constr_2D;
out.power_constr_2D = power_constr_2D;
out.VOPSAR_constr_2d = VOPSAR_constr_2d;

end
