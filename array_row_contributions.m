%% Contribution of the upper and lower array rows to RF shimming
% Reproduces the CP and SAR-constrained comparisons for all channels and
% for each eight-channel row separately. Flip angles are in degrees and
% peak 10 g local SAR is in W/kg.

%% Load dataset
paths = setup_paths();

ld_meas = load(fullfile(paths.data, "meas_DXHB-3KQV.mat"), ...
    "B_plus", "splitix", "mask", "cmask", "bmask");
B_plus = ld_meas.B_plus;
splitix = ld_meas.splitix;
mask = ld_meas.mask;
cmask = ld_meas.cmask;
bmask = ld_meas.bmask;

nT_per_V = B_plus.*1e9;

ld_VOP = load(fullfile(paths.data, "VOPs.mat"));
VOPs = double(ld_VOP.VOPs);

%% Run all anatomical ROI and channel-row combinations
s_masks = {'all', 'brain', 'sc'};
s_chs   = {'all', 'upper', 'lower'};

nMasks = numel(s_masks);
nChs   = numel(s_chs);

out_tmp = cell(nMasks * nChs, 1);

for ii = 1:numel(out_tmp)
    [jj, kk] = ind2sub([nMasks, nChs], ii);

    out_tmp{ii} = run_shim( ...
        nT_per_V, splitix, VOPs, mask, ...
        s_masks{jj}, s_chs{kk});

    fprintf('done: mask=%s, channels=%s\n', ...
        s_masks{jj}, s_chs{kk});
end

out = reshape(out_tmp, [nMasks, nChs]);

%% Plot equal-SAR maps and VOP-constrained L-curves
% The transform preserves the sagittal orientation used for plotting.
imgtr = @(x) flipud(rot90(squeeze(x)));

f = figure;
p = panel();
p.pack(3, 5);

CAX = [0,10];

% CP reference maps.
for jj = 1:3
    p(jj, 1).select();
    MAP = out{jj,1}.FA_CP_all_2D;
    if jj == 1
        MASK = imgtr(cmask);
    elseif jj == 2
        MASK = imgtr(cmask);
        MASK(1:splitix,:) = NaN;
    else
        MASK = imgtr(cmask);
        MASK(splitix+1:end,:) = NaN;
    end
    imagescc(imgtr(abs(MAP)).*bmask, 0, 1, imgtr(MASK));
    clim(CAX);
    colormap(gca,'turbo'); colorbar off;
    title(sprintf('%.2f %%', out{jj,1}.nrmse_CP_2D.*100));
end

% Select the constrained solution closest to the all-channel CP SAR.
for jj = 1:3
    for kk = 1:3
        refSAR = out{jj,1}.VOPSAR_CP_2d;
        [~, ixpower] = min(abs(refSAR - out{jj,kk}.VOPSAR_constr_2d));

        p(jj,kk+1).select();
        MAP = out{jj,kk}.FA_constr_all_2D(:,:,:,ixpower);
        if jj == 1
            MASK = imgtr(cmask);
        elseif jj == 2
            MASK = imgtr(cmask);
            MASK(1:splitix,:) = NaN;
        else
            MASK = imgtr(cmask);
            MASK(splitix+1:end,:) = NaN;
        end
        imagescc(imgtr(abs(MAP)).*bmask, 0, 1, imgtr(MASK));
        clim(CAX);
        colormap(gca,'turbo'); colorbar off;
        title(sprintf('%.2f %%', out{jj,kk}.nrmse_constr_2D(ixpower).*100));
    end
end

% L-curves.
LW = 2;
lcs = lines;
for jj = 1:3
    for kk = 1:3
        p(jj, 5).select();
        hold on;

        if kk == 1
            plot(out{jj,1}.nrmse_CP_2D.*100, out{jj,1}.VOPSAR_CP_2d, ...
                'xk', 'LineWidth', 1.5*LW, 'MarkerSize', 10);
            plot([0, 60], out{jj,1}.VOPSAR_CP_2d.*[1, 1], '--', ...
                'Color', '#AAAAAA', 'LineWidth', LW, 'HandleVisibility', 'off');
        end
        xx = out{jj,kk}.nrmse_constr_2D.*100;
        yy = out{jj,kk}.VOPSAR_constr_2d;
        [yy, sortix] = sort(yy);
        xx = xx(sortix);
        plot(xx, yy, 'LineWidth',LW, 'Color', lcs(kk,:));

        ylim([0, 20]);
        xlim([0, 60]);
        box on;
        xlabel('NRMSE [%]');
        ylabel('pSAR_{10g} [W/kg]');
        box on;
    end
end
p(1, 5).select();
legend('CP mode', 'all ch.', 'upper ch. only', 'lower ch. only');

p(1, 1).select(); set_img_ylabel('brain + C-spine');
p(2, 1).select(); set_img_ylabel('brain');
p(3, 1).select(); set_img_ylabel('C-spine');

all_images_black_background(gcf);
p.ch.margin = 10;
p.fontsize = 14;
font_to_allAxes;

%% Local RF-shimming function
function out = run_shim(nT_per_V, splitix, VOPs, mask, maskmode, channelmode)
%RUN_SHIM Compute CP and SAR-constrained shims for one ROI/channel subset.
%   OUT = RUN_SHIM(NT_PER_V, SPLITIX, VOPS, MASK, MASKMODE, CHANNELMODE)
%   builds a small-tip-angle system matrix from B1+ maps NT_PER_V (nT/V).
%   SPLITIX separates brain from C-spine, VOPS contains local-SAR matrices,
%   MASK defines the anatomical target, MASKMODE selects all/brain/C-spine,
%   and CHANNELMODE selects all/upper/lower transmit channels. OUT contains
%   flip-angle maps (degrees), NRMSE, RF weights, average power (W), and peak
%   10 g local SAR (W/kg).

% The row split follows the acquisition's brain/C-spine boundary.
if strcmp(maskmode, 'brain')
    mask(1:splitix,:) = 0;
elseif strcmp(maskmode, 'sc')
    mask(splitix+1:end,:) = 0;
end

% Build optimization and full-field display system matrices.
B_plus = nT_per_V .* 1e-9;

stepdur = 2e-3;
TR = stepdur;

optmask = double(mask);
optmask(optmask == 0) = NaN;

stepSize = [2,2,3].*1e-3;
dist = createVoxelGrid(optmask, stepSize) .* optmask;
dist = reshape(dist(~isnan(dist)), [], 3);

% The unmasked matrix is used only to display the full sagittal field.
plotmask = ones(size(optmask));

sp = flattenize(B_plus.*plotmask);
distp = ones(size(sp,1), 3);
Ap = createAmat(sp, distp, [0, 0, 0], stepdur);

targFA = 5;


B1m = double(B_plus .* optmask);
s = flattenize(B1m);

A = createAmat(s, dist, [0, 0, 0], stepdur);

ix_upper = 1:2:16; % elements
ix_lower = 2:2:16;
if strcmp(channelmode, 'upper')
    A = A(:,ix_upper);
    Ap = Ap(:,ix_upper);
    VOPs = VOPs(ix_upper,ix_upper,:);
    p_CP = getCPmode(8, 1, 1);
elseif strcmp(channelmode, 'lower')

    A = A(:,ix_lower);
    Ap = Ap(:,ix_lower);
    VOPs = VOPs(ix_lower,ix_lower,:);
    p_CP = getCPmode(8, 1, 1);
else
    p_CP = getCPmode(16, 1, 1);
end

% Scale the CP reference to the target mean flip angle within the ROI.
FAtemp = calcFA(p_CP, A);
scalefactor_cp = targFA ./ mean(abs(FAtemp));
p_CP = p_CP .* scalefactor_cp;
[~, nrmse_CP_2D, power_CP_2D, VOPSAR_CP_2d] = ...
    eval_pulse(p_CP, A, optmask, targFA, VOPs, [], stepdur, TR);
[FA_CP_all_2D, ~] = ...
    eval_pulse(p_CP, Ap, plotmask, targFA, VOPs, [], stepdur, TR);


out = struct();
out.FA_CP_all_2D = FA_CP_all_2D;
out.nrmse_CP_2D = nrmse_CP_2D;
out.power_CP_2D = power_CP_2D;
out.VOPSAR_CP_2d = VOPSAR_CP_2d;
out.p_CP = p_CP;

%% Trace the VOP-constrained L-curve

minSAR = 0.1;
maxSAR = 1000;
NSAR = 256;
dowarmstart = 1;
maxiter = 500;

p0 = p_CP;

[~, ~, VOPSAR0] = calcLcurve(double(A), targFA, p0, VOPs, maxSAR, ...
    optmask, 0, 'interior-point', stepdur, TR, maxiter);

% Sweep both directions because warm starts can follow different local minima.
SARLIMs_L = (linspace(minSAR, VOPSAR0, NSAR));

[pp, err, sar] = calcLcurve(double(A), targFA, p0, VOPs, SARLIMs_L, ...
    optmask, dowarmstart, 'interior-point', stepdur, TR, maxiter);

SARLIMs_L_f = flip(linspace(minSAR, VOPSAR0, NSAR));

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
out.p_constr = p_pareto;

end
