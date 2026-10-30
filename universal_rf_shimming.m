%% Universal RF-shimming analysis across subjects
% Reproduces convergence, leave-one-subject-out (LOSO), and flip-angle-map
% comparisons. Flip angles are in degrees and NRMSE values are fractional
% internally (converted to percent only for display).

%% Load all B1+ datasets
paths = setup_paths();

ld = load(fullfile(paths.data, "universal_shim_datasets.mat"));
B_plus = ld.B_plus;
mask = ld.mask;

% Load virtual observation point (VOP) local-SAR matrices.
ld_VOP = load(fullfile(paths.data, "VOPs.mat"));
VOPs = double(ld_VOP.VOPs);

%% Universal-shim analysis over all possible training combinations
%
% Creates:
%   A) CP-normalized convergence plot
%   B) Absolute LOSO NRMSE comparison: CP vs universal vs tailored
%   C) CP, universal and tailored flip-angle maps for representative subjects
%
% Required variables:
%   mask, B_plus, VOPs
%
% run_UP_experiment must return:
%   out_up.nrmse_indi        [Nsubjects x 1]
%   out_up.nrmse_CP_indi     [Nsubjects x 1]
%   out_up.FA_indi           [Nx x Ny x 1 x Nsubjects]
%   out_up.FA_CP_indi        [Nx x Ny x 1 x Nsubjects]

clear R

%% Settings

ntot = size(B_plus,4);

if ~isfolder(paths.results)
    mkdir(paths.results);
end
resultsFile = fullfile(paths.results, 'US_all_combinations.mat');

% Set false after the results have been generated once.
recompute = true;

% Empty: automatically select median and worst LOSO/CP subjects.
% Alternatively, specify manually, e.g.:
% mapSubjects = [3 7];
%
% Setting mapSubjects = 1:ntot shows all subjects, although the resulting
% figure will be quite tall.
mapSubjects = 1:ntot;

mapCLim = [2.5 7.5]; % Flip-angle display range [degrees].

%% Run all training combinations

if recompute || ~isfile(resultsFile)

    R = struct;

    R.ntot = ntot;
    R.ntrainList = 1:ntot-1;

    % For every training-set size:
    %   rows    = training combinations
    %   columns = evaluated subjects
    R.trainsets = cell(ntot-1,1);
    R.nrmse_US  = cell(ntot-1,1);

    % Subject-specific reference results
    R.nrmse_CP   = [];
    R.nrmse_LOSO = nan(ntot,1);
    R.nrmse_TP   = nan(ntot,1);

    R.FA_CP      = [];
    R.FA_LOSO    = [];
    R.FA_TP      = [];

    for ntrain = 1:ntot-1

        trainsets = nchoosek(1:ntot,ntrain);
        nsets = size(trainsets,1);

        R.trainsets{ntrain} = trainsets;

        % One row per training combination, one column per subject
        nrmse_US_tmp = nan(nsets,ntot);
        nrmse_CP_tmp = nan(nsets,ntot);

        % Only populated for ntrain = 1 or ntot-1
        FA_keep_tmp = cell(nsets,1);

        % Only one copy is needed
        FA_CP_tmp = cell(nsets,1);

        keepMaps = ntrain == 1 || ntrain == ntot-1;

        parfor jj = 1:nsets

            trainix = trainsets(jj,:);

            out_up = run_UP_experiment( ...
                trainix, mask, B_plus, VOPs);

            nrmse_US_tmp(jj,:) = out_up.nrmse_indi(:).';
            nrmse_CP_tmp(jj,:) = out_up.nrmse_CP_indi(:).';

            if keepMaps

                if ntrain == 1
                    % Tailored case: evaluate the training subject itself
                    subjectToKeep = trainix;
                else
                    % LOSO case: evaluate the excluded subject
                    subjectToKeep = setdiff(1:ntot,trainix);
                end

                FA_keep_tmp{jj} = ...
                    out_up.FA_indi(:,:,:,subjectToKeep);
            end

            % Store the CP maps only once
            if ntrain == 1 && jj == 1
                FA_CP_tmp{jj} = out_up.FA_CP_indi;
            end
        end

        R.nrmse_US{ntrain} = nrmse_US_tmp;

        %% Store and check the CP reference

        if isempty(R.nrmse_CP)
            R.nrmse_CP = nrmse_CP_tmp(1,:).';
        end

        maxCPdifference = max( ...
            abs(nrmse_CP_tmp - R.nrmse_CP.'), ...
            [], 'all');

        tolerance = 1e-10 * max(1,max(abs(R.nrmse_CP)));

        if maxCPdifference > tolerance
            warning(['CP NRMSE differs between optimization runs. ' ...
                     'Check run_UP_experiment for run-dependent behavior.']);
        end

        %% Extract tailored results from ntrain = 1

        if ntrain == 1

            sampleFA = FA_keep_tmp{1};
            mapSize = size(sampleFA);
            mapSize(end+1:3) = 1;

            R.FA_TP   = nan(mapSize(1),mapSize(2),mapSize(3),ntot);
            R.FA_LOSO = nan(mapSize(1),mapSize(2),mapSize(3),ntot);

            for jj = 1:nsets

                ss = trainsets(jj,1);

                R.nrmse_TP(ss) = nrmse_US_tmp(jj,ss);
                R.FA_TP(:,:,:,ss) = FA_keep_tmp{jj};
            end

            R.FA_CP = FA_CP_tmp{1};
        end

        %% Extract LOSO results from ntrain = ntot-1

        if ntrain == ntot-1

            for jj = 1:nsets

                heldOutSubject = setdiff( ...
                    1:ntot, trainsets(jj,:));

                R.nrmse_LOSO(heldOutSubject) = ...
                    nrmse_US_tmp(jj,heldOutSubject);

                R.FA_LOSO(:,:,:,heldOutSubject) = ...
                    FA_keep_tmp{jj};
            end
        end

        fprintf( ...
            'Finished ntrain = %d: %d combinations\n', ...
            ntrain, nsets);
    end

    save(resultsFile,'R','-v7.3');

else
    load(resultsFile,'R');
end

%% Basic checks

assert(all(isfinite(R.nrmse_CP)), ...
    'Missing CP NRMSE values.');

assert(all(isfinite(R.nrmse_LOSO)), ...
    'Missing LOSO NRMSE values.');

assert(all(isfinite(R.nrmse_TP)), ...
    'Missing tailored NRMSE values.');

%% Prepare convergence results
%
% Convergence is shown normalized subject-wise to CP:
%
%       NRMSE_US(subject) / NRMSE_CP(subject)
%
% Each plotted point represents the mean over the training or test
% subjects belonging to one particular training combination.

nLevels = ntot-1;

trainSplitMean = cell(nLevels,1);
testSplitMean  = cell(nLevels,1);

qTrain = nan(nLevels,3);
qTest  = nan(nLevels,3);

fractionTestBetterThanCP = nan(nLevels,1);

for ntrain = 1:nLevels

    trainsets = R.trainsets{ntrain};

    % Subject-wise CP normalization
    normalizedNRMSE = ...
        R.nrmse_US{ntrain} ./ R.nrmse_CP.';

    nsets = size(trainsets,1);

    trainSplitMean{ntrain} = nan(nsets,1);
    testSplitMean{ntrain}  = nan(nsets,1);

    allTestValues = [];

    for jj = 1:nsets

        trainix = trainsets(jj,:);
        testix  = setdiff(1:ntot,trainix);

        trainSplitMean{ntrain}(jj) = mean( ...
            normalizedNRMSE(jj,trainix), ...
            'omitnan');

        testSplitMean{ntrain}(jj) = mean( ...
            normalizedNRMSE(jj,testix), ...
            'omitnan');

        allTestValues = [ ...
            allTestValues; ...
            normalizedNRMSE(jj,testix).']; %#ok<AGROW>
    end

    qTrain(ntrain,:) = prctile( ...
        trainSplitMean{ntrain}, [25 50 75]);

    qTest(ntrain,:) = prctile( ...
        testSplitMean{ntrain}, [25 50 75]);

    fractionTestBetterThanCP(ntrain) = ...
        mean(allTestValues < 1,'omitnan');
end

%% Automatically select representative subjects for the maps

relativeLOSO = R.nrmse_LOSO ./ R.nrmse_CP;

if isempty(mapSubjects)

    medianRelativePerformance = median(relativeLOSO);

    [~,medianSubject] = min( ...
        abs(relativeLOSO - medianRelativePerformance));

    [~,worstSubject] = max(relativeLOSO);

    mapSubjects = unique( ...
        [medianSubject worstSubject], ...
        'stable');
end

nMapSubjects = numel(mapSubjects);

%% Create main figure

figure( ...
    'Color','w', ...
    'OuterPosition',[92          -100        1366        1460]);

p = panel();
p.pack('v', 2);
p(1).pack('v', 2);
p(2).pack(3, nMapSubjects);

plotColors = lines(3);

% Panel A: convergence

axA = p(1,1).select();
hold on;

xTrain = [];
yTrain = [];

xTest = [];
yTest = [];

for ntrain = 1:nLevels

    xTrain = [ ...
        xTrain; ...
        repmat( ...
            ntrain-0.12, ...
            numel(trainSplitMean{ntrain}),1)]; %#ok<AGROW>

    yTrain = [ ...
        yTrain; ...
        trainSplitMean{ntrain}(:)]; %#ok<AGROW>

    xTest = [ ...
        xTest; ...
        repmat( ...
            ntrain+0.12, ...
            numel(testSplitMean{ntrain}),1)]; %#ok<AGROW>

    yTest = [ ...
        yTest; ...
        testSplitMean{ntrain}(:)]; %#ok<AGROW>
end

dotalpha = 0.25;
hTrainPoints = swarmchart( ...
    axA, xTrain, yTrain, ...
    10, plotColors(1,:), ...
    'filled', ...
    'MarkerFaceAlpha',dotalpha, ...
    'MarkerEdgeAlpha',dotalpha, ...
    'XJitterWidth',0.14);

hTestPoints = swarmchart( ...
    axA, xTest, yTest, ...
    10, plotColors(2,:), ...
    'filled', ...
    'MarkerFaceAlpha',dotalpha, ...
    'MarkerEdgeAlpha',dotalpha, ...
    'XJitterWidth',0.14);

hTrainSummary = errorbar( ...
    axA, ...
    (1:nLevels)-0.12, ...
    qTrain(:,2), ...
    qTrain(:,2)-qTrain(:,1), ...
    qTrain(:,3)-qTrain(:,2), ...
    'd-', ...
    'Color',plotColors(1,:), ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',plotColors(1,:), ...
    'MarkerSize',8, ...
    'LineWidth',2, ...
    'CapSize',8);

hTestSummary = errorbar( ...
    axA, ...
    (1:nLevels)+0.12, ...
    qTest(:,2), ...
    qTest(:,2)-qTest(:,1), ...
    qTest(:,3)-qTest(:,2), ...
    'd-', ...
    'Color',plotColors(2,:), ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',plotColors(2,:), ...
    'MarkerSize',8, ...
    'LineWidth',2, ...
    'CapSize',8);

uistack([hTrainSummary hTestSummary],'top');

xlabel(axA,'Number of training subjects');
ylabel(axA,'NRMSE / NRMSE_{CP}');

title(axA,'Convergence across training-set sizes');

xticks(axA,1:nLevels);
xlim(axA,[0.5 nLevels+0.5]);

grid(axA,'on');
box(axA,'on');

legend( ...
    axA, ...
    [hTrainPoints hTestPoints ...
     hTrainSummary hTestSummary], ...
    {'Training subsets', ...
     'Test subsets', ...
     'Training median and IQR', ...
     'Test median and IQR'}, ...
    'Location','best');

text( ...
    axA,-0.10,1.06,'A', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',13);

% Panel B: absolute LOSO comparison

axB = p(1,2).select();

barData = 100 .* [ ...
    R.nrmse_CP, ...
    R.nrmse_LOSO, ...
    R.nrmse_TP];

bar(axB,barData,'grouped');

xlabel(axB,'Held-out subject');
subjectTickLabels = {};
for jj = 1:ntot
    subjectTickLabels{jj} = sprintf('S%d', jj);
end
xticklabels(subjectTickLabels);
ylabel(axB,'NRMSE [%]');

title(axB,'Leave-one-subject-out performance');

xticks(axB,1:ntot);
xlim(axB,[0.4 ntot+0.6]);

ylim(axB,[0 1.12*max(barData,[],'all')]);

grid(axB,'on');
box(axB,'on');

legend( ...
    axB, ...
    {'CP', ...
     'Universal (LOSO)', ...
     'Tailored'}, ...
    'Location','best');

text( ...
    axB,-0.10,1.06,'B', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',13);

% Panel C: representative flip-angle maps

mapAxes = gobjects(3,nMapSubjects);

methodNames = { ...
    'CP', ...
    'Universal (LOSO)', ...
    'Tailored'};

for rr = 1:nMapSubjects

    ss = mapSubjects(rr);

    maps = { ...
        getSubjectSlice(R.FA_CP,ss), ...
        getSubjectSlice(R.FA_LOSO,ss), ...
        getSubjectSlice(R.FA_TP,ss)};

    nrmseValues = 100 .* [ ...
        R.nrmse_CP(ss), ...
        R.nrmse_LOSO(ss), ...
        R.nrmse_TP(ss)];

    targetMask = getSubjectSlice(mask,ss) > 0;

    for mm = 1:3 % CP, UP, TP

        ax = p(2,mm,rr).select();
        mapAxes(mm,rr) = ax;

        displayMap = abs(maps{mm});

        % Black background outside the target mask
        displayMap(~targetMask) = NaN;

        hImg = imagesc(ax,displayMap);
        hImg.AlphaData = ~isnan(displayMap);

        axis(ax,'image');
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);

        set(ax,'Color','k');

        clim(ax,mapCLim);

        hold(ax,'on');

        title(ax, sprintf('%.2f', nrmseValues(mm)./nrmseValues(1)), 'Color','red');

        if rr == 1
            set_img_ylabel(methodNames{mm});
        end

        if mm == 1

            if rr == 1
                text( ...
                    ax,-0.10,1.08,'C', ...
                    'Units','normalized', ...
                    'FontWeight','bold', ...
                    'FontSize',13);
            end
        end
    end
end
all_images_black_background;
colormap(gcf,turbo);

% Shared map colorbar at the right side of the layout
cb = colorbar(mapAxes(end,end));
cb.Label.String = 'Flip angle [degrees]';

p.fontsize=14;
p.margintop=10;
p.marginleft = 20;
p(1).marginbottom=30;
p(2).ch.margin = 0;
font_to_allAxes;

%% Display numerical summary

summaryTable = table( ...
    (1:ntot).', ...
    100.*R.nrmse_CP, ...
    100.*R.nrmse_LOSO, ...
    100.*R.nrmse_TP, ...
    relativeLOSO, ...
    'VariableNames',{ ...
        'Subject', ...
        'CP_NRMSE_percent', ...
        'LOSO_NRMSE_percent', ...
        'Tailored_NRMSE_percent', ...
        'LOSO_relative_to_CP'});

disp(summaryTable);

convergenceTable = table( ...
    (1:nLevels).', ...
    qTrain(:,2), ...
    qTrain(:,1), ...
    qTrain(:,3), ...
    qTest(:,2), ...
    qTest(:,1), ...
    qTest(:,3), ...
    fractionTestBetterThanCP, ...
    'VariableNames',{ ...
        'NTrainingSubjects', ...
        'TrainMedian', ...
        'TrainQ25', ...
        'TrainQ75', ...
        'TestMedian', ...
        'TestQ25', ...
        'TestQ75', ...
        'FractionTestBetterThanCP'});

disp(convergenceTable);

%% Local helper function

function img = getSubjectSlice(X,subjectIndex)
%GETSUBJECTSLICE Extract and orient one subject's sagittal map for display.
%   IMG = GETSUBJECTSLICE(X, SUBJECTINDEX) extracts SUBJECTINDEX from a 3-D
%   or 4-D stack X and applies the transpose/rotation used by the figure.
%   IMG has the same values and physical units as X.

    if ndims(X) >= 4 && size(X,4) >= subjectIndex

        img = X(:,:,1,subjectIndex);

    elseif ndims(X) == 3 && ...
            size(X,3) > 1 && ...
            size(X,3) >= subjectIndex

        img = X(:,:,subjectIndex);

    else

        img = X(:,:,1);
    end

    img = rot90(flip(img.',1),-1);
end
