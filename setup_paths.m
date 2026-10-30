function paths = setup_paths(tdmRoot)
%SETUP_PATHS Configure the minimal paths required by the figure scripts.

repoRoot = fileparts(mfilename('fullpath'));

if nargin < 1 || isempty(tdmRoot)
    tdmRoot = fullfile(repoRoot, '..', 'tdm_ptx');
end

[tdmExists, tdmInfo] = fileattrib(tdmRoot);
if tdmExists
    tdmRoot = tdmInfo.Name;
end

paths = struct(...
    'root', repoRoot, ...
    'data', fullfile(repoRoot, 'data'), ...
    'results', fullfile(repoRoot, 'results'), ...
    'tdm_ptx', tdmRoot);

localUtils = fullfile(repoRoot, 'utils');
panelDir = fullfile(repoRoot, 'third_party', 'panel-2.12');
paretoDir = fullfile(repoRoot, 'third_party', 'find_pareto_frontier');
tdmUtils = fullfile(tdmRoot, 'utils');
tdmLossFunctions = fullfile(tdmRoot, 'loss_functions');

requiredFiles = {...
    fullfile(localUtils, 'eval_pulse.m'), ...
    fullfile(panelDir, 'panel.m'), ...
    fullfile(paretoDir, 'find_pareto_frontier.m'), ...
    fullfile(tdmUtils, 'createAmat.m'), ...
    fullfile(tdmLossFunctions, 'constraint_SAR.m')};

missing = requiredFiles(~cellfun(@isfile, requiredFiles));
if ~isempty(missing)
    error('headNeckRfShimming:MissingDependency', ...
        'Missing required dependency: %s', missing{1});
end

% Add public dependencies first and local implementations last so that the
% intentionally compatible local eval_pulse wrapper takes precedence.
addpath(tdmLossFunctions, '-begin');
addpath(tdmUtils, '-begin');
addpath(paretoDir, '-begin');
addpath(panelDir, '-begin');
addpath(localUtils, '-begin');

assertResolvedFrom('eval_pulse', fullfile(localUtils, 'eval_pulse.m'));
assertResolvedFrom('calcLcurve', fullfile(localUtils, 'calcLcurve.m'));
assertResolvedFrom('panel', fullfile(panelDir, 'panel.m'));
assertResolvedFrom('find_pareto_frontier', ...
    fullfile(paretoDir, 'find_pareto_frontier.m'));
assertResolvedFrom('createAmat', fullfile(tdmUtils, 'createAmat.m'));
assertResolvedFrom('calcSarForPulse', ...
    fullfile(tdmUtils, 'calcSarForPulse.m'));
assertResolvedFrom('constraint_SAR', ...
    fullfile(tdmLossFunctions, 'constraint_SAR.m'));
end

function assertResolvedFrom(functionName, expectedFile)
resolvedFile = which(functionName);
if ~strcmpi(resolvedFile, expectedFile)
    error('headNeckRfShimming:PathConflict', ...
        '%s resolves to "%s" instead of "%s".', ...
        functionName, resolvedFile, expectedFile);
end
end
