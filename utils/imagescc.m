function [h, tightZoomInd] = imagescc(im, showstats, tightzoomFlag, ROIs, showROILabels, ROIlabels)
%IMAGESCC Display an image with optional cropping, ROI contours, and stats.
%   [H, TIGHTZOOMIND] = IMAGESCC(IM, SHOWSTATS, TIGHTZOOMFLAG, ROIS,
%   SHOWROILABELS, ROILABELS) rotates IM into the publication orientation,
%   makes NaN pixels transparent, and optionally crops to finite values.
%   ROIS may be a binary mask or a cell array of structs with ROI_def fields.
%   H is the image handle; TIGHTZOOMIND contains crop row/column bounds.
%   Image values and color limits retain the units of IM.

if nargin < 6
    ROIlabels = {};
end

if nargin < 5
    showROILabels = 0;
end

if nargin < 4
    ROIs = {};
end

if ~isempty(ROIs) && ~isstruct(ROIs) % Convert one binary mask to ROI format.
    ROIs(isnan(ROIs)) = 0;
    ROIs = logical(ROIs);
    ROIs_new{1}.ROI_def = ROIs;
    ROIs = ROIs_new;
end

if nargin < 3
    tightzoomFlag = 1;
end

if nargin < 2
    showstats = 0;
end

if tightzoomFlag == 1
    tightZoomInd = tightzoom(im, 1);
    im = tightzoom(im);
end

h = imagesc(squeeze(rot90(im,1)), 'AlphaData', squeeze(~isnan(rot90(im,1))));

if numel(ROIs) > 0
    hold on;

    % Draw contours before labels so the contour does not cover label text.
    for ii = 1:numel(ROIs)
        if tightzoomFlag ~= 1
            ROI{ii} = rot90(ROIs{ii}.ROI_def,1);
        else
            ROI{ii} = rot90(ROIs{ii}.ROI_def(tightZoomInd(1,1):tightZoomInd(1,2), tightZoomInd(2,1):tightZoomInd(2,2)),1);
        end
        contour(ROI{ii}, 1, 'k-', 'LineWidth', 2);
    end

    for ii = 1:numel(ROIs)
        [yy, xx] = find(ROI{ii});
        if showROILabels == 1
            t_h = text(fix(max(xx)), fix(min(yy)), '# ');
            if isempty(ROIlabels)
                set(t_h, 'String', sprintf('#: %d ', ii), ...
                    'BackgroundColor', [1 1 1]);
            else
                set(t_h, 'String', sprintf('%s', ROIlabels{ii}), ...
                    'BackgroundColor', [1 1 1]);
            end
        end
    end
    set(gca,'xtick',[],'ytick',[]);
    clear yy xx
end

colorbar;
try
    caxis(set_caxis2(im,.01));
end

set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

set(gca,'ydir','reverse');

axis image;

if showstats == 1
    values = im(:);
    values(isnan(values)) = [];
    values(isinf(values)) = [];
    MEAN = mean(values);
    STD = std(values);
    hh = xlabel(sprintf('\\mu\\pm\\sigma=%.2g\\pm%.2g', MEAN, STD));
    hh.Interpreter = 'tex';
end
end
