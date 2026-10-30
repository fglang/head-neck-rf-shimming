function yl = set_img_ylabel(text, fontsize)
%SET_IMG_YLABEL Add a centered row label to an image axis.
%   YL = SET_IMG_YLABEL(TEXT, FONTSIZE) enables the current axis, adds TEXT
%   as its y-label, hides tick marks, and returns the label handle YL.
%   FONTSIZE is in points and defaults to 12.

    if nargin < 2
        fontsize = 12;
    end
    axis on;
    yl = ylabel(text, 'FontSize', fontsize, 'HorizontalAlignment', 'center');
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
end
