function default_plot_settings()
%DEFAULT_PLOT_SETTINGS Apply shared formatting to the current figure.
%   Sets a white figure background, 14-point text, 1.5-point plot lines,
%   and black backgrounds for axes containing images. There are no inputs
%   or outputs.

set(gcf, 'color', 'w');
fontsize(gcf, 14, "points");
lineObjects = findobj(gcf, 'Type', 'Line');
for ii = 1:numel(lineObjects)
    lineObjects(ii).LineWidth = 1.5;
end
all_images_black_background(gcf);
