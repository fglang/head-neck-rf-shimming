function all_images_black_background(figh, bgcolor)
%ALL_IMAGES_BLACK_BACKGROUND Set image-axis backgrounds to a uniform color.
%   ALL_IMAGES_BLACK_BACKGROUND(FIGH, BGCOLOR) updates axes containing image
%   objects in figure FIGH. BGCOLOR is a MATLAB color specification and
%   defaults to black. The function has no outputs.

    if nargin < 2
        bgcolor = 'black';
    end

    if nargin < 1
        figh = gcf;
    end

    allAxesInFigure = findall(figh, 'type', 'axes');

    for jj = 1:length(allAxesInFigure)
        if any(arrayfun(@(x) isa(x, 'matlab.graphics.primitive.Image'), ...
                allAxesInFigure(jj).Children))
            set(allAxesInFigure(jj), 'color', bgcolor);
            set(allAxesInFigure(jj), 'XColor', bgcolor);
            set(allAxesInFigure(jj), 'YColor', bgcolor);
        end
    end

    set(gcf, 'color', 'w');
end
