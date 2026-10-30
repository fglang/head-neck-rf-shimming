function out2d = tightzoom(inp2d, indicesOut)
%TIGHTZOOM Crop a 2-D array to the bounding box of its finite-valued region.
%   OUT2D = TIGHTZOOM(INP2D) removes all-NaN borders. With INDICESOUT = 1,
%   OUT2D is [[rowStart,rowEnd]; [columnStart,columnEnd]] instead of the
%   cropped array. The cropped values retain the units of INP2D.

    if nargin < 2
        indicesOut = 0;
    end
    
    inp2d = squeeze(inp2d);

    xmask = any(~isnan(inp2d),1);
    x1 = find(xmask, 1, 'first');
    x2 = find(xmask, 1, 'last');
    
    ymask = any(~isnan(inp2d),2);
    y1 = find(ymask, 1, 'first');
    y2 = find(ymask, 1, 'last');
    
    if indicesOut == 1
        out2d = [[y1, y2]; [x1, x2]];
    else
        out2d = inp2d(y1:y2, x1:x2);
    end
end
