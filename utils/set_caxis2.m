function cax = set_caxis2(im, t)
%SET_CAXIS2 Compute display limits from lower and upper image quantiles.
%   CAX = SET_CAXIS2(IM, T) ignores NaN and Inf values and returns the
%   T and 1-T quantiles. T is a fraction and defaults to 1. CAX has the
%   same units as IM.

    if nargin < 2
        t = 1;
    end

    values = im(:);
    values(isnan(values)) = [];
    values(isinf(values)) = [];

    lower = quantile(values, t);
    upper = quantile(values, 1-t);
    
    cax = [lower, upper];
end
