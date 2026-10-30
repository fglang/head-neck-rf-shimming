function [xcs, ycs] = convhull_sorted(x, y)
%CONVHULL_SORTED Return convex-hull points ordered by ascending y-value.
%   [XCS, YCS] = CONVHULL_SORTED(X, Y) computes the planar convex hull of
%   paired coordinates X and Y, then sorts its vertices by Y. Output units
%   match the corresponding input units.

k = convhull(x, y);
xc = x(k);
yc = y(k);

[ycs, sortix] = sort(yc);
xcs = xc(sortix);
end
